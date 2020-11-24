//////////////////////////
// 1. external includes //
//////////////////////////

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <omp.h>
#include "mex.h"
#include "matrix.h"


//////////////////////////
// 2. define statements //
//////////////////////////

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define BOUND(X,A,B) MIN(MAX(X,A),B)
#define THREADS MIN(MAXTHREADS, omp_get_max_threads()-1)


//////////////////////////
// 3. internal includes //
//////////////////////////

// a. generic
#include "includes\HighResTimer_class.hpp" // class for timing
#include "includes\assert.cpp"             // assert() function
#include "includes\logs.cpp"               // log:: functions
#include "includes\linear_interp.cpp"      // linear_interp:: functions
#include "includes\index.cpp"              // index:: functions
#include "includes\mymex.cpp"              // functions to interact with mex

// b. basic
#include "includes\par_struct.cpp"  // define par_struct + setup/destroy functions
#include "includes\utilities.cpp"   // transformation and utility functions

////////////////
// 4. gateway //
////////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    logs::simulate(-1,"Starting simulation...\n");


    //////////////////
    // 5. setup par //
    //////////////////
    
    par_struct* par = new par_struct;
    par::setup(par,plhs,prhs,2);

        logs::simulate(0,"Setup completed.\n");

        // easy access to outputs
        auto d = par->d;
        auto z = par->z;
        
        auto M = par->M;
        auto N = par->N;
        auto P = par->P;
        auto X = par->X;                
        auto Y = par->Y;    
        auto H = par->H;
        auto C = par->C;
        auto B = par->B; 
        auto A = par->A; 

        auto MPC  = par->MPC; 
        auto simV = par->simV;


    /////////////////
    // 6. parallel //
    /////////////////
        
    #pragma omp parallel num_threads(THREADS)
    {

        // allocate memory for interpolants 
        auto interp = new linear_interp::set_struct*[par->Nd*par->Nz*par->NP];             

    //////////////////
    // 7. time loop //
    //////////////////

    for(int t = 0; t < par->Tendsim; t++){ // forward through time

        #pragma omp master
        logs::simulate(0,"t = %d\n",t);
        

    ////////////////////
    // 8. intepolants //
    ////////////////////

        logs::simulate(1,"interpolants setup\n",t);

    for(int d_now = 0; d_now < par->Nd; d_now++){
    for(int z_now = 0; z_now < par->Nz; z_now++){
    for(int i_P   = 0; i_P   < par->NP; i_P++){
         
        // b. grids
    
            // i. dimensions
            int dimx;
            int Nx[2];   
            double *x[2];                     
            if(z_now == 1){ // adjustment, X
                
                dimx  = 1; 
                Nx[0] = par->NX;                     
                x[0]  = &par->grid_X[index::d3(t, i_P, 0, par->NP, par->NX)];

            } else if(par->NN == 1){ // no-adjustment, M

                dimx  = 1; 
                Nx[0] = par->NM;                     
                x[0]  = &par->grid_M[index::d3(t, i_P, 0, par->NP, par->NM)];

            } else {  // no-adjustment, M and N

                dimx  = 2;                
                Nx[0] = par->NM;
                Nx[1] = par->NN;
                x[0]  = &par->grid_M[index::d3(t, i_P, 0, par->NP, par->NM)];
                x[1]  = &par->grid_N[index::d3(t, i_P, 0, par->NP, par->NN)];

            }

        // c. values
        int i_cell = index::d3(t,d_now,z_now,par->Nd,par->Nz);
        int i_grid;
        if(z_now == 0){
            i_grid = index::d3(i_P,0,0,par->NN,par->NM);
        } else {
            i_grid = index::d2(i_P,0,par->NX);        
        }

        int dimy = 2;
        if(z_now == 1){
            dimy = 3;
        }

        double *y[3];
        y[0] = &par->invV[i_cell][i_grid];
        y[1] = &par->C_ast[i_cell][i_grid];
        y[2] = &par->B_ast[i_cell][i_grid];                

        // d. create
        int i_interp = index::d3(d_now, z_now, i_P, par->Nz, par->NP);
        interp[i_interp] = new linear_interp::set_struct;
        linear_interp::create(interp[i_interp], dimx, dimy, Nx, x, y, 1);

    } } }


    /////////////////////////
    // 9. individuals loop //
    /////////////////////////

        logs::simulate(1,"individuals loop\n",t);

    #pragma omp for
    for(int i = 0; i < par->Nsim; i++){

            logs::simulate(2,"i = %d\n",i);

        // a. index
        int index      = t*par->Nsim + i;
        int index_plus = (t+1)*par->Nsim + i;

            logs::simulate(2,"index = %d, index_plus = %d\n", index, index_plus);

            if(t < par->T0sim){
                M[index] = NAN;
                N[index] = NAN;
                P[index] = NAN;
                Y[index] = NAN;                                
                X[index] = NAN;
                H[index] = NAN;                
                C[index] = NAN;                
                B[index] = NAN;
                A[index] = NAN;
                MPC[index] = NAN;
                continue;
            }

        // b. initial values
        if(t == 0 && par->sim_IRF == 0){

            // i. inheritance status
            d[i] = 0;

            // ii. income
            P[i] = par->psi_sim[i]*par->P_ini_dist[i];
            Y[i] = P[i]*par->xi_sim[i];

            // iii. states        
            if(par->NN > 1){
                M[i] = Y[i];
                N[i] = par->N_ini_dist[i];
            } else {
                M[i] = par->N_ini_dist[i] + Y[i];
                N[i] = 0.0;
            }
            X[i] = M[i] + N[i] - par->adjcost + par->omega[t]*P[i];

                logs::simulate(2,"initial values\n",i);

        } else if(t == par->T0sim && par->sim_IRF == 1){

            M[index] = NAN;
            N[index] = NAN;
            Y[index] = NAN;                                
            X[index] = NAN;               
            C[index] = NAN;                
            MPC[index] = NAN;

            // i. base period
            d[index] = 1;
            P[index] = par->P_IRF[i];
            A[index] = par->A_IRF[i];
            B[index] = par->B_IRF[i];                        

            // ii. states for next period
            double Rnow  = par->R;
            if(A[index]< 0){
                Rnow = par->Rneg;
            }
            P[index_plus] = par->G[t]*P[index]*par->psi_sim[index_plus];
            Y[index_plus] = P[index_plus]*par->xi_sim[index_plus];
            M[index_plus] = Rnow*A[index] + Y[index_plus] + H[index];
            N[index_plus] = par->Rb*B[index];
            X[index_plus] = M[index_plus] + N[index_plus] - par->adjcost + P[index_plus]*par->omega[t+1];

            continue;

        }

        // c. forget inheritance? 
        int d_now = d[index];
        if(par->forget_inh == 1){
            d_now = 1; // act as if inheritance has been received
        }


    /////////////////////////
    // 10. optimal choices //
    /////////////////////////

        logs::simulate(3,"d = %d, M = %g, N = %g, P = %g, X = %g\n",
            d[index],M[index],N[index],P[index],X[index]);   

    // i. adjustment cost shock
    int z_max;
    if(par->adjcost_pval[index] < par->adjcostinf){ 
        z_max = 1; 
    } else {
        z_max = par->Nz;
    }

    // ii. loop over adjustment choices
    double invV = -HUGE_VAL;
    for(int z_now = 0; z_now < z_max; z_now++){

            logs::simulate(4,"z = %d:\n",z_now);   

        // a. states
        double xi[2];
        if(z_now == 1){
            xi[0] = X[index];
        } else if(par->NN == 1){
            xi[0] = M[index];
        } else {        
            xi[0] = M[index];
            xi[1] = N[index];                
        } 

        // b. P indexes and relative difference
        auto P_vec       = &par->grid_P[t*par->NP];
        int i_P_left     = linear_interp::binary_search(0, par->NP, P_vec, P[index]);
        double P_reldiff = (P[index] - P_vec[i_P_left])/(P_vec[i_P_left+1]-P_vec[i_P_left]);
       
            logs::simulate(5,"P = %g, P_reldiff = %g, i_P_left = %d\n",P[index],P_reldiff,i_P_left);  
            
        // c. interpoaltion
        double results[2*3];
        double results_MPC[2*3];            
        for(int add = 0; add < 2; add++){
            int i_P      = i_P_left+add;
            int i_interp = index::d3(d_now, z_now, i_P, par->Nz, par->NP);
            linear_interp::evaluate(interp[i_interp],&results[add*3],xi,nullptr);
        }      

        xi[0] += par->MPC_add;
        for(int add = 0; add < 2; add++){
            int i_P      = i_P_left+add;
            int i_interp = index::d3(d_now, z_now, i_P, par->Nz, par->NP);
            linear_interp::evaluate(interp[i_interp],&results_MPC[add*3],xi,nullptr);
        }

            // print
            if(z_now == 0){ 
                logs::simulate(5," invV = [%g %g], C = [%g %g]\n"
                    ,results[0],results[3],results[1],results[4]);  
            } else {
                logs::simulate(5," invV = [%g %g], C = [%g %g], B = [%g %g]\n"
                    ,results[0],results[3],results[1],results[4],results[2],results[5]);  
            }

        // d. interpolation in P dimension
        double invV_now   = results[0] + P_reldiff * (results[3]-results[0]);
        double C_now      = results[1] + P_reldiff * (results[4]-results[1]);
        double C_MPC  = results_MPC[1] + P_reldiff * (results_MPC[4]-results_MPC[1]);        

        // e. implied choices and MPC
        double A_now, B_now;
        if(z_now == 1){
            B_now = results[2] + P_reldiff * (results[5]-results[2]);
            A_now = X[index] - C_now - B_now - P[index]*par->omega[t];         
        } else {
            if(t < par->TR-1){
                B_now = N[index];
                A_now = M[index] - C_now;          
            } else {
                B_now = 0;
                A_now = M[index] + N[index] - C_now;          
            }

        }

        logs::simulate(4," invV = %g, C = %g, B = %g, A = %g\n",invV_now,C_now,B_now,A_now);            

        // f. update
        if(invV_now > invV){

            invV     = invV_now;            
            z[index] = z_now;
            C[index] = C_now;
            B[index] = B_now;
            A[index] = A_now;            
            MPC[index] = (C_MPC - C_now) / par->MPC_add;

            simV[index] = utilities::inv_trans(invV,par);

        }

    }

        logs::simulate(3,"optimal z = %d\n", z[index]);
        logs::simulate(3,"optimal choices, C = %g, B = %g, A = %g\n", C[index], B[index], A[index]);
        
    /////////////////////
    // 11. next period //
    /////////////////////

    if (t < par->Tendsim-1) {

        // a. add wealth shock
        H[index] = par->A_shock[index];

        // ii. add inheritance
        if(par->Nd == 1){
            d[index_plus] = 0;
        } else {
            d[index_plus] = 1;
            if(d[index] == 0){ // chance of receiving
                if(t+1 == par->Tendsim-1 || par->inh_pval[index_plus] < par->inh_p[t+1]){ // receives
                    H[index] += par->inh_A*par->inh_fac[t];
                    d[index_plus] = 1;
                } else { // don't receives
                    d[index_plus] = 0;
                }
            }
        }

        // iii. interest rate
        double Rnow  = par->R;
        if(A[index] < 0){
            Rnow = par->Rneg;
        }

        // iv. next period variables
        P[index_plus] = par->G[t]*P[index]*par->psi_sim[index_plus];
        Y[index_plus] = P[index_plus]*par->xi_sim[index_plus];
        M[index_plus] = Rnow*A[index] + Y[index_plus] + H[index];
        N[index_plus] = par->Rb*B[index];
        X[index_plus] = M[index_plus] + N[index_plus] - par->adjcost + P[index_plus]*par->omega[t+1];

            logs::simulate(3,"next period\n");

    }


    } // N
    #pragma omp barrier

        // clean up workers
        for(int d_now = 0; d_now < par->Nd; d_now++){
        for(int z_now = 0; z_now < par->Nz; z_now++){
        for(int i_P   = 0; i_P   < par->NP; i_P++){
            int i_interp = index::d3(d_now, z_now, i_P, par->Nz, par->NP);
            linear_interp::destroy(interp[i_interp]);
        } } }

    } // T        
    } // parallel


    //////////////////////
    // 12. save to disc //
    //////////////////////

    if(par->save_to_disc > 0){
    
        logs::simulate(0,"Saving to disc...\n");

        FILE* sim_out = fopen("sim_output.txt","w");

        for(int i = 0; i < par->save_to_disc; i++){

            // initial values
            if(par->NN > 1) {
                fprintf(sim_out,"%d,%g,%g,%g,%g,%g,%g\n",24,
                    NAN,
                    0.0,
                    par->N_ini_dist[i]/par->R,
                    0.0,
                    0.0,
                    0.0);
            } else {
                fprintf(sim_out,"%d,%g,%g,%g,%g,%g,%g\n",24,
                        NAN,
                        par->N_ini_dist[i]/par->R,
                        0.0,
                        0.0,
                        0.0,
                        0.0);            
            }

        for(int t = 0; t < par->Tendsim; t += 1){ // forward through time

            int index = t*par->Nsim + i;
            fprintf(sim_out,"%d,%g,%g,%g,%g,%g,%g\n",(int)(25+(double)t),
                C[index]/P[index],
                (A[index]+H[index])/P[index],
                B[index]/P[index],
                H[index]/P[index],
                Y[index]/P[index],
                P[index]);

        } }
        fclose(sim_out);

    }


    //////////////////
    // 13. clean up //
    //////////////////

    par::destroy(par,2);
    delete par;

        logs::simulate(0,"Done.\n");

        // clean assertions file
        FILE* log_file = fopen("log_assert.txt","w");
        fprintf(log_file,"\n");
        fclose(log_file);     

} // simulate