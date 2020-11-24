// 22-08-2017

/////////////////////////
// 1. external includes //
//////////////////////////

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <omp.h>
#include "nlopt-2.4.2-dll64/nlopt.h"
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
#include "includes\sol_struct.cpp"  // define sol_struct + setup/destroy functions
#include "includes\utilities.cpp"   // transformation and utility functions

// c. solve
#include "includes\last_period.cpp"     // last period
#include "includes\interpolate.cpp"     // setup and interpolant gatewys
#include "includes\compute_WQ.cpp"      // compute W and Q
#include "includes\G2EGM_1D.cpp"        // G2EGM upper envelope
#include "includes\model.cpp"           // model functions
#include "includes\vfi.cpp"             // value function iteration 
#include "includes\model_vfi_pure.cpp"  // model functions for pure vfi
#include "includes\vfi_pure.cpp"        // pure value function iteration


////////////////
// 4. gateway //
////////////////

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
     
        logs::solve(-1,"Solving model...\n"); // reset log file
        
        // for timing
        double time_compute_WQ, time_G2EGM_1D, time_vfi;
        HighResTimer timer, timer_all, timer_sec;
        timer.StartTimer();


    //////////////////
    // 5. setup par //
    //////////////////
    
    par_struct* par = new par_struct;
    par::setup(par,plhs,prhs,1);
    if(par->vfi_pure == 0){
        logs::solve(0," using NEGM\n");
    } else {
        logs::solve(0," using VFI\n");
    }


    /////////////////
    // 6. parallel //
    /////////////////

    #pragma omp parallel num_threads(MAXTHREADS)
    {

        // a. setup workers
        auto sol = new sol_struct;
        sol::setup(par,sol);


        // b. last period

            #pragma omp master
            logs::solve(0,"last period\n");

        last_period(par->TR-1, par, sol);

        // c. time loop
        for(int t = par->TR-2; t >= par->tmin; t--){ // start from next to last

                #pragma omp master
                logs::solve(0,"t = %d\n", t);

            sol->t = t;

            //////////
            // NEGM //
            //////////

            if(par->vfi_pure == 0){ // NEGM

                    #pragma omp master 
                    timer_all.StartTimer();

                // i. create interpolants
                interpolate::setup(t,par,sol);
                    
                // ii. precomputations
                    #pragma omp master
                    timer_sec.StartTimer();
                compute_WQ::all(t,par,sol);
                    #pragma omp master 
                    time_compute_WQ = timer_sec.StopTimer();

                // iii. G2EGM
                    #pragma omp master
                    timer_sec.StartTimer();          
                G2EGM_1D::all(t,par,sol);
                    #pragma omp master 
                    time_G2EGM_1D = timer_sec.StopTimer();

                // iv. solve rest by vfi
                    #pragma omp master
                    timer_sec.StartTimer();        
                vfi::solve(t,par,sol);
                    #pragma omp master 
                    time_vfi = timer_sec.StopTimer();
                
                    #pragma omp master 
                    {
                    double time_all = timer_all.StopTimer();

                        logs::solve(1,"time        = %3.2f secs\n", time_all);     
                        logs::solve(1," compute_WQ = %3.2f secs\n", time_compute_WQ); 
                        logs::solve(1," G2EGM_1D   = %3.2f secs\n", time_G2EGM_1D);                     
                        logs::solve(1," vfi        = %3.2f secs\n", time_vfi);                                                            
                    
                    }
                    #pragma omp barrier


            //////////////
            // PURE VFI //
            //////////////
            
            } else { // pure VFI

                    #pragma omp master 
                    timer_all.StartTimer();

                // i. create interpolants
                interpolate::setup_next_period(t,par,sol,1);

                // ii. solve rest by vfi        
                vfi_pure::solve(t,par,sol);

                    #pragma omp master 
                    {
                    double time_all = timer_all.StopTimer();

                        logs::solve(1,"time = %3.2f secs\n", time_all);                                                                
                    
                    }
                    #pragma omp barrier

            }
            
            //////////////////////////
            // DESTROY INTERPOLANTS //            
            //////////////////////////

            if(par->vfi_pure == 0){ // NEGM
                interpolate::destroy(t,par,sol);                
            } else {                
                interpolate::destroy_next_period(t,par,sol);                                
            }
            

        } // t

        // c. clean up workers
        sol::destroy(sol);
        delete sol;

    } // parallel


    /////////////////
    // 7. clean up //
    /////////////////
    
    par::destroy(par,1);
    delete par;

        double time = timer.StopTimer();
        logs::solve(0,"Time: %5.2f secs\n",time);
        logs::solve(0,"Done.\n");

        // clean assertions file
        FILE* log_file = fopen("log_assert.txt","w");
        fprintf(log_file,"\n");
        fclose(log_file);
        
} // mex gateway