// 22-08-2017

namespace compute_WQ {

    // forward declarations
    double weigth_d_func(sol_struct *sol);
    void calc(sol_struct *sol);

    void next_period_values(sol_struct *sol);
    void update(sol_struct *sol);

    void next_period_values_vfi_pure(sol_struct *sol);
    void update_vfi_pure(sol_struct *sol);

////////////////
// 1. profile //
////////////////

void all(int t, par_struct *par, sol_struct *sol){
	
    int level = 2;
        
        logs::solve(level,"compute_WQ::profile\n");

    // parallel
    #pragma omp for collapse(3)
    for(int d   = 0; d   < par->Nd; d++){
    for(int i_P = 0; i_P < par->NP; i_P++){
    for(int i_B = 0; i_B < par->NN; i_B++){

            logs::solve(level+1,"(d,i_P,i_N) = (%1d,%3d,%3d)\n",d,i_P,i_B);

        int z = 0; // z is not a post-decision state

    	// a. pack
    	sol->t = t;
    	sol->d = d;
    	sol->z = z;
    	sol->P = par->grid_P[index::d2(t,i_P,par->NP)];
        sol->B = par->grid_N[index::d3(t,i_P,i_B,par->NP,par->NN)];

        // b. post-decision vectors
        sol->A_vec = &par->grid_A[index::d3(t,i_P,0,par->NP,par->NA)];

        int i_cell = index::d3(0,d,0,par->Nd,1);           // 1 as z is not a post-decision state
        int i_grid = index::d3(i_P,i_B,0,par->NN,par->NA); // i_A = 0, first element

        sol->W_vec    = &par->W[i_cell][i_grid];
        sol->invW_vec = &par->invW[i_cell][i_grid];
        sol->Q_vec    = &par->Q[i_cell][i_grid];

        // c. post-decision vectors W and Q
        calc(sol);

    } } }

}


//////////////////////
// 2. post-decision //
//////////////////////

void calc(sol_struct *sol){

    int level = 3;

        logs::solve(level,"compute_WQ::do\n");    

    //////////////
    // 1. setup //
    //////////////

    // a. unpack
	auto par      = sol->par;
	auto t        = sol->t;
    auto W_vec    = sol->W_vec;
    auto invW_vec = sol->invW_vec;
    auto Q_vec    = sol->Q_vec;

	// b. initialize
    if(par->vfi_pure == 0){
        for(int i_A = 0; i_A < par->NA; i_A++){
           W_vec[i_A] = 0.0;
           Q_vec[i_A] = 0.0;
        }
    } else {
        sol->W = 0;
    }


    /////////////
    // 2. loop //
    /////////////
    
    auto W_shocks = sol->W_shocks;
    auto Q_shocks = sol->Q_shocks;
    for(int d_plus = 0; d_plus < par->Nd; d_plus++){
            
            logs::solve(level+1,"d_plus = %d\n",d_plus);

        // a. weight
        sol->d_plus = d_plus;
    	double weigth_d_plus = weigth_d_func(sol);
    	if(weigth_d_plus <= 1e-8){ continue;}

            logs::solve(level+1,"weight = %g\n",weigth_d_plus);

        // b. initialize
        if(par->vfi_pure == 0){
    	    for(int i_A = 0; i_A < par->NA; i_A++){
    	       W_shocks[i_A] = 0.0;
    	       Q_shocks[i_A] = 0.0;
    	    }
        } else {
            W_shocks[0] = 0.0; 
        }

        // c. loop over shocks 
    	int Nshocks  = par->Nshocks[t+1];
        sol->psi     = &par->psi[(t+1)*par->Nshocks_max];       
        sol->xi      = &par->xi[(t+1)*par->Nshocks_max];
        sol->weights = &par->weights[(t+1)*par->Nshocks_max];
    	for(int i_shock = 0; i_shock < Nshocks; i_shock++){

                logs::solve(level+2,"i_shock = %d\n",i_shock);

            sol->i_shock = i_shock;

            if(par->vfi_pure == 0){

                // i. next-period states
        		next_period_values(sol);

                // ii. interpolaion
        		interpolate::next_period(sol,sol->invV_plus_vec,sol->C_plus_vec);

                // iii. update
        		update(sol);

            } else {

                // i. next-period states
                next_period_values_vfi_pure(sol);

                // ii. interpolaion
                interpolate::next_period_vfi_pure(sol,sol->invV_plus_vec);

                // iii. update
                update_vfi_pure(sol);

            }    	
        }

        // d. accumulate 
        if(par->vfi_pure == 0){
        	for(int i_A = 0; i_A < par->NA; i_A++){
        		W_vec[i_A] += weigth_d_plus*W_shocks[i_A];
        		Q_vec[i_A] += weigth_d_plus*Q_shocks[i_A];
        	}
        } else {
            sol->W += weigth_d_plus*W_shocks[0];
        }

    }


    ////////////////////
    // 4. transform W //
    ////////////////////

    if(par->vfi_pure == 0){
        if(par->epstein_zin == 1){
            for(int i_A = 0; i_A < par->NA; i_A++){
                W_vec[i_A] = pow(W_vec[i_A],1.0/(1.0-par->rho));
                Q_vec[i_A] = Q_vec[i_A]*pow(W_vec[i_A],par->rho-par->sigma);
            }
        }
        for(int i_A = 0; i_A < par->NA; i_A++){
            invW_vec[i_A] = utilities::trans(W_vec[i_A],par);
        }
    } else {
        if(par->epstein_zin == 1){
            sol->W = pow(sol->W,1.0/(1.0-par->rho));
        }
    }
    
}

double weigth_d_func(sol_struct *sol){

    auto par    = sol->par;
    auto t      = sol->t;
    auto d      = sol->d;
    auto d_plus = sol->d_plus;    

    if(par->Nd == 1){
        return 1.0;
    }

    if(d == 0){ // some chance of receiving inheritance
        if(d_plus == 1){ // receive inheritance
            if(t == par->TR-2){
                return 1.0;
            } else {
                return par->inh_p[t];
            }
        } else { // don't receive inheritance
            if(t == par->TR-2){
                return 0.0;
            } else {
                return 1.0-par->inh_p[t];
            }
        }
    } else {
        if(d_plus == 1){ // parent stay death
            return 1.0;
        } else { // resurrection
            return 0.0;
        }
    }

}

void next_period_values(sol_struct *sol){

    int level = 5;

        logs::solve(level,"next_period_values\n");

    // a. unpack
    auto par = sol->par;
    auto t   = sol->t;

    double psi = sol->psi[sol->i_shock];
    double xi  = sol->xi[sol->i_shock];
        
        logs::solve(level+1,"psi = %g, xi = %g\n",psi,xi);

    // b. inheritance
    double inh_A = 0.0;
    if(sol->d == 0 && sol->d_plus == 1){ // receive inheritance
        inh_A = par->inh_A*par->inh_fac[t];        
    }

        logs::solve(level+1,"inh_A = %g\n",inh_A);

    // c. income
        
        double P_min = par->grid_P[(t+1)*par->NP];
        double P_max = par->grid_P[(t+1)*par->NP + par->NP-1];

    sol->P_plus = par->G[t]*sol->P*psi;
    sol->P_plus = BOUND(sol->P_plus,P_min,P_max);
    sol->Y_plus = sol->P_plus*xi;
    sol->inh_A  = inh_A;
    
    // d. illiquid assets
    sol->N_plus = par->Rb*sol->B;

        logs::solve(level+1,"P_plus = %g, Y_plus = %g, N_plus = %g\n",sol->P_plus,sol->Y_plus,sol->N_plus);

    // e. cash-on-hand
    for(int i_A = 0; i_A < par->NA; i_A++){

        // i. interest rate
        double Rnow = par->R;
        if(sol->A_vec[i_A] < 0){ Rnow = par->Rneg;}        

        // ii. cash-on-hand without adjustment
        sol->M_plus_vec[i_A] = Rnow*sol->A_vec[i_A] + sol->Y_plus + inh_A;

        // iii. cash-on-hand after adjustment
        sol->X_plus_vec[i_A] = sol->M_plus_vec[i_A] + sol->N_plus - par->adjcost + sol->P_plus*par->omega[t+1];

            logs::solve(level+2,"i_A = %d: M_plus = %g, X_plus = %g\n",i_A,sol->M_plus_vec[i_A],sol->X_plus_vec[i_A]);

    }
    
}

void next_period_values_vfi_pure(sol_struct *sol){

    int level = 5;

        logs::solve(level,"next_period_values\n");

    // a. unpack
    auto par = sol->par;
    auto t   = sol->t;

    double psi = sol->psi[sol->i_shock];
    double xi  = sol->xi[sol->i_shock];
        
        logs::solve(level+1,"psi = %g, xi = %g\n",psi,xi);

    // b. inheritance
    double inh_A = 0.0;
    if(sol->d == 0 && sol->d_plus == 1){ // receive inheritance
        inh_A = par->inh_A*par->inh_fac[t];        
    }

        logs::solve(level+1,"inh_A = %g\n",inh_A);

    // c. income

        double P_min = par->grid_P[(t+1)*par->NP];
        double P_max = par->grid_P[(t+1)*par->NP + par->NP-1];
    
    sol->P_plus = par->G[t]*sol->P*psi;
    sol->P_plus = BOUND(sol->P_plus,P_min,P_max);
    sol->Y_plus = sol->P_plus*xi;
    sol->inh_A  = inh_A;

        logs::solve(level+1,"P_plus = %g, Y_plus = %g\n",sol->P_plus,sol->Y_plus);
    
    // d. illiquid assets
    sol->N_plus = par->Rb*sol->B;

        logs::solve(level+1,"N_plus = %g\n",sol->N_plus);

    // e. interest rate
    double Rnow = par->R;
    if(sol->A < 0){ Rnow = par->Rneg;}        

    // f. cash-on-hand
    sol->M_plus = Rnow*sol->A + sol->Y_plus + inh_A;

    // g. without adjustment
    sol->X_plus = sol->M_plus + sol->N_plus - par->adjcost + sol->P_plus*par->omega[t+1];

        logs::solve(level+2,"M_plus = %g, X_plus = %g\n",sol->M_plus,sol->X_plus);

}

void update(sol_struct *sol){

    int level = 5;

        logs::solve(level,"update\n");

    auto par = sol->par;

    for(int i_A = 0; i_A < par->NA; i_A++){

        double weight = sol->weights[sol->i_shock];
        
        // a. optimal z_plus
        int z_plus_ast = -1;
        double invV_plus_ast = -HUGE_VAL;
        for(int z_plus = 0; z_plus < par->Nz; z_plus++){
            auto invV_plus_now = sol->invV_plus_vec[z_plus*par->NA + i_A];
            if(invV_plus_now > invV_plus_ast){
                invV_plus_ast = invV_plus_now;
                z_plus_ast    = z_plus;
            }
        }  

        // b. updates
        for(int z_plus = 0; z_plus < z_plus_ast+1; z_plus++){         
            
            // i. weight
            double weight_now;
            if(z_plus_ast == 0){
                weight_now = weight;
            } else if(z_plus == 0){
                weight_now = par->adjcostinf*weight;
            } else if(z_plus == 1){
                weight_now = (1.0-par->adjcostinf)*weight;
            }
            auto invV_plus = sol->invV_plus_vec[z_plus*par->NA + i_A];   
            auto C_plus = sol->C_plus_vec[z_plus*par->NA + i_A];          

                // without positive C_plus the algorithm might break down
                assert(1,invV_plus > 0,"invV_plus = %g, C_plus = %g, A = %g, B = %g",
                    invV_plus,C_plus,sol->A_vec[i_A],sol->B);

            // ii. post-decision value function
            double V_plus = NAN;
            if(par->epstein_zin == 0){ // CRRA
                sol->W_shocks[i_A] += weight_now*utilities::inv_trans(invV_plus,par);
            } else { // Epstein-Zin
                V_plus = utilities::inv_trans(invV_plus,par);
                sol->W_shocks[i_A] += weight_now*pow(V_plus,1.0-par->rho); 
            }
            
            // iii. marginal utility of cash
            double Rnow = par->R;
            if(sol->A_vec[i_A] < 0){ Rnow = par->Rneg;}    

                // without positive C_plus the algorithm might break down
                assert(1,C_plus > 0,"C_plus = %g, A = %g, B = %g",
                    C_plus,sol->A_vec[i_A],sol->B);

            if(par->epstein_zin == 0){ // CRRA
                sol->Q_shocks[i_A] += weight_now*par->beta*Rnow*utilities::marg_u(C_plus,par);            
            } else { // Epstein-Zin
                sol->Q_shocks[i_A] += weight_now*par->beta*Rnow*utilities::marg_u(C_plus,par)*pow(V_plus,par->sigma-par->rho); 
            }   
        } 
    }

}

void update_vfi_pure(sol_struct *sol){

    int level = 5;

        logs::solve(level,"update\n");

    auto par      = sol->par;
    double weight = sol->weights[sol->i_shock];

    // a. optimal z_plus
    int z_plus_ast = -1;
    double invV_plus_ast = -HUGE_VAL;
    for(int z_plus = 0; z_plus < par->Nz; z_plus++){
        auto invV_plus_now = sol->invV_plus_vec[z_plus];
        if(invV_plus_now > invV_plus_ast){
            invV_plus_ast = invV_plus_now;
            z_plus_ast    = z_plus;
        }
    }  

    // b. updates
    for(int z_plus = 0; z_plus < z_plus_ast+1; z_plus++){         
        
        // i. weight
        double weight_now;
        if(z_plus_ast == 0){
            weight_now = weight;
        } else if(z_plus == 0){
            weight_now = par->adjcostinf*weight;
        } else if(z_plus == 1){
            weight_now = (1-par->adjcostinf)*weight;
        }
        auto invV_plus = sol->invV_plus_vec[z_plus];   

            logs::solve(level,"invV_plus = %g, weight = %g\n",invV_plus, weight_now);

        // ii. post-decision value function
        double V_plus = NAN;
        if(par->epstein_zin == 0){ // CRRA
            sol->W_shocks[0] += weight_now*utilities::inv_trans(invV_plus,par);
        } else { // Epstein-Zin
            V_plus = utilities::inv_trans(invV_plus,par);
            sol->W_shocks[0] += weight_now*pow(V_plus,1.0-par->rho); 
        }
           
    } 

}

} // namespace