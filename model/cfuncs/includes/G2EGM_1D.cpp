namespace G2EGM_1D {

	// forward declaration
	void choice_and_value(sol_struct *sol);

////////////
// 1. all //
////////////

void all(int t, par_struct *par, sol_struct *sol){

    int level = 2;
        
        logs::solve(level,"G2EGM_1D::profile\n");

    // parallel
    #pragma omp for collapse(3)
    for(int d   = 0; d   < par->Nd; d++){
    for(int i_P = 0; i_P < par->NP; i_P++){
    for(int i_B = 0; i_B < par->NN; i_B++){

            logs::solve(level+1,"(d,i_P,i_B) = (%1d,%3d,%3d)\n",d,i_P,i_B);

    	int z = 0; // z is not a post-decision state

    	// a. pack
    	sol->t = t;
    	sol->d = d;
    	sol->z = z;
    	sol->P = par->grid_P[index::d2(t,i_P,par->NP)];
    	sol->B = par->grid_N[index::d3(t,i_P,i_B,par->NP,par->NN)];

        // b. indexes
        int i_cell    = index::d3(t,d,z,par->Nd,par->Nz);
        int i_pd_cell = index::d3(0,d,0,par->Nd,1);

        int i_states    = index::d3(i_P,i_B,0,par->NN,par->NM);
        int i_pd_states = index::d3(i_P,i_B,0,par->NN,par->NA);
	        
        // c. post-decision vectors
        sol->A_vec = &par->grid_A[index::d3(t,i_P,0,par->NP,par->NA)];

        sol->W_vec    = &par->W[i_pd_cell][i_pd_states];
        sol->invW_vec = &par->invW[i_pd_cell][i_pd_states];
        sol->C_vec    = &par->C_EGM[i_pd_cell][i_pd_states]; 
        sol->Q_vec    = &par->Q[i_pd_cell][i_pd_states];
        sol->M_vec    = &par->M_EGM[i_pd_cell][i_pd_states];

        // d. profile vectors
        sol->M_common  = &par->grid_M[index::d3(t,i_P,0,par->NP,par->NM)];

        sol->V_vec     = &par->V[i_cell][i_states];       
        sol->invV_vec  = &par->invV[i_cell][i_states];
        sol->C_ast_vec = &par->C_ast[i_cell][i_states];     
 		
        // e. G2EGM_1D
        G2EGM_1D::choice_and_value(sol);

    } } }

}

/////////////////////////
// 2. choice and value //
/////////////////////////

void choice_and_value(sol_struct *sol){

	int level = 3;

		logs::solve(level,"G2EGM_1D::choice_and_value\n");

	///////////////
	// 1. unpack //
	///////////////

	// a. input
	auto par = sol->par;
	auto t   = sol->t;

	auto A_vec    = sol->A_vec;
	auto invW_vec = sol->invW_vec;
	auto Q_vec    = sol->Q_vec;

	auto M_common = sol->M_common;

	// b. working memory
	auto C_vec = sol->C_vec;
	auto M_vec = sol->M_vec;
	
	// c. output
	auto V_vec     = sol->V_vec;
	auto invV_vec  = sol->invV_vec;	
	auto C_ast_vec = sol->C_ast_vec;


	///////////////////
	// 2. initialize //
	///////////////////

	for(int i_M = 0; i_M < par->NM; i_M++){
		C_ast_vec[i_M] = 0.0;
		if(par->epstein_zin == 0){
			V_vec[i_M]    = -HUGE_VAL;	
    		invV_vec[i_M] = 0.0;			
		} else {
			V_vec[i_M]    = 0.0;	
    		invV_vec[i_M] = 0.0;
		}
	}		
	

	//////////////
	// 3. EGM ////
	//////////////

		logs::solve(level+1,"EGM\n");	

	double M_min = HUGE_VAL;
	double M_max = -HUGE_VAL;
	for(int i_A = 0; i_A < par->NA; i_A++){

		// a. invert Euler-equation and endogenous M
		C_vec[i_A] = utilities::inv_marg_u(Q_vec[i_A],par);			
			
		// b. endogenous M
		M_vec[i_A] = A_vec[i_A] + C_vec[i_A];

		// c. update min/max M
		M_min = MIN(M_vec[i_A],M_min);
		M_max = MAX(M_vec[i_A],M_max);

	}


	///////////////////
	// 4. constraint //
	///////////////////

		logs::solve(level+1,"constraint\n");

	// a. next-lowest nodes
	int i_M_con = 1;
	double W = utilities::inv_trans(invW_vec[0],par);
    if(par->epstein_zin == 1){
        W = pow(W,1.0-par->sigma);
    }   	
    
    // b. all nodes on constraint
	while(i_M_con < par->NM && M_common[i_M_con] < M_min){

		// i. make the borrowing constraint bind
		C_ast_vec[i_M_con] = M_common[i_M_con] + sol->P*par->omega[t];

		// ii. value-of-choice 
		V_vec[i_M_con] = utilities::u(C_ast_vec[i_M_con],par) + par->beta*W;
	    if(par->epstein_zin == 1){
	        V_vec[i_M_con] = pow(V_vec[i_M_con],1.0/(1.0-par->sigma));
	    } 

		i_M_con++;

	}


	/////////////////////////
	// 5. upper envelope ////
	/////////////////////////
	
		logs::solve(level+1,"upper envelope\n");

	for(int i_A = 0; i_A < par->NA-1; i_A++){

        // a. M inteval and C slope
		double M_low  = M_vec[i_A];
		double M_high = M_vec[i_A+1];

		double C_low  = C_vec[i_A];
		double C_high = C_vec[i_A+1];

		double C_slope = (C_high-C_low)/(M_high-M_low);

        // a. A inteval and invW slope
		double invW_low  = invW_vec[i_A];
		double invW_high = invW_vec[i_A+1];

		double A_low  = A_vec[i_A];
		double A_high = A_vec[i_A+1];

		if(A_low > A_high){continue;}
		
		double invW_slope = (invW_high-invW_low)/(A_high-A_low);

        // c. loop through common grid
		for(int i_M = i_M_con; i_M < par->NM; i_M++){

			// i. current m
			double M = M_common[i_M];

			// ii. interpolate?
			bool interp = (M >= M_low) && (M <= M_high);			
			bool extrap_above = i_A == par->NA-1 && M > M_max;

			// iii. interpolation
			if(interp || extrap_above){

				// o. implied guess
				double C_guess = C_low + C_slope * (M - M_low);
				double A_guess = M - C_guess;

				// oo. implied post-decision value function
				double invW = invW_low + invW_slope * (A_guess - A_low);;				
				double W = utilities::inv_trans(invW,par);
			    if(par->epstein_zin == 1){
			        W = pow(W,1.0-par->sigma);
			    }   	

				// ooo. value-of-choice
				double V_guess = utilities::u(C_guess,par) + par->beta * W;
			    if(par->epstein_zin == 1){
			        V_guess = pow(V_guess,1.0/(1.0-par->sigma));
			    } 

				// oooo. update
				if(V_guess > V_vec[i_M]){
					V_vec[i_M]     = V_guess;
					C_ast_vec[i_M] = C_guess;
				}

			}

		} // common 

	} // endogenous


    //////////////////////////////////////
    // 6. force positive C slope at end //
    //////////////////////////////////////

	if(C_ast_vec[par->NM-1] <= C_ast_vec[par->NM-2]){ // fall at end of line
		
		// a. latest index with positive C slope
		int i_M_up = par->NM-2;
		while(i_M_up >= 0 && C_ast_vec[i_M_up-1] >= C_ast_vec[i_M_up]){i_M_up--;}

		// b. C slope
		double C_slope = (C_ast_vec[i_M_up] - C_ast_vec[i_M_up-1]) / (M_common[i_M_up] - M_common[i_M_up-1]); 

		// c. apply C slope at end
		for(int i_M = i_M_up+1; i_M < par->NM; i_M++){
			C_ast_vec[i_M] = C_ast_vec[i_M_up] + C_slope * (M_common[i_M] - M_common[i_M_up]);
		}

	}


    ////////////////////
    // 7. transform V //
    ////////////////////

    for(int i_M = 1; i_M < par->NM; i_M++){
        invV_vec[i_M] = utilities::trans(V_vec[i_M],par);
    }

}


} // namespace