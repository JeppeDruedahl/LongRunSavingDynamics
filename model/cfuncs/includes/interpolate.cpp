namespace interpolate {

//////////////
// 1. setup //
//////////////

void setup_next_period(int t, par_struct *par, sol_struct *sol, int NA){

	int level = 3;

	////////////////////
	// 1. next-period //
	////////////////////

	for(int d   = 0; d   < par->Nd; d++){ 	// discrete states
	for(int z   = 0; z   < par->Nz; z++){ 	// discrete choices
	for(int i_P = 0; i_P < par->NP; i_P++){ // continuous state (determine M,N,X grids)

		// a. grids
			
			// i. dimensions
			int dimx;
			int Nx[2];
			double *x[2];

			if(z == 1){
				
				dimx = 1; 
				Nx[0] = par->NX;	
				x[0] = &par->grid_X[index::d3(t+1,i_P,0,par->NP,par->NX)];
			
			} else if(par->NN == 1){
			
				dimx = 1; 
				Nx[0] = par->NM;	
				x[0] = &par->grid_M[index::d3(t+1,i_P,0,par->NP,par->NM)];

			} else {

				dimx = 2;
				Nx[0] = par->NM;
				Nx[1] = par->NN;	
				x[0] = &par->grid_M[index::d3(t+1,i_P,0,par->NP,par->NM)];
				x[1] = &par->grid_N[index::d3(t+1,i_P,0,par->NP,par->NN)];

			}

		// b. values
		int dimy = 2;
		if(par->vfi_pure == 1){ // no need to interpolate consumption
			dimy = 1;
		}
		int i_cell = index::d3(t+1,d,z,par->Nd,par->Nz);
		int i_grid;
		if(z == 0){
			i_grid = index::d3(i_P,0,0,par->NN,par->NM);
		} else {
			i_grid = index::d2(i_P,0,par->NX);				
		}

		int i_interp = index::d3(d, z,i_P,par->Nz,par->NP);
		if(par->no_interp_vec == 1 && par->vfi_pure == 0){
			
			double *y[1];

			y[0] = &par->invV[i_cell][i_grid]; 
			sol->next_period_interp_invV[i_interp] = new linear_interp::set_struct;
			auto interp_invV = sol->next_period_interp_invV[i_interp];
            linear_interp::create(interp_invV, dimx, 1, Nx, x, y, 1);

			y[0] = &par->C_ast[i_cell][i_grid]; 
			sol->next_period_interp_C[i_interp] = new linear_interp::set_struct;
			auto interp_C = sol->next_period_interp_C[i_interp];
            linear_interp::create(interp_C, dimx, 1, Nx, x, y, 1);

		} else {
					
			double *y[2];
			y[0] = &par->invV[i_cell][i_grid];
			y[1] = &par->C_ast[i_cell][i_grid];

			// c. create
			sol->next_period_interp[i_interp] = new linear_interp::set_struct;
			auto interp = sol->next_period_interp[i_interp];
            linear_interp::create(interp, dimx, dimy, Nx, x, y, NA);

		}

        	logs::solve(level+1,"next-period (d=%d,z=%d,i_P=%d)\n",d,z,i_P);

	} } }

}

void setup_post_decison(int t, par_struct *par, sol_struct *sol){

	int level = 3;

	//////////////////////
	// 2. post-decision //
	//////////////////////

    int z = 0; // z is not a post-decision state

	for(int d   = 0; d   < par->Nd; d++){	 // discrete states
	for(int i_P = 0; i_P < par->NP; i_P++){	 // continuous state (determine M,N,X grids)

		// a. grids

			// i. dimensions
			int dimx;
			int Nx[2];
			double *x[2];
			if(par->NN == 1){
			
				dimx = 1;
				Nx[0] = par->NA;					
				x[0] = &par->grid_A[index::d3(t,i_P,0,par->NP,par->NA)];

			} else {

				dimx = 2;
				Nx[0] = par->NA;
				Nx[1] = par->NN;					
				x[0] = &par->grid_A[index::d3(t,i_P,0,par->NP,par->NA)];
				x[1] = &par->grid_N[index::d3(t,i_P,0,par->NP,par->NN)];

			}

		// b. values
		int i_cell = index::d3(0,d,z,par->Nd,1);
		int i_grid = index::d3(i_P,0,0,par->NN,par->NA);
		int dimy = 1;
		double *y[1];
		y[0] = &par->invW[i_cell][i_grid];
			
		// c. create
		int i_interp = index::d3(d,z,i_P,par->Nz,par->NP);
		sol->post_decision_interp[i_interp] = new linear_interp::set_struct; 
		auto interp = sol->post_decision_interp[i_interp];

        linear_interp::create(interp, dimx, dimy, Nx, x, y, 1);           

			logs::solve(level+1,"post-decision (d=%d,i_P=%d)\n",d,i_P);

	} }

}

void setup_profiled(int t, par_struct *par, sol_struct *sol){

	int level = 3;

	/////////////////
	// 3. profiled //
	/////////////////

	for(int d = 0; d < par->Nd; d++){		// discrete state
	for(int i_P = 0; i_P < par->NP; i_P++){ // continuous state (determine M,N,X grids)

		int z = 0;

		// a. grids
	
			// i. dimensions
			int dimx;
			int Nx[2];
			double *x[2];

			if(par->NN == 1){
			
				dimx = 1; 
				Nx[0] = par->NM;					
				x[0]  = &par->grid_M[index::d3(t,i_P,0,par->NP,par->NM)];

			} else {

				dimx = 2;
				Nx[0] = par->NM;
				Nx[1] = par->NN;					
				x[0] = &par->grid_M[index::d3(t,i_P,0,par->NP,par->NM)];
				x[1] = &par->grid_N[index::d3(t,i_P,0,par->NP,par->NN)];

			}

		// b. values
		int i_cell = index::d3(t,d,z,par->Nd,par->Nz);
		int i_grid = index::d3(i_P,0,0,par->NN,par->NM);

		int dimy = 1;
		double *y[1];
		y[0] = &par->invV[i_cell][i_grid];

		// c. create
		int i_interp = index::d3(d,z,i_P,par->Nz,par->NP);
		sol->profilled_interp[i_interp] = new linear_interp::set_struct;
		auto interp = sol->profilled_interp[i_interp];

        linear_interp::create(interp, dimx, dimy, Nx, x, y, 1);

        // d. create for C
		sol->profilled_C_interp[i_interp] = new linear_interp::set_struct;
		auto interp_C = sol->profilled_C_interp[i_interp]; 

		y[0] = &par->C_ast[i_cell][i_grid];
        linear_interp::create(interp_C, dimx, dimy, Nx, x, y, 1);

			logs::solve(level+1,"profilled (d=%d,i_P=%d)\n",d,i_P);

	} }
 
}

void setup(int t, par_struct *par, sol_struct *sol){

	int level = 2;

		logs::solve(level,"interpolate::setup\n");

	setup_next_period(t,par,sol,par->NA);
	setup_post_decison(t,par,sol);
	setup_profiled(t,par,sol);

}

void destroy_next_period(int t, par_struct *par, sol_struct *sol){

	////////////////////
	// 1. next-period //
	////////////////////

	for(int d   = 0; d   < par->Nd; d++){ 	// discrete states
	for(int z   = 0; z   < par->Nz; z++){ 	// discrete choices		                    
	for(int i_P = 0; i_P < par->NP; i_P++){ // continuous state (determine M,N,X grids)

		int i_interp = index::d3(d,z,i_P,par->Nz,par->NP);
		if(par->no_interp_vec == 1 && par->vfi_pure == 0){
			linear_interp::destroy(sol->next_period_interp_invV[i_interp]);
			linear_interp::destroy(sol->next_period_interp_C[i_interp]);
		} else {
			linear_interp::destroy(sol->next_period_interp[i_interp]);
		}

	} } }

}

void destroy_post_decision(int t, par_struct *par, sol_struct *sol){

	//////////////////////
	// 2. post-decision //
	//////////////////////

	for(int d   = 0; d   < par->Nd; d++){	// discrete states
	for(int i_P = 0; i_P < par->NP; i_P++){	// continuous state (determine M,N,X grids)		

		int i_interp = index::d3(d,0,i_P,par->Nz,par->NP);	
		linear_interp::destroy(sol->post_decision_interp[i_interp]);
	}
	}

}

void destroy_profiled(int t, par_struct *par, sol_struct *sol){

	/////////////////
	// 3. profiled //
	/////////////////

	for(int d = 0; d < par->Nd; d++){		// discre states
	for(int i_P = 0; i_P < par->NP; i_P++){	// continuous state (determine M,N,X grids)

		int i_interp = index::d3(d,0,i_P,par->Nz,par->NP);
		linear_interp::destroy(sol->profilled_interp[i_interp]);
		linear_interp::destroy(sol->profilled_C_interp[i_interp]);

	} }
 
}

void destroy(int t, par_struct *par, sol_struct *sol){

	int level = 2;

		logs::solve(level,"interpolate::destroy\n");

	destroy_next_period(t,par,sol);
	destroy_post_decision(t,par,sol);
	destroy_profiled(t,par,sol);
			
}


/////////////
// 2. call //
/////////////

void next_period(sol_struct *sol, double *invV_plus_in, double *C_plus_in){

	// a. unpack
	auto par    = sol->par;
	auto t      = sol->t;
	auto d_plus = sol->d_plus;
	auto P_plus = sol->P_plus;

	auto P_plus_vec = &par->grid_P[index::d2(t+1,0,par->NP)];

	// b. find left index and reldiff for P
	int i_P_left = linear_interp::binary_search(0, par->NP, P_plus_vec, P_plus);
	double P_reldiff = (P_plus - P_plus_vec[i_P_left])/(P_plus_vec[i_P_left+1]-P_plus_vec[i_P_left]);

    // c. interpolate
    for(int z_plus = 0; z_plus < par->Nz; z_plus++){

		for(int add = 0; add < 2; add++){

			// i. P index
			int i_P = i_P_left + add;

	    	// ii. interpolants
			int i_interp = index::d3(d_plus,z_plus,i_P,par->Nz,par->NP);	    	
	    	auto interp = sol->next_period_interp[i_interp];

	    		auto interp_invV = sol->next_period_interp_invV[i_interp];
	    		auto interp_C    = sol->next_period_interp_C[i_interp];

	    	// iii. xi and xi_vec
			if(par->no_interp_vec == 1){
	    	for(int i_A = 0; i_A < par->NA; i_A++){

				double xi[2];
		    	if(z_plus == 1){

		    		xi[0]  = sol->X_plus_vec[i_A];

		    	} else {

		    		xi[0] = sol->M_plus_vec[i_A];
		    		xi[1] = sol->N_plus;

		    	}

		    	// iv. evaluate   
		    	auto yi = &sol->yi_P[add*par->NA*2 + 0*par->NA + i_A];
		    	linear_interp::evaluate(interp_invV,yi,xi,nullptr);
				yi = &sol->yi_P[add*par->NA*2 + 1*par->NA + i_A];
		    	linear_interp::evaluate(interp_C,yi,xi,nullptr);

		    }
	    	} else {
	
				double xi[2];
				double *xi_vec;
		    	if(z_plus == 1){

		    		xi[0]  = sol->X_plus_vec[0];
		    		xi_vec = sol->X_plus_vec;

		    	} else {

		    		xi[0] = sol->M_plus_vec[0];
		    		xi[1] = sol->N_plus;

		    		xi_vec = sol->M_plus_vec;

		    	}

		    	// iv. evaluate   
		    	auto yi = &sol->yi_P[add*par->NA*2];
		    	linear_interp::evaluate(interp,yi,xi,xi_vec);

	    	}

	    }

	    // v. interpolate over P
    	int index_left = 0*par->NA*2;
    	int index_right = 1*par->NA*2;	    	 
    	int index = z_plus*par->NA*2;
	    for(int i = 0; i < 2*par->NA; i++){
	    	double left = sol->yi_P[index_left+i];
	    	double right = sol->yi_P[index_right+i];	    
	    	sol->yi[index + i] = left + P_reldiff*(right-left);
    	}

    }

    // c. move in to output vectors
	for(int z_plus = 0; z_plus < par->Nz; z_plus++){
    for(int i_A = 0; i_A < par->NA; i_A++){
				
		auto invV_plus = &invV_plus_in[z_plus*par->NA];
		auto C_plus = &C_plus_in[z_plus*par->NA];

		invV_plus[i_A] = sol->yi[z_plus*par->NA*2 + i_A];
		C_plus[i_A]    = sol->yi[z_plus*par->NA*2 + par->NA + i_A];

    } }

}

void next_period_vfi_pure(sol_struct *sol, double *invV_plus){

	// a. unpack
	auto par    = sol->par;
	auto t      = sol->t;
	auto d_plus = sol->d_plus;
	auto P_plus = sol->P_plus;

	auto P_plus_vec = &par->grid_P[index::d2(t+1,0,par->NP)];

	// b. find left index and reldiff for P
	int i_P_left = linear_interp::binary_search(0, par->NP, P_plus_vec, P_plus);
	double P_reldiff = (P_plus - P_plus_vec[i_P_left])/(P_plus_vec[i_P_left+1]-P_plus_vec[i_P_left]);

    // c. interpolate
    for(int z_plus = 0; z_plus < par->Nz; z_plus++){

		for(int add = 0; add < 2; add++){

			// i. P index
			int i_P = i_P_left + add;

	    	// ii. interpolants
			int i_interp = index::d3(d_plus,z_plus,i_P,par->Nz,par->NP);	    	
	    	auto interp  = sol->next_period_interp[i_interp];

	    	// iii. xi and xi_vec
			double xi[2];
	    	if(z_plus == 1){

	    		xi[0]  = sol->X_plus;

	    	} else {

	    		xi[0] = sol->M_plus;
	    		xi[1] = sol->N_plus;

	    	}

	    	// iv. evaluate   
	    	linear_interp::evaluate(interp,&sol->yi_P[add],xi,nullptr);

	    }

	    // v. interpolate over P
	    double left = sol->yi_P[0];
		double right = sol->yi_P[1];	    
		invV_plus[z_plus] = left + P_reldiff*(right-left);

    }

}

void post_decision(sol_struct *sol, double *invW){

	auto par   = sol->par;
	auto t     = sol->t;
	auto d     = sol->d;
	auto P     = sol->P;	

	auto P_vec = &par->grid_P[index::d2(t,0,par->NP)];

	// a. find left index and reldiff for P
	int i_P_left = linear_interp::binary_search(0, par->NP, P_vec, P);
	double P_reldiff = (P - P_vec[i_P_left])/(P_vec[i_P_left+1]-P_vec[i_P_left]);

    // b. xi
   	double xi[2];
   	xi[0] = sol->A;
	if(par->NN > 1){
   		xi[1] = sol->B;
   	}

   	// c. evaluate	
	double invW_P[2];
	for(int add = 0; add < 2; add++){

		// i. P index
		int i_P = i_P_left + add;
 		
 		// ii. interpolate
		int i_interp = index::d3(d,0,i_P,par->Nz,par->NP);	    	 		
   		auto interp = sol->post_decision_interp[i_interp];
		linear_interp::evaluate(interp,&invW_P[add],xi,nullptr);	

   }
   // d. interpolate over P
   *invW = invW_P[0] + P_reldiff * (invW_P[1]-invW_P[0]);
   *invW = MAX(*invW, 0);

}

void profilled(sol_struct *sol, double *invV){

	auto par   = sol->par;
	auto t     = sol->t;
	auto d     = sol->d;
	auto P     = sol->P;

	auto P_vec = &par->grid_P[index::d2(t,0,par->NP)];

	// a. find left index and reldiff for P
	int i_P_left = linear_interp::binary_search(0, par->NP, P_vec, P);
	double P_reldiff = (P - P_vec[i_P_left])/(P_vec[i_P_left+1]-P_vec[i_P_left]);

    // b. xi
   	double xi[2];
   	xi[0] = sol->M_profile;
   	if(par->NN > 1){
   		xi[1] = sol->N_profile;
   	}

   	// c. evaluate	
	double invV_P[2];
	for(int add = 0; add < 2; add++){

		// i. P index
		int i_P = i_P_left + add;
 		
 		// ii. interpolate
		int i_interp = index::d3(d,0,i_P,par->Nz,par->NP); 		
   		auto interp = sol->profilled_interp[i_interp];
   		linear_interp::evaluate(interp,&invV_P[add],xi,nullptr);

   }

   // d. interpolate over P
   *invV = invV_P[0] + P_reldiff * (invV_P[1]-invV_P[0]);
   *invV = MAX(*invV, 0);

}

void profilled_C(sol_struct *sol, double *C){

	auto par   = sol->par;
	auto t     = sol->t;
	auto d     = sol->d;
	auto P     = sol->P;	
	auto P_vec = &par->grid_P[index::d2(t,0,par->NP)];

	// a. find left index and reldiff for P
	int i_P_left = linear_interp::binary_search(0, par->NP, P_vec, P);
	double P_reldiff = (P - P_vec[i_P_left])/(P_vec[i_P_left+1]-P_vec[i_P_left]);

    // b. xi
   	double xi[2];
   	xi[0] = sol->M_profile;
   	if(par->NN > 1){
   		xi[1] = sol->N_profile;
   	}

   	// c. evaluate	
	double C_P[2];
	for(int add = 0; add < 2; add++){

		// i. P index
		int i_P = i_P_left + add;
 		
 		// ii. interpolate
		int i_interp = index::d3(d,0,i_P,par->Nz,par->NP); 		
   		auto interp = sol->profilled_C_interp[i_interp];
   		linear_interp::evaluate(interp,&C_P[add],xi,nullptr);

   }

   // d. interpolate over P
   *C = C_P[0] + P_reldiff * (C_P[1]-C_P[0]);
   *C = MAX(*C, 0);

}

} // namespace