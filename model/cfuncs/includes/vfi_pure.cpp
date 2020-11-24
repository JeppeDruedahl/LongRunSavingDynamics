namespace vfi_pure {

    // forward declations
    void solve_given_states(sol_struct *sol);

//////////////
// 1. solve //
//////////////

void solve(int t, par_struct *par, sol_struct *sol){
	
    int level = 2;

        logs::solve(level,"vfi_pure::solve\n");

    // a. non-adjuster
    #pragma omp for collapse(4)
    for(int d   = 0; d   < par->Nd; d++){
    for(int i_P = 0; i_P < par->NP; i_P++){
    for(int i_N = 0; i_N < par->NN; i_N++){
    for(int i_M = 0; i_M < par->NM; i_M++){

        int z = 0;

          logs::solve(level+1,"(z,d,i_P,i_N,i_M) = (%1d, %1d, %3d, %3d, %3d)\n",z,d,i_P,i_N,i_M);

        // a. pack
        sol->t = t;
        sol->d = d;
        sol->z = z; 

        sol->P = par->grid_P[index::d2(t,i_P,par->NP)];                
        sol->N = par->grid_N[index::d3(t,i_P,i_N,par->NP,par->NN)];
        sol->M = par->grid_M[index::d3(t,i_P,i_M,par->NP,par->NM)];

        sol->i_cell = index::d3(t,d,z,par->Nd,par->Nz);
        sol->i_grid = index::d3(i_P,i_N,i_M,par->NN,par->NM);

		// b. solve
		vfi_pure::solve_given_states(sol);

    } } } }

    // b. adjusters
    if(par->Nz == 2){

        #pragma omp for collapse(3)
        for(int d   = 0; d   < par->Nd; d++){
        for(int i_P = 0; i_P < par->NP; i_P++){
        for(int i_X = 0; i_X < par->NX; i_X++){

            int z = 1;

              logs::solve(level+1,"(z,d,i_P,i_X) = (%1d, %1d, %3d, %3d)\n",z,d,i_P,i_X);

            // a. pack
            sol->t = t;
            sol->d = d;
            sol->z = z; 

            sol->P = par->grid_P[index::d2(t,i_P,par->NP)];                
            sol->X = par->grid_X[index::d3(t,i_P,i_X,par->NP,par->NX)];

            sol->i_cell = index::d3(t,d,z,par->Nd,par->Nz);
            sol->i_grid = index::d2(i_P,i_X,par->NX);

            // b. solve
            vfi_pure::solve_given_states(sol);

        } } }
        
    }

}


///////////////////////////
// 2. solve given states //
///////////////////////////

void solve_given_states(sol_struct *sol){

    int level = 3;

        logs::solve(level,"vfi_pure::solve_given_states\n");

    /////////////////////
    // 1. preparations //
    /////////////////////

	auto par = sol->par;

    auto t = sol->t;
    auto z = sol->z;
    auto M = sol->M;
    auto N = sol->N;    
    auto X = sol->X;
    auto P = sol->P;    

    auto i_cell = sol->i_cell;
    auto i_grid = sol->i_grid;  

    auto opt          = sol->opt[z]; 
    auto choices      = sol->choices[z];
    auto lower_bounds = sol->lower_bounds[z]; 
    auto upper_bounds = sol->upper_bounds[z];

    //////////////
    // 2. solve //
    //////////////

    // a. pooled wealth
    if((z == 0 && M <= -par->omega[t]*P) || (z == 1 && X <= 0)){
        par->invV[i_cell][i_grid]  = 0.0;
        par->V[i_cell][i_grid]     = -HUGE_VAL;
        par->B_ast[i_cell][i_grid] = 0.0;
        par->C_ast[i_cell][i_grid] = 0.0;
        return;
    }

    // b. bounds
    if(z == 0) {
        lower_bounds[0] = 0.0;
        upper_bounds[0] = M+par->omega[t]*P;
    } else {
        lower_bounds[0] = 0.0;
        upper_bounds[0] = X;
        lower_bounds[1] = 0.0;
        upper_bounds[1] = X;                
    }

    nlopt_set_lower_bounds(opt, lower_bounds);
    nlopt_set_upper_bounds(opt, upper_bounds);

    // c. solve with multistart
	double V_ast = -HUGE_VAL;
    double C_ast = NAN;
    double B_ast = NAN;

        // multistart setings
        int Nk_C = 3;
        int Nk_B = 1;
        if(z == 1){
            Nk_C = 2;
            Nk_B = 3;    
        }

        if(par->no_multistart == 1){
            Nk_C = 1;
            Nk_B = 1;
        }

        if(z == 0){
            logs::solve(level+1,"M = %g, N = %g, P = %g\n",M,N,P);
        } else {
            logs::solve(level+1,"X = %g, P = %g\n",X,P);
        }


    for(int k_C = 0; k_C < Nk_C; k_C++){
    for(int k_B = 0; k_B < Nk_B; k_B++){ 

        // i. initial guess
        if(par->no_multistart == 1){

            if(z == 0){
                choices[0] = 0.95*(M+par->omega[t]*P);;
            } else {
                choices[0] = 0.95*X;
                choices[1] = 0.01*X;
            }

        } else {

            double phi_C = (double)(k_C+1)/(double)Nk_C;
            if(z == 0){
                choices[0] = phi_C*(M+par->omega[t]*P);
            } else {
                double phi_B = (double)k_B/(double)(Nk_B-1);         
                choices[0] = ((phi_C*phi_C)/(phi_C + phi_B))*X;
                choices[1] = ((phi_B*phi_B)/(phi_C + phi_B))*X;  
            }

        }

        // ii. optimizer
    	double minf;

            if(z == 0){
                logs::solve(level+1,"(%d):  C = %g\n",k_C,choices[0]);
            } else {
                logs::solve(level+1,"(%d,%d): C = %g, B = %g\n",k_C,k_B,choices[0],choices[1]);
            }

    	auto nlopt_result = nlopt_optimize(opt, choices, &minf);

            if(z == 0){
                logs::solve(level+1," C = %g, V = %g\n",choices[0],-minf);
            } else {
                logs::solve(level+1," C = %g, B = %g, V = %g\n",choices[0],choices[1],-minf);
            }

            assert(0,nlopt_result >= 0,"optimizer failed! %d", nlopt_result);

        // iii. check for max
    	double V_now = -minf;
    	if(V_now > V_ast){
    		V_ast = V_now;
            C_ast = choices[0];
            if(z == 1){
                B_ast = choices[1];            
            } else{
                B_ast = sol->N;
            }
    	}

    } }

    // d. save value and choice   
    par->V[i_cell][i_grid]    = V_ast;
    par->invV[i_cell][i_grid] = utilities::trans(par->V[i_cell][i_grid],par);

    par->C_ast[i_cell][i_grid] = C_ast;
    par->B_ast[i_cell][i_grid] = B_ast;

}

} // namespace