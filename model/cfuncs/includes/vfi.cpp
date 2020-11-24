namespace vfi {

    // forward declations
    void solve_given_states(sol_struct *sol);

//////////////
// 1. solve //
//////////////

void solve(int t, par_struct *par, sol_struct *sol){
	
    int level = 2;

        logs::solve(level,"vfi::solve\n");

    // parallel
    #pragma omp for collapse(3)
    for(int d   = 0; d   < par->Nd; d++){
    for(int z   = 1; z   < par->Nz; z++){   // z = 0 is solved by EGM alone 
    for(int i_P = 0; i_P < par->NP; i_P++){

        sol->choice_ast_lag = NAN;

    for(int i_X = 0; i_X < par->NX; i_X++){

            logs::solve(level+1,"(d,z,i_P,i_X) = (%1d, %1d, %3d,%3d)\n",d,z,i_P,i_X);

    	// a. pack
        sol->t = t;
        sol->d = d;
        sol->z = z;

        sol->P = par->grid_P[index::d2(t,i_P,par->NP)];
        sol->X = par->grid_X[index::d3(t,i_P,i_X,par->NP,par->NX)];
        
        sol->i_cell = index::d3(t,d,z,par->Nd,par->Nz);
        sol->i_grid = index::d2(i_P,i_X,par->NX);

		// b. solve
		vfi::solve_given_states(sol);

    } } } }

}


///////////////////////////
// 2. solve given states //
///////////////////////////

void solve_given_states(sol_struct *sol)
{

    int level = 3;

        logs::solve(level,"vfi::solve_given_states\n");

    /////////////////////
    // 1. preparations //
    /////////////////////

	auto par = sol->par;

    auto z = sol->z;

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
    if(sol->X <= 0){
        par->invV[i_cell][i_grid]   = 0.0;
        par->V[i_cell][i_grid]      = -HUGE_VAL;
        par->B_ast[i_cell][i_grid]  = 0.0;
        par->C_ast[i_cell][i_grid]  = 0.0;
        return;
    }

    // b. bounds
    lower_bounds[0] = 0.0;
    upper_bounds[0] = 1.0;

    nlopt_set_lower_bounds(opt, lower_bounds);
    nlopt_set_upper_bounds(opt, upper_bounds);

    // c. solve with multistart
	double V_ast = -HUGE_VAL;
	double choice_ast = NAN;
    for(int k = 0; k < par->Nk+2; k++){

        // i. guess
        sol->vfi_A0 = 0;
        if(k < par->Nk){
            choices[0] = (double)k / (double)par->Nk;
        } else if(k == par->Nk){
            choices[0] = 0.99;
        } else {
            if(std::isnan(sol->choice_ast_lag) == false){
                choices[0] = sol->choice_ast_lag;
            } else {
                continue;
            }
        }
         
            logs::solve(level+1,"choice = %g, X = %g, P = %g\n",choices[0],sol->X,sol->P);

        // ii. optimizer
    	double minf;
    	auto nlopt_result = nlopt_optimize(opt, choices, &minf);

            logs::solve(level+1,"choice = %g\n",choices[0]);

        logs::solve(level+1,"done");

            assert(0,nlopt_result >= 0,"optimizer failed!");

        // iii. check for max
    	double V_now = -minf;
    	if(V_now > V_ast){
    		V_ast   = V_now;
            choice_ast = choices[0];
    	}

    }

    // d. save value and choice
 
        logs::solve(level+1,"saving");
    
    par->V[i_cell][i_grid]    = V_ast;
    par->invV[i_cell][i_grid] = utilities::trans(par->V[i_cell][i_grid],par);

    sol->B                     = choice_ast*sol->X;
    par->B_ast[i_cell][i_grid] = sol->B;
    par->C_ast[i_cell][i_grid] = model::consmption_choice(sol);

    sol->choice_ast_lag = choice_ast;

}

} // namespace