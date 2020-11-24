void calc_last_period(par_struct *par, sol_struct *sol);

void last_period(int t, par_struct *par, sol_struct *sol)
{
    
    int level = 2;
	
	// a. z == 0
    #pragma omp for collapse(4)
    for(int d   = 0; d   < par->Nd; d++){
    for(int i_P = 0; i_P < par->NP; i_P++){
    for(int i_N = 0; i_N < par->NN; i_N++){
    for(int i_M = 0; i_M < par->NM; i_M++){

    	int z = 0;
    	
          logs::solve(level+1,"(z,d,i_P,i_N,i_M) = (%1d, %1d, %3d, %3d, %3d)\n",z,d,i_P,i_N,i_M);

        // i. pack
        sol->t = t;
        sol->d = d;
        sol->z = z; 

        sol->P = par->grid_P[index::d2(t,i_P,par->NP)];                
        double N = par->grid_N[index::d3(t,i_P,i_N,par->NP,par->NN)];
        double M = par->grid_M[index::d3(t,i_P,i_M,par->NP,par->NM)];
        sol->X = M+N;

        sol->i_cell = index::d3(t,d,z,par->Nd,par->Nz);
        sol->i_grid = index::d3(i_P,i_N,i_M,par->NN,par->NM);

		// ii. solve
		calc_last_period(par,sol);

    } } } }

    // b. z == 1
    if(par->Nz == 2){

        #pragma omp for collapse(3)
        for(int d   = 0; d   < par->Nd; d++){
        for(int i_P = 0; i_P < par->NP; i_P++){
        for(int i_X = 0; i_X < par->NX; i_X++){

        	int z = 1;

                logs::solve(level+1,"(d,z,i_P,i_X) = (%1d, %1d, %3d,%3d)\n",d,z,i_P,i_X);

        	// i. pack
            sol->t = t;
            sol->d = d;
            sol->z = z;

            sol->P = par->grid_P[index::d2(t,i_P,par->NP)];
            sol->X = par->grid_X[index::d3(t,i_P,i_X,par->NP,par->NX)];

            sol->i_cell = index::d3(t,d,z,par->Nd,par->Nz);
            sol->i_grid = index::d2(i_P,i_X,par->NX);

    		// ii. solve
    		calc_last_period(par,sol);

        } } }
    }

}

void calc_last_period(par_struct *par, sol_struct *sol)
{

    int level = 2;

    double C_TR, C_TR_plus, V_TR, V_TR_plus;

    // a. corner solution, zeta = 0
    if(par->zeta == 0){

        C_TR = sol->X;

        if(par->epstein_zin == 0){
            V_TR = pow(C_TR,1.0-par->sigma)/(1.0-par->sigma);
        } else {
            V_TR = pow((1.0-par->beta),1.0/(1.0-par->sigma))*C_TR;
        }

    // b. interior solution, zeta > 0
    } else {

        // i. consumption
        
            double nom = par->gamma_1*(sol->X+(1.0/par->Rret+par->gamma_0/par->Rret)*par->kappa*sol->P);
            double denom = pow(par->Rret,-1.0)*pow(par->Rret*par->beta*par->zeta,1.0/par->sigma) + par->gamma_1;

        C_TR = nom/denom;
        if(C_TR < sol->X){
            C_TR_plus = pow(par->beta*par->Rret*par->zeta,1.0/par->sigma)*C_TR;
        } else {
            C_TR = sol->X;
            C_TR_plus = par->gamma_1*(1.0+par->gamma_0)*par->kappa*sol->P;
        }

        // ii. value
        if(par->epstein_zin == 0){
            V_TR_plus = par->zeta*par->gamma_2*pow(C_TR_plus,1.0-par->sigma)/(1.0-par->sigma);
            V_TR = pow(C_TR,1.0-par->sigma)/(1.0-par->sigma) + par->beta*V_TR_plus;
        } else {
            V_TR_plus = pow((1.0-par->beta)*par->zeta*par->gamma_2,1.0/(1.0-par->sigma))*C_TR_plus;
            V_TR = (1.0-par->beta)*pow(C_TR,1.0-par->sigma) + par->beta*pow(V_TR_plus,1.0-par->sigma);
            V_TR = pow(V_TR,1.0/(1.0-par->sigma));
        }
        //V_TR = -10.0;

    }

        logs::solve(level+1,"X = %g, C_TR = %g, V_TR = %g\n",sol->X,C_TR,V_TR);

	// c. save
	int i_cell = sol->i_cell;
	int i_grid = sol->i_grid;
	 
    par->C_ast[i_cell][i_grid] = C_TR;
    par->B_ast[i_cell][i_grid] = 0.0;
    
    par->V[i_cell][i_grid]    = V_TR;
    par->invV[i_cell][i_grid] = utilities::trans(V_TR,par);

}