
    // forward declarations
    namespace interpolate {
        void profilled(sol_struct *sol, double *invV);
        void profilled_C(sol_struct *sol, double *C);    
        void post_decision(sol_struct *sol, double *invW);
    }

namespace model {

////////////////////////
// 1. value-of-choice //
////////////////////////

double calc_value_of_choice(sol_struct *sol){

    int level = 5;

    auto par = sol->par;

    // a. consumption utility
    double u = utilities::u(sol->C,par);

        logs::solve(level,"u = %g\n",u);

    // b. continuation value
    double invW;
    interpolate::post_decision(sol, &invW);
    double W = utilities::inv_trans(invW,par);
    if(par->epstein_zin == 1){
        W = pow(W,1.0-par->sigma);
    }

        logs::solve(level,"W = %g\n",W);

    // c. total
    double val = u + par->beta*W;
    if(par->epstein_zin == 1){
        val = pow(val,1.0/(1.0-par->sigma));
    } 

        logs::solve(level,"val = %g\n", val);

    return val;

}

double value_of_choice(sol_struct *sol){

    int level = 5;

    auto par = sol->par;
    auto t   = sol->t;

    // a. profile values
    sol->M_profile = sol->X - sol->B - sol->P*par->omega[t];
    sol->P_profile = sol->P;
    sol->N_profile = sol->B;

        logs::solve(level,"M = %g, N = %g\n",sol->M_profile,sol->N_profile);

    // b. consumption choice
    if(sol->vfi_A0 == 0){
        interpolate::profilled_C(sol, &sol->C);
    } else {
        sol->C = sol->M_profile;
    }

        logs::solve(level,"C = %g\n",sol->C);

    // c. post-decision states
    sol->A = sol->M_profile - sol->C;

    // d. value of choice
    return calc_value_of_choice(sol);

}


///////////////////////////
// 2. consumption choice //
///////////////////////////
 
double consmption_choice(sol_struct *sol){
 
    auto par = sol->par;
    auto t = sol->t;
 
    // a. profile values
    sol->M_profile = sol->X - sol->B - sol->P*par->omega[t];
    sol->P_profile = sol->P;
    sol->N_profile = sol->B;
 
    // b. calculate
    double val;
    interpolate::profilled_C(sol, &val);
    return val;
 
}


////////////////
// 3. objfunc //
////////////////

double objfunc(unsigned n, const double *choices, double *grad, void *sol_in){

    int level = 5;

        logs::solve(level,"objfunc\n");

	// a. unpack
    auto sol = (sol_struct *) sol_in;
    auto par = sol->par;

    // b. choices
   	sol->B = choices[0]*sol->X;
	
	// c. value-of-choice    
    double obj = -value_of_choice(sol);

    // d. gradient
    if(grad){

        // i. increase b
    	sol->B = sol->X*(choices[0]+par->eps);

        // ii. new value-of-choice
    	double forward = -value_of_choice(sol);

        // iii. gradient
    	grad[0] = (forward - obj)/par->eps;

    }

    // e. output
    return obj;

}

} // namespace