
    // forward declarations
    namespace compute_WQ {
        void calc(sol_struct *sol);
    }

namespace model_vfi_pure {

////////////////////////
// 1. value-of-choice //
////////////////////////

double endofperiod_A(sol_struct *sol){

    auto par = sol->par;
    auto t   = sol->t;

    if(sol->z == 0){
        return sol->M-sol->C;
    } else {
        return sol->X-sol->C-sol->B-par->omega[t]*sol->P;
    }

}

double value_of_choice(sol_struct *sol){

    int level = 4;

        logs::solve(level,"value_of_choice\n");

    auto par = sol->par;
    
    // a. consumption utility
    double u = utilities::u(sol->C,par);

        logs::solve(level,"u = %g\n",u);

    // b. continuation value
    sol->A = endofperiod_A(sol);
    compute_WQ::calc(sol);
    double W = sol->W;
    if(par->epstein_zin == 1){
        W = pow(W,1.0-par->sigma);
    }

        logs::solve(level,"A = %g, W = %g\n",sol->A,W);

    // c. total
    double val = u + par->beta*W;
    if(par->epstein_zin == 1){
        val = pow(val,1.0/(1.0-par->sigma));
    } 

        logs::solve(level,"val = %g\n", val);

    return val;

}


//////////////////////////////
// 2. inequality constraint //
//////////////////////////////

double ineqcon(unsigned n, const double *choices, double *grad, void *sol_in)
{

    int level = 4;

        logs::solve(level,"ineqcon\n");

    // a. unpack
    auto sol = (sol_struct *) sol_in;
    auto par = sol->par;
    auto t   = sol->t;
    auto z   = sol->z;

    // b. choices
    sol->C = choices[0];
    if(z == 0){
        sol->B = sol->N;
    } else {
        sol->B = choices[1];
    }

    // c. gradient
    if (grad) {
        grad[0] = 1.0;
        if(z == 1){
            grad[1] = 1.0;
        }
    }

    // d. distance to constraint
    double Acon = -par->omega[t]*sol->P;
    double A    = endofperiod_A(sol);

    return Acon-A; // positive -> violated

}


////////////////
// 3. objfunc //
////////////////

double objfunc(unsigned n, const double *choices, double *grad, void *sol_in)
{

    int level = 3;

        logs::solve(level,"objfunc\n");

	// a. unpack
    auto sol = (sol_struct *) sol_in;
    auto par = sol->par;
    auto z   = sol->z;

    // b. choices
    sol->C = choices[0];
    if(z == 0){
        sol->B = sol->N;
    } else {
        sol->B = choices[1];
    }

    // c. value-of-choice        
    double obj = -value_of_choice(sol);

        logs::solve(level," C = %g, obj = %g\n",sol->C,-obj);

    // d. gradient
    if(grad){

        // C

            // i. increase C
        	sol->C = choices[0]+par->eps;
            if(z == 0){
                sol->B = sol->N;
            } else {
                sol->B = choices[1];
            }

            // ii. new value-of-choice
        	double forward = -value_of_choice(sol);

            // iii. gradient
        	grad[0] = (forward - obj)/par->eps;

        // B
        if(z == 1){

            // i. increase C
            sol->C = choices[0];
            sol->B = choices[1]+par->eps;

            // ii. new value-of-choice
            double forward = -value_of_choice(sol);

            // iii. gradient
            grad[1] = (forward - obj)/par->eps;

        }

    }

    // e. output
    return obj;

}


} // namespace