//////////////////
// 1. variables //
//////////////////

typedef struct
{
 
    // a. par
    par_struct *par;

    // b. interpolants
    linear_interp::set_struct **next_period_interp;
    linear_interp::set_struct **next_period_interp_invV;
    linear_interp::set_struct **next_period_interp_C;          
    linear_interp::set_struct **post_decision_interp;
    linear_interp::set_struct **profilled_interp;
    linear_interp::set_struct **profilled_C_interp;

    // c. info
    int t, d, z, d_plus, z_plus, i_grid;
    double M, P, N, A, B, X, C, W;
    double P_plus, Y_plus, N_plus, M_plus, X_plus;
    double M_profile, P_profile, N_profile;
    double *psi, *xi, *weights;
    double choice_ast_lag, inh_A;    
    
    // d. optimization
    nlopt_opt* opt; 
    double **lower_bounds, **upper_bounds, **choices;
    int vfi_A0;

    // e. working memory
        
        // next-period interpolation
        double *yi_P, *yi, *M_plus_vec, *X_plus_vec; 
        
        // post-decision functions
        int i_cell, i_shock;
        double *A_vec, *invV_plus_vec, *C_plus_vec;
        double *W_shocks, *Q_shocks;
        double *W_vec, *invW_vec, *Q_vec;

        // EGM
        double *M_vec, *C_vec;
        double *V_vec, *invV_vec, *C_ast_vec;
        double *M_common;

} sol_struct;

// forward declarations
namespace model {
    double objfunc(unsigned n, const double *choices, double *grad, void *sol_in);
}
namespace model_vfi_pure {
    double objfunc(unsigned n, const double *choices, double *grad, void *sol_in);
    double ineqcon(unsigned n, const double *choices, double *grad, void *sol_in);    
}

namespace sol {

//////////////
// 2. setup //
//////////////

void setup(par_struct *par, sol_struct *sol){

    // a. par
    sol->par = par;

     // b. interpolants
    sol->next_period_interp         = new linear_interp::set_struct*[par->Nd*par->Nz*par->NP];
    sol->next_period_interp_invV    = new linear_interp::set_struct*[par->Nd*par->Nz*par->NP];  
    sol->next_period_interp_C       = new linear_interp::set_struct*[par->Nd*par->Nz*par->NP];  
    sol->post_decision_interp       = new linear_interp::set_struct*[par->Nd*par->Nz*par->NP];
    sol->profilled_interp           = new linear_interp::set_struct*[par->Nd*par->Nz*par->NP];
    sol->profilled_C_interp         = new linear_interp::set_struct*[par->Nd*par->Nz*par->NP]; 

    // d. optimization
    sol->opt = new nlopt_opt[par->Nz];
    
        sol->choices      = new double*[par->Nz];
        sol->lower_bounds = new double*[par->Nz];
        sol->upper_bounds = new double*[par->Nz];

    for(int z = 0; z < par->Nz; z++){ // z = 0 -> no choices except consumption

        // i. create
        int Nchoices = 1;
        if(par->vfi_pure == 1){
            if(z == 0){
                Nchoices = 1;
            } else {
                Nchoices = 2;
            }
        }
        sol->opt[z] = nlopt_create(NLOPT_LD_MMA, Nchoices);

            sol->choices[z]      = new double[Nchoices];
            sol->lower_bounds[z] = new double[Nchoices];
            sol->upper_bounds[z] = new double[Nchoices];

        // ii. settings
        if(par->vfi_pure == 0){
            nlopt_set_min_objective(sol->opt[z], model::objfunc, sol);
        } else {
            if(z == 1){
                nlopt_add_inequality_constraint(sol->opt[z], model_vfi_pure::ineqcon, sol, 1e-8);            
            }            
            nlopt_set_min_objective(sol->opt[z], model_vfi_pure::objfunc, sol);            

        }
        nlopt_set_xtol_rel(sol->opt[z], 1e-6);
        nlopt_set_ftol_rel(sol->opt[z], 1e-6);
        nlopt_set_maxeval(sol->opt[z], 500);

    }

    // e. working memory

        // next-period interpolation
        sol->X_plus_vec = new double[par->NA];
        sol->M_plus_vec = new double[par->NA];

        sol->yi_P = new double[2*par->Nz*par->NA*3]; 
        sol->yi = new double[par->Nz*par->NA*3]; 
        
        // compute_WQ::do
        sol->invV_plus_vec = new double[par->Nz*par->NA];
        sol->C_plus_vec    = new double[par->Nz*par->NA];

        sol->W_shocks = new double[par->NA];
        sol->Q_shocks = new double[par->NA];
                 
}


////////////////
// 3. destroy //
////////////////

void destroy(sol_struct *sol){

    auto par = sol->par;

    // b. interpolants
    delete[] sol->next_period_interp;
    delete[] sol->next_period_interp_invV;
    delete[] sol->next_period_interp_C;        
    delete[] sol->post_decision_interp;   
    delete[] sol->profilled_interp; 
    delete[] sol->profilled_C_interp; 

    // d. optimizations
    for(int z = 0; z < par->Nz; z++){
        nlopt_destroy(sol->opt[z]);
        delete[] sol->choices[z];
        delete[] sol->lower_bounds[z];
        delete[] sol->upper_bounds[z];
    }
    delete[] sol->opt;
    delete[] sol->choices;
    delete[] sol->lower_bounds;
    delete[] sol->upper_bounds;

    // e. working memory

        // next-period interpolation
        delete[] sol->X_plus_vec;
        delete[] sol->M_plus_vec;
        delete[] sol->yi_P;
        delete[] sol->yi;

        // compute_WQ::do
        delete[] sol->invV_plus_vec;
        delete[] sol->C_plus_vec;

        delete[] sol->W_shocks;
        delete[] sol->Q_shocks;
        
}

} // namespace