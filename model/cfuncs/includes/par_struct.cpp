//////////////////
// 1. variables //
//////////////////

typedef struct{

    // a. discrete states and choices
    int Nd, Nz, Nk;

    // b. grids
    int NM, NP, NN, NX, NA;
    double *grid_M, *grid_P, *grid_N, *grid_X, *grid_A;

    // c. demographics
    int T, TR;

    // d. preferences
    int epstein_zin;
    double beta, sigma, rho, zeta;
    double gamma_0, gamma_1, gamma_2;

    // e. income process
    double *G, kappa;

    // f. shocks
    int Nshocks_max, *Nshocks;
    double *psi, *xi, *weights;
    
    // g. assets
    double R, Rret, Rneg, Rb, *omega, adjcost, adjcostinf;

    // h. inheritance
    double inh_A, *inh_p, *inh_fac;

    // i. simulate
    int T0sim, Tendsim, forget_inh, Nsim, save_to_disc;
    double *xi_sim, *psi_sim, *A_shock, *inh_pval, *adjcost_pval, *P_ini_dist, *N_ini_dist;
    double MPC_add;
    int Neuler, euler_t;
    double *euler_M, *euler_N, *euler_P;

    // j. technical
    int tmin, sim_IRF, vfi_pure, no_multistart, no_interp_vec;
    double eps, cutoff;

    // output
        
        // solve
        double **C_ast, **B_ast, **A_ast, **invV;
        double **V, **W, **invW, **Q;
        double **C_EGM, **M_EGM;

        // simulate
        int *d, *z;
        double *M, *P, *Y, *N, *X, *H, *C, *A, *B, *MPC, *simV;
        double *P_IRF, *A_IRF, *B_IRF;

        // euler
        double *euler_C, *euler_error;

} par_struct;


namespace par {

//////////////
// 2. setup //
//////////////

void setup(par_struct *par, mxArray *plhs[], const mxArray *prhs[], int type){

    par->eps = 1e-5;

    ///////////////
    // 1. inputs //
    ///////////////

        int in = 0; // first input

    // a. discrete states and choices
    par->Nd = (int) mxGetScalar(mxGetField(prhs[in],0,"Nd"));
    par->Nz = (int) mxGetScalar(mxGetField(prhs[in],0,"Nz"));        
    par->Nk = (int) mxGetScalar(mxGetField(prhs[in],0,"Nk")); 

    // b. grids
    par->NM = (int) mxGetScalar(mxGetField(prhs[in],0,"NM"));
    par->NP = (int) mxGetScalar(mxGetField(prhs[in],0,"NP"));  
    par->NN = (int) mxGetScalar(mxGetField(prhs[in],0,"NN"));
    par->NX = (int) mxGetScalar(mxGetField(prhs[in],0,"NX"));    
    par->NA = (int) mxGetScalar(mxGetField(prhs[in],0,"NA")); 

    par->grid_M = (double*) mxGetPr(mxGetField(prhs[in],0,"grid_M")); 
    par->grid_P = (double*) mxGetPr(mxGetField(prhs[in],0,"grid_P"));
    par->grid_N = (double*) mxGetPr(mxGetField(prhs[in],0,"grid_N"));
    par->grid_X = (double*) mxGetPr(mxGetField(prhs[in],0,"grid_X")); 
    par->grid_A = (double*) mxGetPr(mxGetField(prhs[in],0,"grid_A"));

    // c. demographics
    par->T   = (int) mxGetScalar(mxGetField(prhs[in],0,"T"));
    par->TR  = (int) mxGetScalar(mxGetField(prhs[in],0,"TR"));

    // d. preferences
    par->epstein_zin = (int) mxGetScalar(mxGetField(prhs[in],0,"epstein_zin")); 
    
    par->beta  = (double) mxGetScalar(mxGetField(prhs[in],0,"beta"));
    par->sigma = (double) mxGetScalar(mxGetField(prhs[in],0,"sigma"));
    par->rho   = (double) mxGetScalar(mxGetField(prhs[in],0,"rho"));
    par->zeta  = (double) mxGetScalar(mxGetField(prhs[in],0,"zeta"));

    // e. income process
    par->G = (double*) mxGetPr(mxGetField(prhs[in],0,"G"));
    par->kappa = (double) mxGetScalar(mxGetField(prhs[in],0,"kappa"));

    // f. shocks
    par->Nshocks_max = (int) mxGetScalar(mxGetField(prhs[in],0,"Nshocks_max"));

    par->Nshocks  = (int*) mxGetData(mxGetField(prhs[in],0,"Nshocks"));
    par->psi      = (double*) mxGetPr(mxGetField(prhs[in],0,"psi"));
    par->xi       = (double*) mxGetPr(mxGetField(prhs[in],0,"xi"));      
    par->weights  = (double*) mxGetPr(mxGetField(prhs[in],0,"weights"));

    // g. assets
    par->R    = (double) mxGetScalar(mxGetField(prhs[in],0,"R"));
    par->Rneg = (double) mxGetScalar(mxGetField(prhs[in],0,"Rneg"));
    par->Rb   = (double) mxGetScalar(mxGetField(prhs[in],0,"Rb"));

    par->omega      = (double*) mxGetPr(mxGetField(prhs[in],0,"omega"));
    par->adjcost    = (double) mxGetScalar(mxGetField(prhs[in],0,"adjcost"));
    par->adjcostinf = (double) mxGetScalar(mxGetField(prhs[in],0,"adjcostinf"));

    // h. inheritance
    par->inh_A   = (double) mxGetScalar(mxGetField(prhs[in],0,"inh_A"));
    par->inh_fac = (double*) mxGetPr(mxGetField(prhs[in],0,"inh_fac"));      
    par->inh_p   = (double*) mxGetPr(mxGetField(prhs[in],0,"inh_p"));          

    // i. simualate
    if(type == 2){

        // i. time and number of households
        par->T0sim        = (int) mxGetScalar(mxGetField(prhs[in],0,"T0sim"));
        par->Tendsim      = (int) mxGetScalar(mxGetField(prhs[in],0,"Tendsim"));    
        par->Nsim         = (int) mxGetScalar(mxGetField(prhs[in],0,"Nsim"));

        // ii. settings
        par->forget_inh   = (int) mxGetScalar(mxGetField(prhs[in],0,"forget_inh"));          
        par->save_to_disc = (int) mxGetScalar(mxGetField(prhs[in],0,"save_to_disc"));
        par->sim_IRF      = (int) mxGetScalar(mxGetField(prhs[in],0,"sim_IRF"));
        par->MPC_add      = (double) mxGetScalar(mxGetField(prhs[in],0,"MPC_add"));

        // iii. shocks
        par->psi_sim      = (double*) mxGetPr(mxGetField(prhs[in],0,"psi_sim"));
        par->xi_sim       = (double*) mxGetPr(mxGetField(prhs[in],0,"xi_sim"));
        par->A_shock      = (double*) mxGetPr(mxGetField(prhs[in],0,"A_shock"));
        par->inh_pval     = (double*) mxGetPr(mxGetField(prhs[in],0,"inh_pval"));
        par->adjcost_pval = (double*) mxGetPr(mxGetField(prhs[in],0,"adjcost_pval"));        
        par->P_IRF        = (double*) mxGetPr(mxGetField(prhs[in],0,"P_IRF"));
        par->A_IRF        = (double*) mxGetPr(mxGetField(prhs[in],0,"A_IRF"));
        par->B_IRF        = (double*) mxGetPr(mxGetField(prhs[in],0,"B_IRF"));

        // iv. initial distribution
        par->P_ini_dist      = (double*) mxGetPr(mxGetField(prhs[in],0,"P_ini_dist"));
        par->N_ini_dist      = (double*) mxGetPr(mxGetField(prhs[in],0,"N_ini_dist"));

    }

    if(type == 3){

        par->cutoff  = (double) mxGetScalar(mxGetField(prhs[in],0,"cutoff"));
        par->euler_t = (int) mxGetScalar(mxGetField(prhs[in],0,"euler_t"));
        par->Neuler  = (int) mxGetScalar(mxGetField(prhs[in],0,"Neuler"));

        par->euler_M = (double*) mxGetPr(mxGetField(prhs[in],0,"euler_M"));
        par->euler_N = (double*) mxGetPr(mxGetField(prhs[in],0,"euler_N"));
        par->euler_P = (double*) mxGetPr(mxGetField(prhs[in],0,"euler_P"));

    }

    // j. technical
    par->tmin       = (int) mxGetScalar(mxGetField(prhs[in],0,"tmin"));
    par->vfi_pure   = (int) mxGetScalar(mxGetField(prhs[in],0,"vfi_pure")); 
    par->no_multistart = (int) mxGetScalar(mxGetField(prhs[in],0,"no_multistart")); 
    par->no_interp_vec = (int) mxGetScalar(mxGetField(prhs[in],0,"no_interp_vec")); 
        
        in++;

        // last period
        par->Rret = MAX(par->R,par->Rb);
        par->gamma_0 = (1.0-pow(par->Rret,-(par->T-par->TR))) / (1.0-pow(par->Rret,-1)) - 1.0;
        
            double fac = pow(par->Rret,-1)*pow(par->Rret*par->beta,1.0/par->sigma);
            double nom = 1.0-fac;
            double denom = 1.0-pow(fac,par->T-par->TR);
           
        par->gamma_1 = nom/denom;
           
            fac = pow(par->beta,1.0/par->sigma)*pow(par->Rret,(1.0-par->sigma)/par->sigma);
            nom = 1.0-pow(fac,par->T-par->TR);
            denom = 1.0-fac;
           
        par->gamma_2 = nom/denom;

    // solutions
    if(type == 2 || type == 3){ // when simulating

        par->C_ast = new double*[par->T*par->Nd*par->Nz];
        par->B_ast = new double*[par->T*par->Nd*par->Nz];
        par->A_ast = new double*[par->T*par->Nd*par->Nz];        
        par->invV  = new double*[par->T*par->Nd*par->Nz];

        for(int t = 0; t < par->TR; t++){
        for(int d = 0; d < par->Nd; d++){
        for(int z = 0; z < par->Nz; z++){                    

            int i_cell = index::d3(t,d,z,par->Nd,par->Nz);
            par->C_ast[i_cell] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[in],0,"C_ast"),i_cell));
            par->B_ast[i_cell] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[in],0,"B_ast"),i_cell));
            par->A_ast[i_cell] = (double*) mxGetPr(mxGetCell(mxGetField(prhs[in],0,"A_ast"),i_cell));                         
            par->invV[i_cell]  = (double*) mxGetPr(mxGetCell(mxGetField(prhs[in],0,"invV"),i_cell));

        } } }

    }


    ////////////////////////
    // 2. outputs - solve //
    ////////////////////////

    if(type == 1){

        // a. struct
        const char *field_names[] = {"C_ast", "B_ast", "invV", 
                                     "V", "W", "invW", "Q", "A", 
                                     "C_EGM", "M_EGM",};
        
        int num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sim_struct = plhs[0];

        // b. cell dimensions
        int ndim_cell  = 3;
        auto dims_cell = new size_t[3];

        // c. array dimensions
        auto ndim = new size_t[2];
        auto dims = new size_t*[2];
                        
        // d. solution

            // cell
            dims_cell[0] = par->Nz;
            dims_cell[1] = par->Nd;
            dims_cell[2] = par->TR;

            // array      
            ndim[0] = 3;
            dims[0] = new size_t[3];           
            dims[0][0] = par->NM;
            dims[0][1] = par->NN;
            dims[0][2] = par->NP;

            ndim[1] = 2;
            dims[1] = new size_t[2];
            dims[1][0] = par->NX;
            dims[1][1] = par->NP;    

        par->C_ast    = mymex::set_field_cell(sim_struct,"C_ast",ndim_cell,dims_cell,ndim,dims); 
        par->B_ast    = mymex::set_field_cell(sim_struct,"B_ast",ndim_cell,dims_cell,ndim,dims); 
        par->invV     = mymex::set_field_cell(sim_struct,"invV",ndim_cell,dims_cell,ndim,dims);                                
        par->V        = mymex::set_field_cell(sim_struct,"V",ndim_cell,dims_cell,ndim,dims);

        // e. post-decision
            
            // cell
            dims_cell[0] = 1; // z is not a post-decision state
            dims_cell[1] = par->Nd;            
            dims_cell[2] = 1; // not saved over T
            

            // array            
            dims[0][0] = par->NA;
            dims[0][1] = par->NN;
            dims[0][2] = par->NP;
            
        par->W      = mymex::set_field_cell(sim_struct,"W",ndim_cell,dims_cell,ndim,dims);
        par->invW   = mymex::set_field_cell(sim_struct,"invW",ndim_cell,dims_cell,ndim,dims); 
        par->Q      = mymex::set_field_cell(sim_struct,"Q",ndim_cell,dims_cell,ndim,dims);
        par->C_EGM  = mymex::set_field_cell(sim_struct,"C_EGM",ndim_cell,dims_cell,ndim,dims); 
        par->M_EGM  = mymex::set_field_cell(sim_struct,"M_EGM",ndim_cell,dims_cell,ndim,dims); 

            delete[] dims_cell;
            delete[] ndim;
            delete[] dims[0];
            delete[] dims[1];
            delete[] dims;

    }


    ///////////////////////////
    // 3. outputs - simulate //
    ///////////////////////////

    if(type == 2){

        // a. struct
        const char *field_names[] = {"d", "z", "M", "P", "Y", "N", "X", "H", "C", "B", "A", "MPC", "simV"};
        
        int num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sim_struct = plhs[0]; 

        // b. dimensions
        int ndim = 2;
        size_t* dims = new size_t[2];

            dims[0] = par->Nsim;
            dims[1] = par->Tendsim;

        // c. elements
        par->d = mymex::set_field_int(sim_struct,"d",ndim,dims);
        par->z = mymex::set_field_int(sim_struct,"z",ndim,dims);
        par->M = mymex::set_field_double(sim_struct,"M",ndim,dims);
        par->N = mymex::set_field_double(sim_struct,"N",ndim,dims);
        par->P = mymex::set_field_double(sim_struct,"P",ndim,dims);
        par->Y = mymex::set_field_double(sim_struct,"Y",ndim,dims);        
        par->X = mymex::set_field_double(sim_struct,"X",ndim,dims);
        par->H = mymex::set_field_double(sim_struct,"H",ndim,dims);
        par->C = mymex::set_field_double(sim_struct,"C",ndim,dims);
        par->B = mymex::set_field_double(sim_struct,"B",ndim,dims);
        par->A = mymex::set_field_double(sim_struct,"A",ndim,dims);
        par->MPC = mymex::set_field_double(sim_struct,"MPC",ndim,dims);
        par->simV = mymex::set_field_double(sim_struct,"simV",ndim,dims);

            delete[] dims;

    }
 

    /////////////////////////
    // 4. outputs - euler  //
    /////////////////////////

    if(type == 3){

        // a. struct
        const char *field_names[] = {"euler_C", "euler_error"};
                
        int num_fields = sizeof(field_names)/sizeof(*field_names);
        plhs[0] = mxCreateStructMatrix(1, 1, num_fields, field_names);
        auto sim_struct = plhs[0]; 

        // b. dimensions
        int ndim = 2;
        size_t* dims = new size_t[2];

        dims[0] = par->Neuler;
        dims[1] = 1;

        // c. elements
        par->euler_C     = mymex::set_field_double(sim_struct,"euler_C",ndim,dims);
        par->euler_error = mymex::set_field_double(sim_struct,"euler_error",ndim,dims);   

        delete[] dims;


    }

}


////////////////
// 3. destroy //
////////////////

void destroy(par_struct *par, int type){

    if(type == 1){ // solving

        delete[] par->C_ast;
        delete[] par->B_ast;
        delete[] par->invV;
        delete[] par->V;

        delete[] par->W;
        delete[] par->invW;
        delete[] par->Q;
        delete[] par->C_EGM;
        delete[] par->M_EGM;

    } else { // simulating

        delete[] par->C_ast;
        delete[] par->B_ast;
        delete[] par->A_ast;        
        delete[] par->invV;

    }

}

} // namespace