namespace utilities {

///////////////////////
// 1. transformation //
///////////////////////

double trans(double x, par_struct *par) {

    if(par->epstein_zin == 0 && par->rho >= 1.0){ // v in ]-inf 0[
        if(par->rho == 1){
            return x;
        }
        if(x == -HUGE_VAL){
            return 0.0;
        } else {
            return -1.0/x; // negative inverse
        }        
    } else { // v in [0 inf[
        if(x == -HUGE_VAL){ // -HUGE_VAL indicate minimumal value
            return 0.0;
        } else {
            return x; // no transformation
        }
    }

}

double inv_trans(double x, par_struct *par){
    
    if(par->epstein_zin == 0 && par->rho >= 1.0){   
        if(par->rho == 1){
            return x;
        }     
        if(x <= 0.0){
            return -HUGE_VAL;
        } else {
            return -1.0/x; // negative inverse
        }
    } else {
        return x; // no transformation
    }

}


////////////////
// 2. utility //
////////////////

double u(double C, par_struct *par){

    if(par->epstein_zin == 0){
        return pow(C,1.0-par->rho)/(1.0-par->rho);
    } else {
        return (1.0-par->beta)*pow(C,1.0-par->sigma);
    }

}

double marg_u(double C, par_struct *par){

    if(par->epstein_zin == 0){    
        return pow(C,-par->rho);
    } else {
        return pow(C,-par->sigma);        
    }

}

double inv_marg_u(double u, par_struct *par){
    
    if(par->epstein_zin == 0){      
        return pow(u,-1.0/par->rho);
    } else {
        return pow(u,-1.0/par->sigma);
    }

}

}