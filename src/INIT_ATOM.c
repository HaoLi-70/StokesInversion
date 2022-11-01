
#include "INIT_ATOM.h"

/*--------------------------------------------------------------------------------*/

extern double Geffect(double Gu, double Gl, double Ju, double Jl){
    
    /*######################################################################
      Purpose:
        Computes the effective Lande factor.
      Record of revisions:
        25 Apr. 2018
      Input parameters:
        Gu, The lande factor of upper level.
        Gl, The lande factor of lower level.
        Ju, The total angular momentum of upper level.
        Jl, The total angular momentum of lower level.
      Output parameters:
        Geffect, The effect Lande factor.
      References:
        LL04 Chapter 3, Equation 3.44.
     ######################################################################*/
    
    double Geffect;
    
    Geffect = 0.5*(Gu+Gl)+0.25*(Gu-Gl)*(Ju*(Ju+1)-Jl*(Jl+1));
    
    return Geffect;
}

/*--------------------------------------------------------------------------------*/

extern double Gfactor(double J, double L, double S){
    
    /*######################################################################
      Purpose:
        Computes the Lande factor.
      Record of revisions:
        25 Apr. 2018
      Input parameters:
        J, The total angular momentum of the electronic cloud.
        L, The total orbital angular momentum of the electronic cloud.
        S, The total spin of the electronic cloud.
      Output parameters:
        Gfactor, The Lande factor.
      Method:
        L-S coupling asumption, and if the number is 0, return 0 
          for convernience.
     ######################################################################*/
    
    
    if(J==0) return 0.0;

    double g = 1.0+(J*(J+1.0)-L*(L+1.0)+S*(S+1.0))/(2.0*J*(J+1.0));
    
    return g;
    
}

/*--------------------------------------------------------------------------------*/

extern int Init_Atom(STRUCT_ATOM *Atom){
    
    /*######################################################################
      Purpose:
        Initial the atomic information.
      Record of revisions:
        30 Otc. 2022.
      Output parameters:
        Atom, a structure the atomic information
     ######################################################################*/
    
    Atom->Lambda = 6302.4932;

    Atom->Ju = 0.0;
    Atom->Lu = 2.0;
    Atom->Su = 2.0;
    
    Atom->Jl = 1.0;
    Atom->Ll = 1.0;
    Atom->Sl = 2.0;
    
    Atom->Gu = Gfactor(Atom->Ju, Atom->Lu, Atom->Su);
    Atom->Gl = Gfactor(Atom->Jl, Atom->Ll, Atom->Sl);
    Atom->Geffect = Geffect(Atom->Gu, Atom->Gl, Atom->Ju, Atom->Jl);
    
    fprintf(stderr, "geffect %e \n", Atom->Geffect);
    
    return 0;    
}

/*--------------------------------------------------------------------------------*/

