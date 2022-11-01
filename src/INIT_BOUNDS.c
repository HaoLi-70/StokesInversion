
#include "INIT_BOUNDS.h"

/*--------------------------------------------------------------------------------*/

extern int Init_Bounds(STRUCT_PAR *Str_Par, STRUCT_OBSERVATION *Observation, \
    STRUCT_ATOM *Atom){
    
    /*######################################################################
      Purpose:
        allocate a matrix m[i][j] with subscript range nrl<=i<=nrh,
            ncl<=j<=nch.
      Record of revisions:
        8 Sept. 2021.
      Input parameters:
        Atom, a structure saved the Landu facor.
        Str_Par, a structure saved the model parameters, bounds.
        Observation, a structure saved the lambda.
      Output parameters:
        Str_Par
      Return:
        return the pointer.
     ######################################################################*/
    
  
    int i;
    double MAX=Observation->Profile_Data[0][0];
    double sum=0, weight=0;

    Str_Par->Num_Par = 9;
    Str_Par->Par_Bounds = (double **)MATRIX(0, Str_Par->Num_Par-1, \
        0, 1, enum_dbl, false);
    Str_Par->Distribution = (double *)VECTOR(0, Str_Par->Num_Par-1, \
        enum_dbl, false);
    Str_Par->type = (enum bounds_type *)malloc(Str_Par->Num_Par \
        *sizeof(enum bounds_type));
    Str_Par->narrower_guess = (bool *)malloc(Str_Par->Num_Par*sizeof(bool));

    //S0+S1: MAX value in Stokes I
    for (i=1; i<Observation->Num_Lam; i++) {
      if (Observation->Profile_Data[0][i]>MAX)        
        MAX = Observation->Profile_Data[0][i];
    }
    
    Str_Par->Par_Bounds[7][0] = MAX*0.7;
    Str_Par->Par_Bounds[7][1] = MAX*1.1;
    
    //Velocity: we use the COG method in Stokes I
    for (i=0; i<Observation->Num_Lam; i++) {
      sum += Observation->Lambda[i]*(MAX-Observation->Profile_Data[0][i]);
      weight += MAX-Observation->Profile_Data[0][i];
    }
    
    double Delta_Velocity;
    //=(Observation->Profile_Data[1][0]-Observation->Profile_Data[0][0])/Atom->Lambda*Par_C;
    Delta_Velocity = 1.026*3;
    double Velocity = (sum/weight-Atom->Lambda)/Atom->Lambda*Par_C/1e3;
       

    Str_Par->Par_Bounds[3][0] = Velocity-Delta_Velocity;
    Str_Par->Par_Bounds[3][1] = Velocity+Delta_Velocity;

    Str_Par->Par_Bounds[0][0] = 0;
    Str_Par->Par_Bounds[0][1] = 4000;

    Str_Par->Par_Bounds[1][0] = 0;
    Str_Par->Par_Bounds[1][1] = 3.1415926;

    //Str_Par->Par_Bounds[1][0]=0;
    //Str_Par->Par_Bounds[1][1]=4000;

    Str_Par->Par_Bounds[2][0] = 0;
    Str_Par->Par_Bounds[2][1] = 3.1415926;
    
    Str_Par->Par_Bounds[4][0] = 15.0;
    Str_Par->Par_Bounds[4][1] = 55.0;
    
    Str_Par->Par_Bounds[5][0] = 0;
    Str_Par->Par_Bounds[5][1] = 0.8;

    Str_Par->Par_Bounds[6][0] = 0;
    Str_Par->Par_Bounds[6][1] = 70;//90
    
    Str_Par->Par_Bounds[8][0] = 0;
    Str_Par->Par_Bounds[8][1] = 1;//1.2

    for(i=0; i<Str_Par->Num_Par; i++){
      Str_Par->Distribution[i] = 1e-3 \
          *(Str_Par->Par_Bounds[i][1]-Str_Par->Par_Bounds[i][0]);
      Str_Par->type[i] = enum_reflect;
      Str_Par->narrower_guess[i] = false;
    }

    Str_Par->narrower_guess[3] = true;
    Str_Par->narrower_guess[7] = true;
    Str_Par->type[2] = enum_fold;

    for(i =0 ; i<Observation->Num_Lam; i++){
      Observation->Sigma[0][i] = MAX*1e-3;
      Observation->Sigma[1][i] = MAX*1e-3;
      Observation->Sigma[2][i] = MAX*1e-3;
      Observation->Sigma[3][i] = MAX*1e-3;
    }
   
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

