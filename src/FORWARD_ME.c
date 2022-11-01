
#include "FORWARD_ME.h"

/*--------------------------------------------------------------------------------*/

extern int Forward_ME(STRUCT_ATOM *Atom, double *Par, STRUCT_OBSERVATION \
    *Observation){
  
    /*######################################################################
    Purpose:
      Calculate the Stokes profiles under M-E atmosphere model (normal Zeeman effect)
    Record of revisions:
      30 Otc. 2022
    Input parameters:
      Atom, a structure saved the Landu facor.
      Par, an array saved the model parameters.
      Observation, a structure saved the lambda.
    Output parameters:
      Stokes[4][] in Observation, the synthesized Stokes profiles.
      L, associated dispersion profile.
    Reference:
      LL04 5.47, LL04 9.4 9.8, del Toro Iniesta 2003
     ######################################################################*/
  
  /*
   double Bmag, Theta, Chi; 
   Bmag=sqrt(Bz*Bz+Bx*Bx+By*By);
   Theta=acos(Bz/Bmag);
   Chi=atan(By/Bx);
   */

    double Bmag, Btheta, Bchi;
    Bmag = Par[0];
    Btheta = Par[1];
    Bchi = Par[2];
    //Bmag = sqrt(Par[0]*Par[0]+Par[1]*Par[1]);
    //Btheta = acos(Par[0]/Bmag);
    //Atom->BShift = 4.6686e-3*Atom->Geffect*Atom->Lambda*Atom->Lambda;
        // Atom->Lambda in A, Par[4] in mA, Atom->BShift in mA
    Atom->BShift = 4.6686e-10*Atom->Geffect*Atom->Lambda*Atom->Lambda;

    double *phi, *psi;
    phi = (double *)VECTOR(-1, 1, enum_dbl, false);
    psi = (double *)VECTOR(-1, 1, enum_dbl, false);
  
    double B0, B1;
    B0 = Par[7]*Par[8];
    B1 = Par[7]-B0;
    
    int Indx_Lamd, q;
    double Delta = 0, eta[4], rho[4], temp, temp1, temp2, temp3;

    double Lam_Shifts, ShiftW;
    double H, L;
   
    double ShiftB = Bmag*Atom->BShift/Par[4];
    double ShiftV = 1e3*(-1e3*Par[3]/Par_C*Atom->Lambda-Atom->Lambda)/Par[4];
    
    double SIN_THB = sin(Btheta);
    double COS_THB = cos(Btheta);
    double SIN_THB_SQ = SIN_THB*SIN_THB;
    double COS_THB_SQ = COS_THB*COS_THB;
   
    double COS2PHI = cos(2.*Bchi);
    double SIN2PHI = sin(2.*Bchi);


    for (Indx_Lamd = 0; Indx_Lamd<Observation->Num_Lam; Indx_Lamd++){

      ShiftW = 1e3*Observation->Lambda[Indx_Lamd]/Par[4]+ShiftV;

      for (q=-1; q<=1; q++){
        // q=Mu-Ml,(1 sigma bule, -1 sigma red, 0 pi)
        Lam_Shifts = ShiftW+q*ShiftB;
        Faddeeva(Lam_Shifts, Par[5], &H, &L, 6);
        phi[q] = H/Par_SqrtPi;
        psi[q] = L/Par_SqrtPi;
      }

      //Stokes profiles LL04 Page 414 Eq 9.110
      eta[0] = 0.5*Par[6]*(phi[0]*SIN_THB_SQ \
          +0.5*(phi[-1]+phi[1])*(1+COS_THB_SQ));
      eta[1] = 0.5*Par[6]*(phi[0]-0.5*(phi[-1]+phi[1])) \
          *SIN_THB_SQ*COS2PHI;
      eta[2] = 0.5*Par[6]*(phi[0]-0.5*(phi[-1]+phi[1])) \
          *SIN_THB_SQ*SIN2PHI;
      eta[3] = 0.5*Par[6]*(phi[-1]-phi[1])*COS_THB;
      rho[1] = 0.5*Par[6]*(psi[0]-0.5*(psi[-1]+psi[1])) \
          *SIN_THB_SQ*COS2PHI;
      rho[2] = 0.5*Par[6]*(psi[0]-0.5*(psi[-1]+psi[1])) \
          *SIN_THB_SQ*SIN2PHI;
      rho[3] = 0.5*Par[6]*(psi[-1]-psi[1])*COS_THB;

      temp = 1+eta[0];
      temp1 = eta[1]*eta[1]+eta[2]*eta[2]+eta[3]*eta[3];
      temp2 = rho[1]*rho[1]+rho[2]*rho[2]+rho[3]*rho[3];
      temp3 = eta[1]*rho[1]+eta[2]*rho[2]+eta[3]*rho[3];
      Delta = (temp*temp*(temp*temp-temp1+temp2)-temp3*temp3);

      Observation->Stokes[0][Indx_Lamd] = B0+B1/Delta*(temp*(temp*temp+temp2));
      Observation->Stokes[1][Indx_Lamd] = -B1/Delta*(temp*temp*eta[1] \
          +temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3);
      Observation->Stokes[2][Indx_Lamd] = -B1/Delta*(temp*temp*eta[2] \
          +temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3);
      Observation->Stokes[3][Indx_Lamd] = -B1/Delta*(temp*temp*eta[3] \
          +temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3);

    }
    FREE_VECTOR(phi, -1, enum_dbl);
    FREE_VECTOR(psi, -1, enum_dbl);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/
