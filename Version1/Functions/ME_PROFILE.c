#include "ME_PROFILE.h"


void Stokes_profile(double Bmag, double Theta, double Chi, double Velocity, double Doppler_width, double Source_function, double Coeffi_Abs, double Beta, double a, double Lambda, TRANSITION_INFORMATION Atom, double *Stokes){
    
    /***********************************************************************************
     Purpose:
     Calculate the Stokes profiles under M-E atmosphere model (hyperfine structure)
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     Bmag, strenth of magnetic field vector.
     Theta, inclination of the magnetic field vector.
     Chi, azimuth of the magnetic field.
     Velocity, line of sight velocity.
     Doppler_width, Doppler width in wavelength (meter).
     Source function, the sum of source function $S_0$ and gradient of the source function. S0+S1 in the M-E atmosphere.
     Coeffi_Abs, ratio between the line and continuous absorption coefficients.
     a, the damping parameter a (in the unit of Doppler width)
     Lambda, the wavelength to calculate.
     Atom, a sturcture which saves Atomic information.
     Output:
     Stokes[4], the Stokes parameters i=0,1,2,3  I,Q,U,V respectively
     Reference:
     LL04 5.47, LL04 9.4 9.8, del Toro Iniesta 2003
     ***********************************************************************************/

    //Speed of light
    double Con_c=299792458.0;
    
    //Ratio of circumference to diameter  3.1415926535897932384626433832795028841971
    double Con_Pi=3.141592653589793;
    
    double *phi, *psi;
    
    phi=VECTOR_DOUBLE(-1, 1);
    
    psi=VECTOR_DOUBLE(-1, 1);
    
    double *H, *L;
    
    H=(double *)malloc((sizeof(double)));
    
    L=(double *)malloc((sizeof(double)));
    
    
    float Mu, Ml;
    
    int q;
    
    double Deltanu=0, Delta=0, eta[4], rho[4];
    
    double temp, temp1, temp2, temp3, B0, B1;
    
    B0=Source_function*Beta;
    
    B1=Source_function-B0;
 
    for (q=-1; q<=1; q++) {
        // q=Mu-Ml,(1 sigma bule, -1 sigma red, 0 pi)

        phi[q]=0;
        
        psi[q]=0;
        
        for (Mu=Atom.Ju; Mu>=-Atom.Ju; Mu--) {
            
            Ml=Mu+q;
            
            if (Ml<=Atom.Jl&&Ml>=-Atom.Jl) {
                
                //3j symble
                temp = TJ(Atom.Ju, Atom.Jl, 1, -Mu, Ml, -q);
                
                //LL04 3.14 5.47
                if (temp!=0) {
                    
                    Deltanu=(Lambda-Velocity/Con_c*Atom.Lambda-Atom.Lambda+4.6686e-3*Bmag*(Atom.Gu*Mu-Atom.Gl*Ml)*Atom.Lambda*Atom.Lambda)/Doppler_width;
                    
                    Faddeeva_Function(Deltanu, a, H, L);
                    
                    //LL04 9.8
                    temp = 3*temp*temp/sqrt(Con_Pi);
                    
                    phi[q] += temp*(*H);
                    
                    psi[q] += temp*(*L);
                    
                }
            }
        }
    }
    
    //coefficients of the absorption matrix LL04 Page 386 Eq 9.27
    eta[0] = 0.5*Coeffi_Abs*(phi[0]*sin(Theta)*sin(Theta)+(phi[-1]+phi[1])*0.5*(1+cos(Theta)*cos(Theta)));
    
    eta[1] = 0.5*Coeffi_Abs*(phi[0]-(phi[-1]+phi[1])/2.)*sin(Theta)*sin(Theta)*cos(2.*Chi);
    
    eta[2] = 0.5*Coeffi_Abs*(phi[0]-(phi[-1]+phi[1])/2.)*sin(Theta)*sin(Theta)*sin(2.*Chi);
    
    eta[3] = 0.5*Coeffi_Abs*(phi[1]-phi[-1])*cos(Theta);
    
    rho[1] = 0.5*Coeffi_Abs*(psi[0]-(psi[-1]+psi[1])/2.)*sin(Theta)*sin(Theta)*cos(2.*Chi);
    
    rho[2] = 0.5*Coeffi_Abs*(psi[0]-(psi[-1]+psi[1])/2.)*sin(Theta)*sin(Theta)*sin(2.*Chi);
    
    rho[3] = 0.5*Coeffi_Abs*(psi[1]-psi[-1])*cos(Theta);
    
    //Stokes profiles LL04 Page 414 Eq 9.110
    temp=1+eta[0];
    
    temp1=eta[1]*eta[1]+eta[2]*eta[2]+eta[3]*eta[3];
    
    temp2=rho[1]*rho[1]+rho[2]*rho[2]+rho[3]*rho[3];
    
    temp3=eta[1]*rho[1]+eta[2]*rho[2]+eta[3]*rho[3];
    
    Delta=temp*temp*(temp*temp-temp1+temp2)-temp3*temp3;
    
    Stokes[0]=B0+B1/Delta*(temp*(temp*temp+temp2));
    
    Stokes[1]=-B1/Delta*(temp*temp*eta[1]+temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3);
    
    Stokes[2]=-B1/Delta*(temp*temp*eta[2]+temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3);
    
    Stokes[3]=-B1/Delta*(temp*temp*eta[3]+temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3);
    
    free(H);
    
    free(L);
    
    FREE_VECTOR_DOUBLE(phi, -1);
    
    FREE_VECTOR_DOUBLE(psi, -1);
    
    return;
    
}


void Stokes_profile_Normal(double Bmag, double Theta, double Chi, double Velocity, double Doppler_width, double Source_function, double Coeffi_Abs, double Beta, double a, double Lambda, TRANSITION_INFORMATION Atom, double *Stokes){
    
    /***********************************************************************************
     Purpose:
     Calculate the Stokes profiles under M-E atmosphere model (normal Zeeman effect)
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     Bmag, strenth of magnetic field vector.
     Theta, inclination of the magnetic field vector.
     Chi, azimuth of the magnetic field.
     Velocity, line of sight velocity.
     Doppler_width, Doppler width in wavelength (meter).
     Source function, the sum of source function $S_0$ and gradient of the source function. S0+S1 in the M-E atmosphere.
     Coeffi_Abs, ratio between the line and continuous absorption coefficients.
     a, the damping parameter a (in the unit of Doppler width)
     Lambda, the wavelength to calculate.
     Atom, a sturcture which saves Atomic information.
     Output:
     Stokes[4], the Stokes parameters i=0,1,2,3  I,Q,U,V respectively
     Reference:
     LL04 5.47, LL04 9.4 9.8, del Toro Iniesta 2003
     ***********************************************************************************/

    //Speed of light
    double Con_c=299792458.0;
    
    //Ratio of circumference to diameter 3.1415926535897932384626433832795028841971
    double Con_Pi=3.141592653589793;
   
    double *phi, *psi;
    
    phi=VECTOR_DOUBLE(-1, 1);
    
    psi=VECTOR_DOUBLE(-1, 1);
    
    double *H, *L;
    
    H=(double *)malloc((sizeof(double)));
    
    L=(double *)malloc((sizeof(double)));
    
    double Deltanu=0, Delta=0, eta[4], rho[4];
    
    int q;
    
    double Lambda0, LambdaB, temp=1, temp1, temp2, temp3, B0, B1;
    
    B0=Source_function*Beta;
    
    B1=Source_function-B0;
    
    Lambda0=(Lambda-Velocity/Con_c*Atom.Lambda-Atom.Lambda)/Doppler_width;
    
    LambdaB=4.6686e-3*Bmag*Atom.Geffect*Atom.Lambda*Atom.Lambda/Doppler_width;
    
    for (q=-1; q<=1; q++) {
        // q=Mu-Ml,(1 sigma bule, -1 sigma red, 0 pi)
        
        Deltanu=Lambda0+q*LambdaB;
        
        Faddeeva_Function(Deltanu, a, H, L);
        
        phi[q] = (*H)/sqrt(Con_Pi);
        
        psi[q] = (*L)/sqrt(Con_Pi);
        
    }
    
    //coefficients of the absorption matrix (normal Zeeman effect)
    eta[0] = 0.5*Coeffi_Abs*(phi[0]*sin(Theta)*sin(Theta)+(phi[-1]+phi[1])/2.*(1+cos(Theta)*cos(Theta)));
    
    eta[1] = 0.5*Coeffi_Abs*(phi[0]-(phi[-1]+phi[1])/2.)*sin(Theta)*sin(Theta)*cos(2.*Chi);
    
    eta[2] = 0.5*Coeffi_Abs*(phi[0]-(phi[-1]+phi[1])/2.)*sin(Theta)*sin(Theta)*sin(2.*Chi);
    
    eta[3] = 0.5*Coeffi_Abs*(phi[-1]-phi[1])*cos(Theta);
    
    rho[1] = 0.5*Coeffi_Abs*(psi[0]-(psi[-1]+psi[1])/2.)*sin(Theta)*sin(Theta)*cos(2.*Chi);
    
    rho[2] = 0.5*Coeffi_Abs*(psi[0]-(psi[-1]+psi[1])/2.)*sin(Theta)*sin(Theta)*sin(2.*Chi);
    
    rho[3] = 0.5*Coeffi_Abs*(psi[-1]-psi[1])*cos(Theta);
    
    //Stokes profiles LL04 Page 414 Eq 9.110
    temp=1+eta[0];
    
    temp1=eta[1]*eta[1]+eta[2]*eta[2]+eta[3]*eta[3];
    
    temp2=rho[1]*rho[1]+rho[2]*rho[2]+rho[3]*rho[3];
    
    temp3=eta[1]*rho[1]+eta[2]*rho[2]+eta[3]*rho[3];
    
    Delta=(temp*temp*(temp*temp-temp1+temp2)-temp3*temp3);
    
    Stokes[0]=B0+B1/Delta*(temp*(temp*temp+temp2));
    
    Stokes[1]=-B1/Delta*(temp*temp*eta[1]+temp*(eta[3]*rho[2]-eta[2]*rho[3])+rho[1]*temp3);
    
    Stokes[2]=-B1/Delta*(temp*temp*eta[2]+temp*(eta[1]*rho[3]-eta[3]*rho[1])+rho[2]*temp3);
    
    Stokes[3]=-B1/Delta*(temp*temp*eta[3]+temp*(eta[2]*rho[1]-eta[1]*rho[2])+rho[3]*temp3);
    
    free(H);
    
    free(L);
    
    FREE_VECTOR_DOUBLE(phi, -1);
    
    FREE_VECTOR_DOUBLE(psi, -1);
    
    return;
    
}
