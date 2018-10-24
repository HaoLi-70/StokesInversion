#include "ANGULAR_MOMENTUM.h"


extern int Factorial(int N){
    
    /***********************************************************************************
     Purpose:
     Computes the factorial N!.
     Modified:
     25 Apr. 2018, Hao Li
     Input parameters:
     N, The int number.
     Output parameters:
     F, The farcorial.
     ***********************************************************************************/
    
    int i,F=1;
    
    for (i=1; i<=N; i++)
        
        F=F*i;
    
    return F;
    
}


extern double Geffect(float Gu, float Gl, float Ju, float Jl){
    
    /***********************************************************************************
     Purpose:
     Computes the effective Lande factor.
     Modified:
     25 Apr. 2018, Hao Li
     Input parameters:
     Gu, The lande factor of upper level.
     Gl, The lande factor of lower level.
     Ju, The total angular momentum of upper level.
     Jl, The total angular momentum of lower level.
     Output parameters:
     Geffect, The effect Lande factor.
     References:
     LL04 Chapter 3, Equation 3.44.
     ***********************************************************************************/
    
    double Geffect;
    
    Geffect=0.5*(Gu+Gl)+0.25*(Gu-Gl)*(Ju*(Ju+1)-Jl*(Jl+1));
    
    return Geffect;
    
}


extern double Gfactor(float J,float L,float S){
    
    /***********************************************************************************
     Purpose:
     Computes the Lande factor.
     Modified:
     25 Apr. 2018, Hao Li
     Input parameters:
     J, The total angular momentum of the electronic cloud.
     L, The total orbital angular momentum of the electronic cloud.
     S, The total spin of the electronic cloud.
     Output parameters:
     Gfactor, The Lande factor.
     Method:
     L-S coupling asumption, and if the number is 0, return 0 for convernience.
     ***********************************************************************************/
    
    double g;
    
    if (J==0)
        
        return 0;

    g=1.0+(J*(J+1.0)-L*(L+1.0)+S*(S+1.0))/(2.0*J*(J+1.0));
    
    return g;
    
}


extern double TJ(float J1,float J2,float J3,float M1,float M2,float M3){
    
    /***********************************************************************************
     Purpose:
     Computes the 3j symbol.
     Modified:
     25 Apr. 2018, Hao Li
     Input parameters:
     J1, J2, J3, The total angular momentum of the electronic cloud.
     M1, M2, M3, The magnetic quantum number.
     Output parameters:
     F, The 3j symbol.
     References:
     LL04 Chapter 2, Page 36, Equation 2.19 and Page 38, Equation 2.22.
     ***********************************************************************************/
    
    if((M1+M2+M3)!=0)
        
        return 0;
    
    if((J1+J2)<J3||fabsf(J1-J2)>J3)
        
        return 0;
    
    if(fabsf(M1)>J1)
        
        return 0;
    
    if(fabsf(M2)>J2)
        
        return 0;
    
    if(fabsf(M3)>J3)
        
        return 0;
    
    if(fmod((J1+J2+J3)*2.0,2.0)!=0.0)
        
        return 0;
    
    if(fmod((J1-M1)*2.0,2.0)!=0.0)
        
        return 0;
    
    if(fmod((J2-M2)*2.0,2.0)!=0.0)
        
        return 0;
    
    double a,c,d,e,F;
    
    float t,t1,t2;
    
    e=0;
    
    a=pow(-1,J1-J2-M3);
    
    c=sqrt(1.0*Factorial(J1+J2-J3)*Factorial(J1-J2+J3)*Factorial(-J1+J2+J3)/Factorial(J1+J2+J3+1));
    
    d=sqrt(1.0*Factorial(J1+M1)*Factorial(J1-M1)*Factorial(J2+M2)*Factorial(J2-M2)*Factorial(J3+M3)*Factorial(J3-M3));
    
    t1=0;
    
    if((J2-J3-M1)>t1)
        
        t1=J2-J3-M1;
    
    if((J1+M2-J3)>t1)
        
        t1=J1+M2-J3;
    
    t2=J1+J2-J3;
    
    if((J1-M1)<t2)
        
        t2=J1-M1;
    
    if ((J2+M2)<t2)
        
        t2=J2+M2;
    
    for (t=t1; t<=t2; t++)
        
        e=e+pow(-1,t)/(1.0*Factorial(t)*Factorial(J1+J2-J3-t)*Factorial(J1-M1-t)*Factorial(J2+M2-t)*Factorial(J3-J2+M1+t)*Factorial(J3-J1-M2+t));
    
    F=a*c*d*e;
    
    return F;
    
}


extern double SJ(float J1,float J2,float J3,float J4,float J5,float J6){
    
    /***********************************************************************************
     Purpose:
     Computes the 6j symbol.
     Modified:
     25 Apr. 2018, Hao Li
     Input parameters:
     J1, J2, J3, J4, J5, J6, The total angular momentum of the electronic cloud.
     Output parameters:
     C, The 6j symbol.
     References:
     LL04 Chapter 2, Page 42, Equation 2.35.
     ***********************************************************************************/
    
    if ((J1+J2)<J3||fabs(J1-J2)>J3)
        
        return 0;
    
    if ((J1+J5)<J6||fabs(J1-J5)>J6)
        
        return 0;
    
    if ((J4+J2)<J6||fabs(J4-J2)>J6)
        
        return 0;
    
    if ((J4+J5)<J3||fabs(J4-J5)>J3)
        
        return 0;
    
    if(fmodf((J1+J2+J3)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmodf((J1+J5+J6)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmodf((J4+J2+J6)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmodf((J4+J5+J3)*2.0,2.0)!=0)
        
        return 0;
    
    double a,b,C;
    
    float t1,t2,t;
    
    a=sqrt(1.0*Factorial(J1+J2-J3)*Factorial(J1-J2+J3)*Factorial(-J1+J2+J3)/Factorial(J1+J2+J3+1.0))*sqrt(1.0*Factorial(J1+J5-J6)*Factorial(J1-J5+J6)*Factorial(-J1+J5+J6)/Factorial(J1+J5+J6+1.0))*sqrt(1.0*Factorial(J4+J2-J6)*Factorial(J4-J2+J6)*Factorial(-J4+J2+J6)/Factorial(J4+J2+J6+1.0))*sqrt(1.0*Factorial(J4+J5-J3)*Factorial(J4-J5+J3)*Factorial(-J4+J5+J3)/Factorial(J4+J5+J3+1.0));
    
    b=0;
    
    t1=J1+J2+J3;
    
    if((J1+J5+J6)>t1)
        
        t1=J1+J5+J6;
    
    if((J4+J2+J6)>t1)
        
        t1=J4+J2+J6;
    
    if((J4+J5+J3)>t1)
        
        t1=J4+J5+J3;
    
    t=t1;
    
    t2=J1+J2+J4+J5;
    
    if((J2+J3+J5+J6)<t2)
        
        t2=J2+J3+J5+J6;
    
    if((J1+J3+J4+J6)<t2)
        
        t2=J1+J3+J4+J6;
    
    for (t=t1; t<=t2; t++)
        
        b=b+pow(-1, t)*Factorial(t+1)/Factorial(t-J1-J2-J3)/Factorial(t-J1-J5-J6)/Factorial(t-J4-J2-J6)/Factorial(t-J4-J5-J3)/(Factorial(J1+J2+J4+J5-t)*Factorial(J2+J3+J5+J6-t)*Factorial(J1+J3+J4+J6-t));
    
    C=a*b;
    
    return C;
    
}


extern double NJ(float J1,float J2,float J3,float J4,float J5,float J6,float J7,float J8,float J9){
    
    /***********************************************************************************
     Purpose:
     Computes the 9j symbol.
     Modified:
     25 Apr. 2018, Hao Li
     Input parameters:
     J1, J2, J3, J4, J5, J6, J7, J8, J9, The total angular momentum of the electronic cloud.
     Output parameters:
     C, The 9j symbol.
     References:
     LL04 Chapter 2, Page 47, Equation 2.48.
     ***********************************************************************************/
    
    if ((J1+J2)<J3||fabs(J1-J2)>J3)
        
        return 0;
    
    if ((J4+J5)<J6||fabs(J4-J5)>J6)
        
        return 0;
    
    if ((J7+J8)<J9||fabs(J7-J8)>J9)
        
        return 0;
    
    if ((J1+J4)<J7||fabs(J1-J4)>J7)
        
        return 0;
    
    if ((J2+J5)<J8||fabs(J2-J5)>J8)
        
        return 0;
    
    if ((J3+J6)<J9||fabs(J3-J6)>J9)
        
        return 0;
    
    if(fmod((J1+J2+J3)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmod((J4+J5+J6)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmod((J7+J8+J9)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmod((J1+J4+J7)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmod((J2+J5+J8)*2.0,2.0)!=0)
        
        return 0;
    
    if(fmod((J3+J6+J9)*2.0,2.0)!=0)
        
        return 0;
    
    float t1,t2,t;
    
    double a=0.0;
    
    t2=J1+J9;
    
    if((t2>J2+J6))
        
        t2=J2+J6;
    
    if(t2>(J4+J8))
        
        t2=J4+J8;
    
    t1=fabs(J1-J9);
    
    if(t1<fabs(J2-J6))
        
        t1=fabs(J2-J6);
    
    if(t1<fabs(J8-J4))
        
        t1=fabs(J8-J4);
    
    t=t1;
    
    for (t=t1; t<=t2; t++)
        
        a=a+pow(-1, 2*t)*(2*t+1)*SJ(J1, J9, t, J8, J4, J7)*SJ(J2, J6, t, J4, J8, J5)*SJ(J1, J9, t, J6, J2, J3);
    
    return a;
    
}
