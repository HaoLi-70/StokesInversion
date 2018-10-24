#include "FADDEEVA.h"


extern double DAWSON(double x){
    
    /***********************************************************************************
     Purpose:
     Dawson’s Integral
     Modified:
     7 May 2018, Hao Li
     Input parameters:
     x, the real number.
     Output parameters:
     ans,the integral.
     Reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    double A1=2.0/3.0;
    
    double A2=0.4;
    
    double A3=2.0/7.0;
    
    double H=0.4;
    
    int NMAX=6;
    
    int i,n0;
    
    double d1,d2,e1,e2,sum,x2,xp,xx,ans;
    
    double c[NMAX+1];
    
    double sqrarg, TMP;
    
    xx=fabs(x);
    
    if (xx < 0.2) {
        
        x2=x*x;
        
        ans=x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));
        
    } else {
        
        for (i=1;i<=NMAX;i++){
            
            TMP=(2.0*i-1.0)*H;
            
            sqrarg=(sqrarg=(TMP)) == 0.0 ? 0.0 : sqrarg*sqrarg;
            
            c[i]=exp(-sqrarg);
            
        }
        
        n0=2*(int)(0.5*xx/H+0.5);
        
        xp=xx-n0*H;
        
        e1=exp(2.0*xp*H);
        
        e2=e1*e1;
        
        d1=n0+1;
        
        d2=d1-2.0;
        
        sum=0.0;
        
        for (i=1;i<=NMAX;i++,d1+=2.0,d2-=2.0,e1*=e2)
            
            sum += c[i]*(e1/d1+1.0/(d2*e1));
        
        TMP=((x) >= 0.0 ? fabs(exp(-xp*xp)) : -fabs(exp(-xp*xp)));
        
        ans=0.5641895835*TMP*sum;
        
    }
    
    return ans;
}


extern double Fadd_erfcx(double a){
    
    /***********************************************************************************
     Purpose:
     Scaled complementary error function
     Modified:
     7 May 2018, Hao Li
     Input parameters:
     a, the real number.
     Output parameters:
     erfcx, the result.
     Reference:
     Mofreh R. Zaghloul 2007 MNRAS.
     ***********************************************************************************/
    
    double erfcx=0;
    
    double C_Pi=3.141592653589793;
    
    if (a>26.6)
        
        erfcx=((((((162.421875/a/a-29.53125)/a/a+6.5625)/a/a-1.875)/a/a+0.75)/a/a-0.5)/a/a+1)/sqrt(C_Pi)/a;
    
    else
        
        erfcx=exp(a*a)*erfc(a);
    
    return erfcx;
    
}


extern void Faddeeva_Function_V1(double Nu, double y, double *H, double *L){
    
    /***********************************************************************************
     Purpose:
     Scaled complementary error function (old version)
     Modified:
     7 May 2018, Hao Li
     Input parameters:
     a, the real number.
     Output parameters:
     erfcx, the result.
     Reference:
     Mofreh R. Zaghloul 2011 ACM Trans. Math. Soft..
     ***********************************************************************************/
    
    double C_Pi=3.141592653589793;
    
    double x=fabs(Nu);
    
    double Re=0,Im=0;
    
    double Rmin=1e-20;
    
    double a=sqrt(-C_Pi*C_Pi/log(Rmin/2));
    
    double Max1=-log(Rmin)-x*x>0? (sqrt(-log(Rmin)-x*x)/a>1? sqrt(-log(Rmin)-x*x)/a :1) :1;
    
    double Max2=(sqrt(-log(Rmin))-x)/a>1? (sqrt(-log(Rmin))-x)/a :1;
    
    double Max=Max1>Max2? Max1: Max2;
    
    int n,m,i,j;
    
    double EXP1,EXP2,EXP3,delta,fu1;
    
    double Sigma1=0,Sigma2=0,Sigma3=0,Sigma4=0,Sigma5=0;

    EXP1=exp(-x*x);
    
    EXP2=exp(2*a*x);
    
    fu1=exp(-x*x)*Fadd_erfcx(y);
    
    if (y==0) {
        
        Re=fu1*cos(2*x*y)+2*a*x*sin(x*y)*EXP1/C_Pi;
        
        Im=-fu1*sin(2*x*y)+2*a*x*EXP1/C_Pi;
        
    }else{
        
        Re=fu1*cos(2*x*y)+2*a*x*sin(x*y)*EXP1*sin(x*y)/x/y/C_Pi;
        
        Im=-fu1*sin(2*x*y)+2*a*x*EXP1*sin(2*x*y)/2/x/y/C_Pi;
        
    }
    
    for (n=1; n<=Max; n++) {
        
        EXP3=exp(-a*a*n*n)/(a*a*n*n+y*y);
        
        delta=EXP1*EXP3;
        
        for (m=1; m<=n; m++)
            
            delta/=EXP2;
    
        Sigma1+=EXP1*EXP3;
        
        Sigma2+=delta;
        
        Sigma4+=delta*a*n;
        
    }
    
    j=x/a>1? x/a:1;
    
    i=j+1;
    
    do{
        
        delta=exp(-(a*i-x)*(a*i-x))/(a*a*i*i+y*y);
        
        Sigma3+=delta;
        
        Sigma5+=delta*a*i;
        
        i++;
        
    }while (fabs(delta)>Rmin);
    
    do {
        
        delta=exp(-(a*j-x)*(a*j-x))/(a*a*j*j+y*y);
        
        Sigma3+=delta;
        
        Sigma5+=delta*a*j;
        
        j--;
        
    } while ((fabs(delta)>Rmin&&j>0));
    
    Re+=2*a*(-y*cos(2*x*y)*Sigma1+y/2*Sigma2+y/2*Sigma3)/C_Pi;
    
    Im+=2*a*(y*sin(2*x*y)*Sigma1-0.5*Sigma4+0.5*Sigma5)/C_Pi;
    
    if (Nu==0){
        
        *H=Fadd_erfcx(y);
        
        *L=0;
        
    }else if (Nu>0){
        
        *H=Re;
        
        *L=Im;
        
    }else{
        
        *H=Re;
        
        *L=-Im;
        
    }
    
    return ;
    
}


extern void Faddeeva_Function(double Nu, double y, double *H, double *L){
    
    /***********************************************************************************
     Purpose:
     Scaled complementary error function
     Modified:
     7 May 2018, Hao Li
     Input parameters:
     a, the real number.
     Output parameters:
     erfcx, the result.
     Reference:
     Mofreh R. Zaghloul 2011 ACM Trans. Math. Soft..
     ***********************************************************************************/
    
    double C_Pi=3.141592653;
    
    if (Nu==0){
        
        *H=Fadd_erfcx(y);
        
        *L=0;
        
        return;
        
    }
    
    if (y==0) {
        
        *H=exp(-Nu*Nu);
        
        *L=2.0/sqrt(C_Pi)*DAWSON(Nu);
        
        return;
        
    }
    
    double x=fabs(Nu);
    
    double Re=0,Im=0;
    
    double Tiny=1e-8;
    
    double a=sqrt(-C_Pi*C_Pi/log(Tiny/2));
    
    int n,m,i,j;
    
    double EXP1,EXP2,EXP3,delta,fu1;
    
    double Sigma1=0,Sigma2=0,Sigma3=0,Sigma4=0,Sigma5=0;
    
    float Max1=-log(Tiny)-x*x>0? (sqrt(-log(Tiny)-x*x)/a>1? sqrt(-log(Tiny)-x*x)/a :1) :1;
    
    float Max2=(sqrt(-log(Tiny))-x)/a>1? (sqrt(-log(Tiny))-x)/a :1;
    
    float Max=Max1>Max2? Max1: Max2;
    
    EXP1=exp(-x*x);
    
    EXP2=exp(2*a*x);
    
    fu1=EXP1*Fadd_erfcx(y);
    
    j=(x/a>1? x/a:1);
    
    Re=fu1*cos(2*x*y)+2*a*x*sin(x*y)*EXP1*sin(x*y)/x/y/C_Pi;
    
    Im=-fu1*sin(2*x*y)+2*a*x*EXP1*sin(2*x*y)/2/x/y/C_Pi;
    
    for (n=1; n<=Max; n++) {
        
        EXP3=exp(-a*a*n*n)/(a*a*n*n+y*y);
        
        delta=EXP1*EXP3;
        
        for (m=1; m<=n; m++)
            
            delta/=EXP2;
        
        Sigma1+=EXP1*EXP3;
        
        Sigma2+=delta;
        
        Sigma4+=delta*a*n;
        
    }
    
    j=x/a>1? x/a:1;
    
    i=j+1;
    
    do{
        
        delta=exp(-(a*i-x)*(a*i-x))/(a*a*i*i+y*y);
        
        Sigma3+=delta;
        
        Sigma5+=delta*a*i;
        
        i++;
        
    }while (fabs(delta)>Tiny);
    
    do {
        
        delta=exp(-(a*j-x)*(a*j-x))/(a*a*j*j+y*y);
        
        Sigma3+=delta;
        
        Sigma5+=delta*a*j;
        
        j--;
        
    } while ((fabs(delta)>Tiny&&j>0));
    
    Re+=2*a*(-y*cos(2*x*y)*Sigma1+y/2*Sigma2+y/2*Sigma3)/C_Pi;
    
    Im+=2*a*(y*sin(2*x*y)*Sigma1-0.5*Sigma4+0.5*Sigma5)/C_Pi;
    
    if (Nu>0){
        
        *H=Re;
        
        *L=Im;
        
    }else{
        
        *H=Re;
        
        *L=-Im;
        
    }
    
    return;
    
}
