#include "FADDEEVA.h"

/*--------------------------------------------------------------------------------*/

    /*######################################################################
     
      revision log:
        30 Oct. 2022.
     
     ######################################################################*/

/*--------------------------------------------------------------------------------*/

static double DAWSON(double x);

static double Fadd_erfcx(double a);

static int Hui_p6(double Nu, double y, double *H, double *L);

static int Hum_W4(double Nu, double y, double *H, double *L);

static int Exp_calculation(double *Rexp, double *Rmin, int digits);

static int Faddeeva_4digits(double Nu, double y, double *H, double *L);

static int Faddeeva_5digits(double Nu, double y, double *H, double *L);

static int Faddeeva_6digits(double Nu, double y, double *H, double *L);

static int Faddeeva_7digits(double Nu, double y, double *H, double *L);

static int Faddeeva_8digits(double Nu, double y, double *H, double *L);

/*--------------------------------------------------------------------------------*/

static double DAWSON(double x){
    
    /*######################################################################
      Purpose:
        Dawson’s Integral
      Record of revisions:
        7 May 2018
      Input parameters:
        x, the real number.
      Output parameters:
        ans,the integral.
      Reference:
        Numerical recipes in C 2ed.
     ######################################################################*/
    
    double A1 = 2.0/3.0, A2 = 0.4, A3 = 2.0/7.0, H = 0.4;
    int NMAX = 6;
    int i, n0;
    double d1, d2, e1, e2, sum, x2, xp, xx, ans, sqrarg, TMP;
    double c[NMAX+1];
    
    xx = fabs(x);
    
    if(xx<0.2){  
      x2 = x*x; 
      ans = x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));    
    }else {
      for(i=1;i<=NMAX;i++){ 
        TMP=(2.0*i-1.0)*H;    
        sqrarg=(sqrarg=(TMP)) == 0.0 ? 0.0 : sqrarg*sqrarg;    
        c[i]=exp(-sqrarg);    
      }
        
      n0 = 2*(int)(0.5*xx/H+0.5);
      xp = xx-n0*H;
      e1 = exp(2.0*xp*H);
      e2 = e1*e1;  
      d1 = n0+1;
      d2 = d1-2.0;  
      sum = 0.0;
        
      for(i=1; i<=NMAX; i++,d1+=2.0,d2-=2.0,e1*=e2) \
          sum += c[i]*(e1/d1+1.0/(d2*e1));
      TMP = (x) >= 0.0 ? fabs(exp(-xp*xp)) : -fabs(exp(-xp*xp));
      ans = 0.5641895835*TMP*sum;
    }
    
    return ans;
}

/*--------------------------------------------------------------------------------*/

static double Fadd_erfcx(double a){
    
    /*######################################################################
      Purpose:
        Scaled complementary error function
      Record of revisions::
        7 May 2018
      Input parameters:
        a, the real number.
      Output parameters:
        erfcx, the result.
      Reference:
        aghloul 2007 MNRAS.
     ######################################################################*/
    
    double erfcx = 0;
    
    if(a>26.6){
      erfcx = ((((((162.421875/a/a-29.53125)/a/a+6.5625)/a/a-1.875) \
          /a/a+0.75)/a/a-0.5)/a/a+1)/sqrt(Par_Pi)/a;
    }else{
      erfcx = exp(a*a)*erfc(a);
    }
    
    return erfcx; 
}

/*--------------------------------------------------------------------------------*/

static int Hui_p6(double Nu, double y, double *H, double *L){
    
    /*######################################################################
      Purpose:
        Voigt function.
      Record of revisions:
        07 May 2023
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Hui 1978 Journal of Quantitative Spectroscopy and Radiative Transfer
     ######################################################################*/

    double A[7] = {122.607931777104326, 214.382388694706425, \
        181.928533092181549, 93.155580458138441, 30.180142196210589, \
        5.912626209773153, 0.564189583562615};
    double B[7] = {122.60793177387535, 352.730625110963558, \
        457.334478783897737, 348.703917719495792, 170.354001821091472, \
        53.992906912940207, 10.479857114260399};
    
    double complex Z = y-Nu*I;
    double complex sum1=A[6], sum2=B[6]+Z;
    double complex tmp;
    int i;

    for(i=5; i>-1; i--){
      sum1 = sum1*Z+A[i];
      sum2 = sum2*Z+B[i];
    }

    tmp = sum1/sum2;
    *H = creal(tmp);
    *L = cimag(tmp);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Hum_W4(double Nu, double y, double *H, double *L){
    
    /*######################################################################
      Purpose:
        Voigt function.
      Record of revisions:
        28 Dec, 2018
      Input parameters:
        Nu, wavelength.
        y, damping parameter.
      Output parameters:
        erfcx, the result.
      Reference:
        HUMLÍČEK 1982 J. Quant. Spectrosc. Radiat. Transfer.
     ######################################################################*/
    
    complex double W4, U, T = y-Nu*I;
    double S = fabs(Nu)+y;
    
    if(S >= 15){
      W4 = T*0.56418958355/(0.5+T*T);
    }else if(S >= 5.5){
      U = T*T;
      W4 = T*(1.4104739589+U*0.5641896)/(0.75+U*(3.+U));
    }else if(y >= 1.95*fabs(Nu)-0.176){
      W4 = (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+T*0.5642236)))) \
        /(16.4955+T*(38.82363+T*(39.27121+T*(21.69274+T*(6.699398+T)))));
    }else {
      U = T*T;
      W4 = cexp(U)-T*(36183.30536-U*(3321.990492-U*(1540.786893-U \
          *(219.0312964-U*(35.7668278-U*(1.320521697-U*0.5641900381)))))) \
          /(32066.59372-U*(24322.84021-U*(9022.227659-U*(2186.181081 \
          -U*(364.2190727-U*(61.57036588-U*(1.841438936-U))))))); \
    }
    *H = creal(W4);
    *L = cimag(W4);
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Exp_calculation(double *Rexp, double *Rmin, int digits){

    /*######################################################################
      Purpose:
        Exponent values used in algorithm 916.
      Record of revisions:
        7 May 2018
      Input parameters:
        digits, the accuracy.
      Output parameters:
        Rexp, exponent values.
        Rmin, cut value.
     ######################################################################*/

    switch(digits){
      case 5:
        *Rmin = 1e-5;
        Rexp[0] = 1;
        Rexp[1] = 4.45489550e-01;
        Rexp[2] = 3.93867443e-02;
        Rexp[3] = 6.91094900e-04;
        Rexp[4] = 2.40658018e-06;
        Rexp[5] = 1.66317529e-09;
        Rexp[6] = 2.28113377e-13;
        Rexp[7] = 6.20924113e-18;
        Rexp[8] = 3.35429565e-23;
        Rexp[9] = 3.59616156e-29;
        Rexp[10] = 7.65159708e-36;
        Rexp[11] = 3.23102287e-43;
        Rexp[12] = 2.70771539e-51;
        Rexp[13] = 4.50340532e-60;
        Rexp[14] = 1.48646299e-69;
        Rexp[15] = 9.73738201e-80;
        break;

      case 6:
        *Rmin = 1e-6;
        Rexp[0] = 1;
        Rexp[1] = 5.06487213e-01;
        Rexp[2] = 6.58072802e-02;
        Rexp[3] = 2.19339257e-03;
        Rexp[4] = 1.87540801e-05;
        Rexp[5] = 4.11350601e-08;
        Rexp[6] = 2.31454418e-11;
        Rexp[7] = 3.34084084e-15;
        Rexp[8] = 1.23703808e-19;
        Rexp[9] = 1.17502558e-24;
        Rexp[10] = 2.86317929e-30;
        Rexp[11] = 1.78972679e-36;
        Rexp[12] = 2.86986788e-43;
        Rexp[13] = 1.18052187e-50;
        Rexp[14] = 1.24572777e-58;
        Rexp[15] = 3.37216816e-67;
        break;
            
      case 7:
        *Rmin = 1e-7;
        Rexp[0] = 1;
        Rexp[1] = 5.55946303e-01;
        Rexp[2] = 9.55281542e-02;
        Rexp[3] = 5.07335929e-03;
        Rexp[4] = 8.32770909e-05;
        Rexp[5] = 4.22494602e-07;
        Rexp[6] = 6.62494806e-10;
        Rexp[7] = 3.21077212e-13;
        Rexp[8] = 4.80952495e-17;
        Rexp[9] = 2.22669418e-21;
        Rexp[10] = 3.18628530e-26;
        Rexp[11] = 1.40920585e-31;
        Rexp[12] = 1.92632655e-37;
        Rexp[13] = 8.13862568e-44;
        Rexp[14] = 1.06276669e-50;
        Rexp[15] = 4.28934011e-58;
        break;
            
      case 8:
        *Rmin = 1e-8;
        Rexp[0] = 1;
        Rexp[1] = 5.96688914e-01;
        Rexp[2] = 1.26762815e-01;
        Rexp[3] = 9.58808161e-03;
        Rexp[4] = 2.58206699e-04;
        Rexp[5] = 2.47570691e-06;
        Rexp[6] = 8.45136558e-09;
        Rexp[7] = 1.02718930e-11;
        Rexp[8] = 4.44498217e-15;
        Rexp[9] = 6.84834292e-19;
        Rexp[10] = 3.75661694e-23;
        Rexp[11] = 7.33675917e-28;
        Rexp[12] = 5.10161350e-33;
        Rexp[13] = 1.26300998e-38;
        Rexp[14] = 1.11327370e-44;
        Rexp[15] = 3.49375987e-51;
        break;
            
      default:
        fprintf(stderr, "Digit error in Rexp");
        exit(1);
        break;
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_4digits(double Nu, double y, double *H, double *L){
    
    /*######################################################################
      Purpose:
        Faddeyeva function with 4-digit accuracy.
      Record of revisions:
        28 Dec, 2018
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H=exp(-Nu*Nu);
      *L=2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y, TMP_y = y*y;
    
    if(TMP>=1.6e4){
      W = I/Z/Par_SqrtPi;
    }else if(TMP>=160.){
      W = I*Z/Par_SqrtPi/(TMP1-0.5);
    }else if(TMP>=107){
      W = (TMP1-1)/(TMP1-1.5)*I/Z/Par_SqrtPi;
    }else if(TMP>=28.5){
      if(TMP_y>=6e-14){
        W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
      }else {
        Hum_W4(Nu, y, H, L);
        return 0;
      }
    }else if((TMP>=3.5)&&(TMP_y<0.026)){
      Hum_W4(Nu, y, H, L);
      return 0;
    }else {
      Hui_p6(Nu, y, H, L);
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_5digits(double Nu, double y, double *H, double *L){
    
    /*######################################################################
      Purpose:
        Faddeyeva function with 5-digit accuracy.
      Record of revisions:
        28 Dec, 2018
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y, TMP_y = y*y;
    
    if(TMP>=1.5e5){
      W = I/Z/Par_SqrtPi;
    }else if(TMP>=510.){
      W = I*Z/Par_SqrtPi/(TMP1-0.5);
    }else if(TMP>=110){
      W = (TMP1-1)/(TMP1-1.5)*I/Z/Par_SqrtPi;
    }else if(TMP>=109){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP>=39){
      if(TMP_y>=1e-9){
        W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
      }else {
        Hum_W4(Nu, y, H, L);
        return 0;
      }
    }else {
      if(TMP_y>=0.27){
        Hui_p6(Nu, y, H, L);
      }else if(TMP_y>=1e-9){
        Faddeeva916(Nu, y, H, L, 5);
      }else{
        Hum_W4(Nu, y, H, L);
      }
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_6digits(double Nu, double y, double *H, double *L){
    
    /*######################################################################
      Purpose:
        Faddeyeva function with 6-digit accuracy.
      Record of revisions:
        28 Dec, 2018
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I, TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y, TMP_y = y*y;
    
    if(TMP>=1.451e6){
      W = I/(Z*Par_SqrtPi);
    }else if(TMP>=1.6e3){
      W = I*Z/(Par_SqrtPi*(TMP1-0.5));
    }else if(TMP>=180){
      W = (TMP1-1)/(TMP1-1.5)*I/(Z*Par_SqrtPi);
    }else if(TMP>=110){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP_y<1){
      Faddeeva916(Nu, y, H, L, 6);
      return 0;
    }else {
      Hui_p6(Nu, y, H, L);
      return 0;
    }
    
    *H=creal(W);
    *L=cimag(W);
    return 0;   
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_7digits(double Nu, double y, double *H, double *L){
    
    /*######################################################################
      Purpose:
        Faddeyeva function with 7-digit accuracy.
      Record of revisions:
        28 Dec, 2018
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }
    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y;
    
    if(TMP>=1.5e7){
      W = I/(Z*Par_SqrtPi);
    }else if(TMP>=5.01e3){
      W = I*Z/(Par_SqrtPi*(TMP1-0.5));
    }else if(TMP>=380){
      W = (TMP1-1)/(TMP1-1.5)*I/(Z*Par_SqrtPi);
    }else if(TMP>=115){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP>=114){
      W = (TMP1*(TMP1-4.5)+2.0)/(TMP1*(TMP1-5)+3.75)/Par_SqrtPi/Z*I;
    }else {
      Faddeeva916(Nu, y, H, L, 7);
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);
    return 0;
    
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva_8digits(double Nu, double y, double *H, double *L){
    
    /*######################################################################
      Purpose:
        Faddeyeva function with 8-digit accuracy.
      Record of revisions:
        28 Dec, 2018
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft..
     ######################################################################*/
    
    if(Nu==0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }

    if(y==0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    complex double W, Z = Nu+y*I;
    complex double TMP1 = Z*Z;
    double TMP = Nu*Nu+y*y;
    
    if(TMP>=1.3e8){
      W = I/(Z*Par_SqrtPi);
    }else if(TMP>=1.6e4){
      W = I*Z/(Par_SqrtPi*(TMP1-0.5));
    }else if(TMP>=810){
      W = (TMP1-1)/(TMP1-1.5)*I/(Z*Par_SqrtPi);
    }else if(TMP>=195){
      W = (TMP1-2.5)/(TMP1*(TMP1-3)+0.75)/Par_SqrtPi*Z*I;
    }else if(TMP>=116){
      W = (TMP1*(TMP1-4.5)+2.0)/(TMP1*(TMP1-5)+3.75)/Par_SqrtPi/Z*I;
    }else {
      Faddeeva916(Nu, y, H, L, 8);
      return 0;
    }
    
    *H = creal(W);
    *L = cimag(W);

    return 0;    
}

/*--------------------------------------------------------------------------------*/

extern int Faddeeva(double Nu, double y, double *H, double *L, int digits){

    /*######################################################################
      Purpose:
        Faddeyeva function.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        digits, the accuracy.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2018 ACM Trans. Math. Soft.
     ######################################################################*/

    switch(digits){
      case 4:
        Faddeeva_4digits(Nu, y, H, L);
        break;
      case 5:
        Faddeeva_5digits(Nu, y, H, L);
        break;
      case 6:
        Faddeeva_6digits(Nu, y, H, L);
        break;
      case 7:
        Faddeeva_7digits(Nu, y, H, L);
        break;
      case 8:
        Faddeeva_8digits(Nu, y, H, L);
        break;
      default:
        Faddeeva916(Nu, y, H, L, digits);
        break;
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

extern int Faddeeva916(double Nu, double y, double *H, double *L, int digits){
    
    /*######################################################################
      Purpose:
        Faddeyeva function.
      Record of revisions:
        7 May 2018
      Input parameters:
        Nu, reduced wavelength or frequency shift.
        y, damping parameter.
        digits, the accuracy.
      Output parameters:
        H, Voigt function.
        L, associated dispersion profile.
      Reference:
        Zaghloul 2011.
     ######################################################################*/
    
    if(Nu == 0){
      *H = Fadd_erfcx(y);
      *L = 0;
      return 0;
    }else if(y == 0){
      *H = exp(-Nu*Nu);
      *L = 2.0/Par_SqrtPi*DAWSON(Nu);
      return 0;
    }
    
    int n;
    double Rexp[16];
    double Rmin = 1;
    Exp_calculation(Rexp, &Rmin, digits);
    
    double tmp1 = -log(Rmin);
    double tmp2 = sqrt(tmp1);
    double a = Par_Pi/sqrt(-log(0.5*Rmin));
    
    double x = fabs(Nu);
    int Max = (int)(sqrt(tmp1-x*x)/a),Max1, Max2;
    
    double Sigma1 = 0, Sigma2 = 0, Sigma3 = 0, Sigma4 = 0, Sigma5 = 0;
    double Re = 0, Im = 0, EXP1, EXP2, EXP3, delta, fu1;
    EXP1 = exp(-x*x);
    EXP2 = exp(2*a*x);
    EXP3 = 1/EXP2;
    fu1 = EXP1*Fadd_erfcx(y);
    Re = fu1*cos(2*x*y)+2*a*x*sin(x*y)*EXP1*sin(x*y)/x/y/Par_Pi;
    Im = -fu1*sin(2*x*y)+2*a*x*EXP1*sin(2*x*y)/2/x/y/Par_Pi;
    
    int Backward, Forward;
    Backward = (int)(x/a);
    Forward = Backward+1;
    
    double dx1 = a*Forward-x;
    double dx2 = x-a*Backward;
    Max1 = (int)((tmp2-dx1)/a);
    Max2 = (tmp2-dx2)/a>1?(int)((tmp2-dx2)/a):1;
    if(Backward<=Max2) Max2 = Backward-1;
    
    double EP2 = 1, EP3 = 1;
    if(x < tmp2){
      for(n = 1; n <= Max; n++){
        delta = EXP1*Rexp[n]/(a*a*n*n+y*y);
        Sigma1 += delta;
        EP2 *= EXP2;
        EP3 *= EXP3;
        Sigma2 += delta*EP3;
        Sigma4 += delta*EP3*a*n;
      }
      EXP1 = exp(-dx1*dx1);
      EXP2 = exp(-2*a*dx1);
      EP2 = 1;
      for(n = Forward; n <= Forward+Max1; n++){
        delta = EXP1*Rexp[n-Forward]/(a*a*n*n+y*y);
        Sigma3 += delta*EP2;
        Sigma5 += delta*EP2*a*n;
        EP2 *= EXP2;
      }
        
      EXP1 = exp(-dx2*dx2);
      EXP2 = exp(-2*a*dx2);
      EP2 = 1;
        
      for(n = Backward; n >= Backward-Max2; n--){
        delta = EXP1*Rexp[Backward-n]/(a*a*n*n+y*y);
        Sigma3 += delta*EP2;
        Sigma5 += delta*EP2*a*n;
        EP2 *= EXP2;
      }
    }else{
      EXP1 = exp(-dx1*dx1);
      EXP2 = exp(-2*a*dx1);
      EP2 = 1;
      for(n = Forward; n <= Forward+Max1; n++){
        delta = EXP1*Rexp[n-Forward]/(a*a*n*n+y*y);
        Sigma3 += delta*EP2;
        Sigma5 += delta*EP2*a*n;
        EP2 *= EXP2;            
      }
        
      EXP1 = exp(-dx2*dx2);
      EXP2 = exp(-2*a*dx2);
      EP2 = 1;
        
      for(n = Backward; n >= Backward-Max2; n--){
        delta = EXP1*Rexp[Backward-n]/(a*a*n*n+y*y);
        Sigma3 += delta*EP2;
        Sigma5 += delta*EP2*a*n;
        EP2 *= EXP2;
      }
    }
    Re += 2*a*(-y*cos(2*x*y)*Sigma1+y/2*Sigma2+y/2*Sigma3)/Par_Pi;
    Im += 2*a*(y*sin(2*x*y)*Sigma1-0.5*Sigma4+0.5*Sigma5)/Par_Pi;
    
    *H = Re;
    if(Nu>0){
      *L = Im;
    }else{
      *L = -Im;
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

