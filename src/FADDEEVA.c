
#include "faddeeva.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        Real-valued Faddeeva and Dawson evaluations with reusable
        coefficients and configurable numerical precision.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

// square root of Pi
static const double L_SqrtPi = 1.77245385090551588191;

/*--------------------------------------------------------------------------------*/

static inline double SINC(double x){

    /*######################################################################
      Purpose:
        Evaluate sin(x)/x without losing the removable limit at zero.
      Input parameters:
        x, real argument.
      Return:
        sin(x)/x, with its limiting value used near zero.
     ######################################################################*/

    if(fabs(x) < 1e-4){
      const double x2 = x*x;
      return 1.0-x2/6.0+x2*x2/120.0;
    }

    return sin(x)/x;
}

/*--------------------------------------------------------------------------------*/

static inline double DAWSON(double x){

    /*######################################################################
      Purpose:
        Evaluate Dawson's integral for a real argument.
      Input parameters:
        x, real argument.
      Return:
        Dawson's integral at x.
      Note:
        Numerical recipes in C, 2nd Edition.
     ######################################################################*/

    const double A1 = 2.0/3.0, A2 = 0.4, A3 = 2.0/7.0, h = 0.4;
    const int NMAX = 6;
    double c[NMAX+1];
    double xx = fabs(x), ans;
    
    if(xx<0.2){  
      double x2 = x*x; 
      ans = x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));    
    }else {
      int i;
      double tmp;
      for(i=1;i<=NMAX;i++){ 
        tmp = (2.0*i-1.0)*h;    
        c[i] = exp(-tmp*tmp);    
      }
        
      int n0 = 2*(int)round(xx/(2*h));
      
      double xp = xx-n0*h;
      double e1 = exp(2.0*xp*h);
      double e2 = e1*e1;  
      double d1 = n0+1;
      double d2 = d1-2.0;  
      double sum = 0.0;
        
      for(i=1; i<=NMAX; i++,d1+=2.0,d2-=2.0,e1*=e2){
        sum += c[i]*(e1/d1+1.0/(d2*e1));
      }
      tmp = exp(-xp*xp);
      ans = 0.5641895835*copysign(tmp, x)*sum;
    }
    
    return ans;
}

/*--------------------------------------------------------------------------------*/

static inline double ERFCX(double a){
    
    /*######################################################################
      Purpose:
        Evaluate the scaled complementary error function erfcx(a).
      Input parameters:
        a, real argument.
      Return:
        exp(a*a) multiplied by erfc(a).
      Note:
        Aghloul 2007 MNRAS.
     ######################################################################*/

    if(a>26.6){
      double asqr = a*a;
      return ((((((162.421875/asqr-29.53125)/asqr+6.5625)/asqr-1.875) 
          /asqr+0.75)/asqr-0.5)/asqr+1)/L_SqrtPi/a;

    }else{
      return exp(a*a)*erfc(a);
    }
}

/*--------------------------------------------------------------------------------*/

static inline int Hui_p6(double nu, double y, double *h, double *l){

    /*######################################################################
      Purpose:
        Voigt function（Hui 1978）.
      Input parameters:
        nu, reduced wavelength or frequency shift.
        y, damping parameter.
      Output parameters:
        h, Voigt function.
        l, associated dispersion profile.
      Return:
        0 on success.
      Note:
        Hui 1978 Journal of Quantitative Spectroscopy and Radiative Transfer
     ######################################################################*/

    static const double A[7] = {122.607931777104326,                   
        214.382388694706425, 181.928533092181549, 93.155580458138441,  
        30.180142196210589, 5.912626209773153, 0.564189583562615};
    static const double B[7] = {122.60793177387535,                    
        352.730625110963558, 457.334478783897737, 348.703917719495792, 
        170.354001821091472, 53.992906912940207, 10.479857114260399};

    double z_real = y, z_imag = -nu;
    double sum1_re = A[6], sum1_im = 0.0;
    double sum2_re = B[6]+z_real, sum2_im = z_imag;
    double tmp_re, tmp_im;

    for(int i=5; i>=0; i--){
      // sum1 = sum1*Z+A[i]
      tmp_re = sum1_re*z_real-sum1_im*z_imag+A[i];
      tmp_im = sum1_re*z_imag+sum1_im*z_real;
      sum1_re = tmp_re; 
      sum1_im = tmp_im;

      // sum2 = sum2*Z+B[i]
      tmp_re = sum2_re*z_real-sum2_im*z_imag+B[i];
      tmp_im = sum2_re*z_imag+sum2_im*z_real;
      sum2_re = tmp_re; 
      sum2_im = tmp_im;
    }

    // sum1/sum2;
    double denom = sum2_re*sum2_re+sum2_im*sum2_im;
    *h = (sum1_re*sum2_re+sum1_im*sum2_im)/denom;
    *l = (sum1_im*sum2_re-sum1_re*sum2_im)/denom;

    return 0;
}

/*--------------------------------------------------------------------------------*/

static inline int Hum_W4(double nu, double y, double *h, double *l){

    /*######################################################################
      Purpose:
        Voigt function (HUMLÍČEK 1982).
      Input parameters:
        nu, wavelength.
        y, damping parameter.
      Output parameters:
        h, Voigt function.
        l, associated dispersion profile.
      Note:
        HUMLÍČEK 1982 j. Quant. Spectrosc. Radiat. Transfer.
     ######################################################################*/

    const double abs_nu = fabs(nu);
    const double s = abs_nu +y;
    const double t_real = y, t_imag = -nu;

    if(s >= 15){
      const double u_real = t_real*t_real-t_imag*t_imag;
      const double u_imag = 2.0*t_real*t_imag;
      double d_real = 0.5+u_real, d_imag = u_imag;
      double denom = d_real*d_real+ d_imag*d_imag;

      *h = 0.56418958355*(t_real*d_real+t_imag*d_imag)/denom;
      *l = 0.56418958355*(t_imag*d_real-t_real*d_imag)/denom;

    }else if(s >= 5.5){
      const double u_real = t_real*t_real - t_imag*t_imag;
      const double u_imag = 2.0*t_real*t_imag;
      double n_real = 1.4104739589 + 0.5641896*u_real;
      double n_imag = 0.5641896*u_imag;
      double d_real = 0.75+3.0*u_real+(u_real*u_real-u_imag*u_imag);
      double d_imag = 3.0*u_imag+2.0*u_real*u_imag;
      double denom = d_real*d_real+d_imag*d_imag;
      double tmp_re = (n_real*d_real+n_imag*d_imag)/denom;
      double tmp_im = (n_imag*d_real-n_real*d_imag)/denom;

      *h = t_real*tmp_re-t_imag*tmp_im;
      *l = t_real*tmp_im+t_imag*tmp_re;

    }else if(y >= 1.95*abs_nu-0.176){
      const double c[5] = {0.5642236, 3.778987, 11.96482, 20.20933,     
          16.4955};  
      const double d[6] = {1.0, 6.699398, 21.69274, 39.27121, 38.82363, 
          16.4955};

      double n_real = c[0], n_imag = 0;
      double d_real = d[0], d_imag = 0;
      double tmp_re, tmp_im;

      for(int i=1; i<5; i++){
        tmp_re = n_real*t_real-n_imag*t_imag+c[i];
        tmp_im = n_real*t_imag+n_imag*t_real;
        n_real = tmp_re;
        n_imag = tmp_im;
      }

      for(int i=1; i<6; i++){
        tmp_re = d_real*t_real-d_imag*t_imag + d[i];
        tmp_im = d_real*t_imag+d_imag*t_real;
        d_real = tmp_re;
        d_imag = tmp_im;
      }
      double denom = d_real*d_real+d_imag*d_imag;

      *h = (n_real*d_real+n_imag*d_imag)/denom;
      *l = (n_imag*d_real-n_real*d_imag)/denom;

    }else{
      const double c[7] = {0.5641900381, 1.320521697, 35.7668278,    
          219.0312964, 1540.786893, 3321.990492, 36183.30536};
      const double d[8] = {1, 1.841438936, 61.57036588, 364.2190727, 
          2186.181081, 9022.227659, 24322.84021, 32066.59372};

      const double u_real = t_real*t_real-t_imag*t_imag;
      const double u_imag = 2.0*t_real*t_imag;
      double exp_u = exp(u_real);
      double exp_re = exp_u*cos(u_imag);
      double exp_im = exp_u*sin(u_imag);
      double n_real = c[0], n_imag = 0;
      double d_real = d[0], d_imag = 0;
      double tmp_re, tmp_im;

      for(int i=1; i<7; i++){
        tmp_re = -(n_real*u_real-n_imag*u_imag)+c[i];
        tmp_im = -(n_real*u_imag+n_imag*u_real);
        n_real = tmp_re;
        n_imag = tmp_im;
      }

      tmp_re = n_real*t_real-n_imag*t_imag;
      tmp_im = n_real*t_imag+n_imag*t_real;
      n_real = tmp_re;
      n_imag = tmp_im;

      for(int i=1; i<8; i++){
        tmp_re = -(d_real*u_real-d_imag*u_imag)+d[i];
        tmp_im = -(d_real*u_imag+d_imag*u_real);
        d_real = tmp_re;
        d_imag = tmp_im;
      }
      double denom = d_real*d_real+d_imag*d_imag;

      *h = exp_re-(n_real*d_real+n_imag*d_imag)/denom;
      *l = exp_im-(n_imag*d_real-n_real*d_imag)/denom;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Faddeeva916(double nu, double y, double *h, double *l, 
    STRUCT_FADDEEVA *faddeeva){
    
    /*######################################################################
      Purpose:
        Faddeyeva function (algorithm 916).
      Input parameters:
        nu, reduced wavelength or frequency shift.
        y, damping parameter.
        faddeeva, a structure with precomputed coefficients .
      Output parameters:
        *h, Voigt function.
        *l, associated dispersion profile.
      Note:
        Zaghloul 2011.
     ######################################################################*/
    
    if(nu == 0){
      *h = ERFCX(y);
      *l = 0;
      return 0;
    }else if(y == 0){
      *h = exp(-nu*nu);
      *l = 2.0/L_SqrtPi*DAWSON(nu);
      return 0;
    }
    
    const double x = fabs(nu), ya = fabs(y);

    if(x+ya>1e7){
      if(x<ya){
        const double xs = x/ya;
        const double scale = (1.0/ya)/L_SqrtPi/(xs*xs+1.0);
        *h = scale;
        *l = copysign(xs*scale, nu);
      }else{
        const double ys = ya/x;
        const double scale = (1.0/x)/L_SqrtPi/(ys*ys+1.0);
        *h = ys*scale;
        *l = copysign(scale, nu);
      }
      return 0;
    }

    if(x+ya<1e-4){
      *h = 1.0-2.0*y/L_SqrtPi-nu*nu+y*y;
      *l = 2.0*nu/L_SqrtPi-2.0*nu*y;
      return 0;
    }

    const double XX = x*x, YY = y*y, XY = x*y;
    const double XY2 = 2*XY;
    const double SINXY = sin(XY), SIN2XY = sin(XY2), COS2XY = cos(XY2);
    const double EXPXX = exp(-XX), EXP2AX = exp(2*faddeeva->a*x);
    const double const1 = EXPXX*ERFCX(y);

    const double remaining_log = faddeeva->log_min_value-XX;
    const int max_terms_1 = remaining_log>0.0
        ? (int)(sqrt(remaining_log)/faddeeva->a+0.5) : 0;
    const int max_terms_2 = (int)((faddeeva->sqrt_log_min_value-x)
        /faddeeva->a+0.5);
        
    int max_terms = max_terms_1>max_terms_2?max_terms_1:max_terms_2;
    if(max_terms<0) max_terms = 0;
    if(max_terms>=faddeeva->max_terms) max_terms = faddeeva->max_terms-1;

    int backward_terms = (int)(x/faddeeva->a+0.5);
    int forward_terms = backward_terms+1;

    const double forward_offset = faddeeva->a*forward_terms-x;
    const double backward_offset = x-faddeeva->a*backward_terms;

    int max_forward_terms = (int)((faddeeva->sqrt_log_min_value 
        -forward_offset)/faddeeva->a);
    int max_backward_terms = (int)((faddeeva->sqrt_log_min_value 
        -backward_offset)/faddeeva->a);

    if(max_forward_terms>=faddeeva->max_terms){ 
      max_forward_terms = faddeeva->max_terms-1;
    }
    if(max_backward_terms>=faddeeva->max_terms){ 
      max_backward_terms = faddeeva->max_terms-1;
    }

    const int istart = (backward_terms-max_backward_terms)>1?(
        backward_terms-max_backward_terms):1;
    
    double Re = const1*COS2XY
        +2*faddeeva->a*x*EXPXX*SINXY*SINC(XY)/M_PI;
    double Im = -const1*SIN2XY
        +2*faddeeva->a*x*EXPXX*SINC(XY2)/M_PI;
    double EXP1, EXP2, ep_tmp, delta; 
    double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;

    EXP1 = exp(-forward_offset*forward_offset);
    EXP2 = exp(-2*faddeeva->a*forward_offset);
    ep_tmp = 1.;
    for(int n = forward_terms; n <= forward_terms+max_forward_terms; n++){
      delta = EXP1*faddeeva->exp_a2n2[n-forward_terms]/(
          faddeeva->a_squared*n*n+YY);
      sum3 += delta*ep_tmp;
      sum5 += delta*ep_tmp*faddeeva->a*n;
      ep_tmp *= EXP2;            
    }

    EXP1 = exp(-backward_offset*backward_offset);
    EXP2 = exp(-2*faddeeva->a*backward_offset);
    ep_tmp = 1.;
    for(int n = backward_terms; n >= istart; n--){
      delta = EXP1*faddeeva->exp_a2n2[backward_terms-n]/(faddeeva->a_squared*n*n+YY);
      sum3 += delta*ep_tmp;
      sum5 += delta*ep_tmp*faddeeva->a*n;
      ep_tmp *= EXP2;
    }

    Re += faddeeva->a*y*sum3/M_PI;
    Im += faddeeva->a*sum5/M_PI;

    if(x < faddeeva->sqrt_log_min_value || x < 10){
      ep_tmp = 1.;
      for(int n = 1; n <= max_terms; n++){
        delta = EXPXX*faddeeva->exp_a2n2[n]/(faddeeva->a_squared*n*n+YY);
        sum1 += delta;
        ep_tmp /= EXP2AX;
        sum2 += delta*ep_tmp;
        sum4 += delta*ep_tmp*faddeeva->a*n;
      }

      Re += faddeeva->a*y*(-2*COS2XY*sum1+sum2)/M_PI;
      Im += faddeeva->a*(2*y*SIN2XY*sum1-sum4)/M_PI;
    }

    *h = Re;
    *l = nu>0 ? Im : -Im;

    return 0;
}

/*--------------------------------------------------------------------------------*/

int Faddeeva_init(STRUCT_FADDEEVA *faddeeva){

    /*######################################################################
      Purpose:
        Precompute some coefficients for Faddeyeva function (algorithm 916).
      Input parameters:
        faddeeva, a structure with precomputed coefficients.
      Output parameters:
        faddeeva, a structure with precomputed coefficients.
     ######################################################################*/

    faddeeva->min_value = pow(10.,-faddeeva->precision_digits);
    faddeeva->a = M_PI/sqrt(-log(0.5*faddeeva->min_value));
    faddeeva->a_squared = faddeeva->a*faddeeva->a;
    faddeeva->log_min_value = -log(faddeeva->min_value);
    faddeeva->sqrt_log_min_value = sqrt(faddeeva->log_min_value);
    faddeeva->max_terms = 50;
    faddeeva->exp_a2n2 = (double *)malloc(faddeeva->max_terms*sizeof(double));
    if(!faddeeva->exp_a2n2) return -1;
    faddeeva->precision_idx = (faddeeva->precision_digits>4? 
        faddeeva->precision_digits-4:0);

    faddeeva->exp_a2n2[0] = 1.;
    for(int i=1; i<faddeeva->max_terms; i++){
      faddeeva->exp_a2n2[i] = exp(-faddeeva->a_squared*i*i);
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

int Faddeeva(double nu, double y, double *h, double *l, 
    STRUCT_FADDEEVA *faddeeva){

    /*######################################################################
      Purpose:
        Faddeyeva function with required accuracy (precision_digits). 
      Input parameters:
        nu, reduced wavelength or frequency shift.
        y, damping parameter.
        faddeeva, a structure storing precomputed coefficients and 
          precision_digits.
      Output parameters:
        h, Voigt function.
        l, associated dispersion profile.
      Note:
        Zaghloul 2018 ACM Trans. Math. Soft.
     ######################################################################*/

    if(nu == 0){
      *h = ERFCX(y);
      *l = 0;
      return 0;
    }else if(y == 0){
      *h = exp(-nu*nu);
      *l = 2.0/L_SqrtPi*DAWSON(nu);
      return 0;
    }

    static const double THRESHOLD[5][5] = {         
        {1.6e4, 160.0, 107.0, 28.5, 3.5 },          
        {1.5e5, 510.0, 110.0, 39., -0.1 },          
        {1.451e6, 1.6e3, 180.0, 110.0, -0.1 },      
        {1.5e7, 5.01e3, 380.0, 115.0, 114.0 },      
        {1.3e8, 1.6e4, 810.0, 195.0, 116.0 }        
    };

    const double z_real = nu; const double z_imag = y;
    const double z2_real = z_real*z_real-z_imag*z_imag; 
    const double z2_imag = 2.0*z_real*z_imag;
    const double tmp = z_real*z_real+z_imag*z_imag;
    double TMP_y = y*y;

    if(tmp >= THRESHOLD[faddeeva->precision_idx][0]){
      //I/Z/L_SqrtPi;
      double denom = L_SqrtPi*(z_real*z_real+z_imag*z_imag);
      *h = z_imag/denom;
      *l = z_real/denom;

    }else if(tmp >= THRESHOLD[faddeeva->precision_idx][1]){
      //I*Z/L_SqrtPi/(Z2-0.5);
      double z2n_real = z2_real-0.5;
      double denom = L_SqrtPi*(z2n_real*z2n_real+z2_imag*z2_imag);
      *h = (-z_imag*z2n_real+z_real*z2_imag)/denom;
      *l = (z_real*z2n_real+z_imag*z2_imag)/denom;

    }else if(tmp >= THRESHOLD[faddeeva->precision_idx][2]){
      //(Z2-1)/(Z2-1.5)*I/Z/L_SqrtPi;
      double z2n_real = z2_real-1.5;
      double Z2Z_re = z2n_real*z_real-z2_imag*z_imag;
      double Z2Z_im = z2n_real*z_imag+z2_imag*z_real; 
      double denom = L_SqrtPi*(Z2Z_re*Z2Z_re+Z2Z_im*Z2Z_im);
      z2n_real = z2_real-1.0;
      *h = (-z2_imag*Z2Z_re+z2n_real*Z2Z_im)/denom;
      *l = (z2n_real*Z2Z_re+z2_imag*Z2Z_im)/denom;

    }else if(tmp >= THRESHOLD[faddeeva->precision_idx][3]){

      if(faddeeva->precision_digits<=4 && TMP_y<6e-14){
        Hum_W4(nu, y, h, l);
        return 0;
      }else if(faddeeva->precision_digits==5 && tmp<39 && TMP_y<1e-9){
        Hum_W4(nu, y, h, l);
        return 0;
      }

      //(Z2-2.5)/(Z2*(Z2-3)+0.75)/L_SqrtPi*Z*I;
      double z2n_real = z2_real-2.5;
      // A = (Z2-2.5)*Z*I;
      double A_re = -(z2n_real*z_imag+z2_imag*z_real);
      double A_im = z2n_real*z_real-z2_imag*z_imag;

      z2n_real = z2_real-3.;
      // B = Z2*(Z2-3)+0.75
      double B_re = z2n_real*z2_real-z2_imag*z2_imag+0.75;
      double B_im = (z2n_real+z2_real)*z2_imag;

      double denom = L_SqrtPi*(B_re*B_re+B_im*B_im);

      *h = (A_re*B_re+A_im*B_im)/denom;
      *l = (-A_re*B_im+A_im*B_re)/denom;

    }else if(tmp >= THRESHOLD[faddeeva->precision_idx][4]){
      if(faddeeva->precision_digits<=4){
        if(TMP_y<0.026){
          Hum_W4(nu, y, h, l);
        }else{
          Hui_p6(nu, y, h, l);
        }
        return 0;
      }else if(faddeeva->precision_digits==5){
        if(TMP_y>=0.27){
          Hui_p6(nu, y, h, l);
        }else if(TMP_y>=1e-9){
          Faddeeva916(nu, y, h, l, faddeeva);
        }else{
          Hum_W4(nu, y, h, l);
        }
        return 0;
      }else if(faddeeva->precision_digits==6){
        if(TMP_y>=1){
          Hui_p6(nu, y, h, l);
        }else{
          Faddeeva916(nu, y, h, l, faddeeva);
        }
        return 0;
      }

      //(Z2*(Z2-4.5)+2.0)/(Z2*(Z2-5)+3.75)/L_SqrtPi/Z*I;
      double z2n_real = z2_real-4.5;
      // A = (Z2*(Z2-4.5)+2.0)*I
      double A_re = -(z2n_real+z2_real)*z2_imag;
      double A_im = z2n_real*z2_real-z2_imag*z2_imag+2.;

      z2n_real = z2_real-5.;
      // B = Z2*(Z2-5)+3.75
      double B_re = z2n_real*z2_real-z2_imag*z2_imag+3.75;
      double B_im = (z2n_real+z2_real)*z2_imag;

      // C = B*Z
      double C_re = B_re*z_real-B_im*z_imag;
      double C_im = B_re*z_imag+B_im*z_real;

      double denom = L_SqrtPi*(C_re*C_re+C_im*C_im);

      *h = (A_re*C_re+A_im*C_im)/denom;
      *l = (-A_re*C_im+A_im*C_re)/denom;

    }else{
      if(faddeeva->precision_digits<=4){
        Hui_p6(nu, y, h, l);
        return 0;
      }else{
        Faddeeva916(nu, y, h, l, faddeeva);
      }
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/
