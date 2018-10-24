#ifndef ANGULAR_MOMENTUM_H

#define ANGULAR_MOMENTUM_H


#include <math.h>


typedef struct {
    
    float Lambda;
    
    float Ju;
    
    float Lu;
    
    float Su;
    
    float Gu;
    
    float Jl;
    
    float Ll;
    
    float Sl;
    
    float Gl;
    
    float Geffect;
    
}TRANSITION_INFORMATION;


extern int Factorial(int N);

extern double Geffect(float Gu, float Gl, float Ju, float Jl);

extern double Gfactor(float J,float L,float S);

extern double TJ(float J1,float J2,float J3,float M1,float M2,float M3);

extern double SJ(float J1,float J2,float J3,float J4,float J5,float J6);

extern double NJ(float J1,float J2,float J3,float J4,float J5,float J6,float J7,float J8,float J9);


#endif /* ANGULAR_MOMENTUM_H */
