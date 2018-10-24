#include <stdio.h>

#include <math.h>

#include <time.h>

#include "ANGULAR_MOMENTUM.h"

#include "ALLOCATION.h"

#include "MCMC.h"

#include "INIT.h"

#include "PATH.h"

#include "INVERSION.h"


int main(int argc, const char * argv[]) {
    
    ///****************** Begin to calculate the running time **********///
    clock_t timea=clock();
    
    
    ///****************************** Initial the parameter ******************************///
    //The wavelength number of the profile and the number of parameter
    
    int Num_Wav;
    
    Init_Wav(&Num_Wav);
    
    int i, j;
    
    FILE *File1, *File2;
    
    File1=fopen(Path_Profile, "r");
    
    File2=fopen(Path_RESULT, "w");
    
    
    ///****************** Parameters for the Inversion Procedure **********///
    /* If Flag_Correlation=1, save the correlation coefficients.
     If Flag_Histdata=1, save the hist data.  */
    
    int Flag_Correction=1, Flag_Histdata=1;
    
    double Accept_Rates;
    

    ///****************************** Memory allocation ******************************///
    
    double **Profile_Data;
    
    Profile_Data=MATRIX_DOUBLE(0, Num_Wav, 0, 3);
    
    double *Wavelength_Data;
    
    Wavelength_Data=VECTOR_DOUBLE(0, Num_Wav);
    
    double **Result;
    
    Result=MATRIX_DOUBLE(0, 11, 0, 1);
    
    
    ///****************************** Input the profile data ******************************///

    for (i=0; i<Num_Wav; i++) {
        fscanf(File1, "%le ",&Wavelength_Data[i]);
        for (j=0; j<4; j++) {
            fscanf(File1, "%le ",&Profile_Data[i][j]);
        }
    }
    
    fclose(File1);
    
    Inversion(Profile_Data, Wavelength_Data, Num_Wav, Flag_Correction, Flag_Histdata, Result, &Accept_Rates);
    
    for (i=0; i<12; i++)
        
        fprintf(File2,"%e %e \n",Result[i][0],Result[i][1]);

    fclose(File2);
    
    
    ///**************************** Free the memory ********************************///
    
    FREE_MATRIX_DOUBLE(Profile_Data, 0, 0);
    
    FREE_MATRIX_DOUBLE(Result, 0, 0);
    
    FREE_VECTOR_DOUBLE(Wavelength_Data, 0);
    
    
    ///**************************** Print the running time ****************************///
    
    clock_t timeb=clock();
    
    long hours=(timeb-timea)/CLOCKS_PER_SEC/3600;
    
    long minutes=((timeb-timea)/CLOCKS_PER_SEC-hours*3600)/60;
    
    double seconds=(timeb-timea)*1.0/CLOCKS_PER_SEC-hours*3600-minutes*60;
    
    printf("\nrunning time= %lu h %lu min %.2lf sec\n",hours,minutes,seconds);
    
    return 0;
}

