#include "INVERSION.h"


int Inversion(double **Profile_Data, double *Wavelength_Data, int Num_Wav, int Flag_Correlation, int Flag_Histdata, double **Result, double *Accept_Rates){
    
    /***********************************************************************************
     Purpose:
     Inversion of Stokes Profile.
     Modified:
     17 Oct 2018, Hao Li
     Input parameters:
     Profile_Data[][4], the Stokes profiles used to inverse the magnetic field.
     Wavelength_Data[Num_Wav], wavelength of each profile points.
     Num_Wav, number of profile points.
     Flag_Correlation,
     Flag_Histdata,
     Result[][], the result of the inversion.
     ***********************************************************************************/
    
    
    ///****************************** Initial the parameter ******************************///
    
    int Num_Parameter, Num_Chain_Gemc, Num_Chain_Dream, Begin_Generation=2000, Num_Generation=7000, Num_Cr=3, Pair_Number=3, Indx_Bounds_Gemc=1, Indx_Bounds_Dream=0;
    
    _Bool Cr_adapt_flag=1;
    
    Init_Parameter(&Num_Parameter, &Num_Chain_Gemc, &Num_Chain_Dream, &Begin_Generation, &Num_Generation, &Num_Cr, &Pair_Number, &Indx_Bounds_Gemc, &Indx_Bounds_Dream, &Cr_adapt_flag);
    
    int Chain_bestfit, i, j, k, m, n, Convergence=0, ITER=0;
    
    _Bool flag;
    
    FILE *File;
    
    
    ///****************************** Memory allocation ******************************///
    
    double **Parameter_Bounds;
    Parameter_Bounds=MATRIX_DOUBLE(0, Num_Parameter-1, 0, 1);
    
    double **RESULT_GEMC;
    RESULT_GEMC=MATRIX_DOUBLE(0, Num_Chain_Gemc-1, 0, Num_Parameter-1);
    
    double *LIKELIHOOD;
    LIKELIHOOD=VECTOR_DOUBLE(0, Num_Chain_Gemc-1);
    
    double *LIKELIHOOD_SORT;
    LIKELIHOOD_SORT=VECTOR_DOUBLE(0, Num_Chain_Gemc-1);
    
    double **RESULT_DREAM;
    RESULT_DREAM=MATRIX_DOUBLE(0, Num_Chain_Dream-1, 0, Num_Parameter*Num_Generation-1);
    
    double **LIKELIHOOD_DREAM;
    LIKELIHOOD_DREAM=MATRIX_DOUBLE(0, Num_Generation-1, 0, Num_Chain_Dream-1);
    
    double *Best_fit;
    Best_fit=VECTOR_DOUBLE(0, Num_Parameter);
    
    double *STD, *Mean;
    
    STD=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    Mean=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    double **MAG;
    
    MAG=MATRIX_DOUBLE(0, Num_Chain_Dream-1, 0, 3*Num_Generation-1);
    
    TRANSITION_INFORMATION *ATOM;
    
    ATOM=malloc(sizeof(TRANSITION_INFORMATION));
    
    ///****************************** Initial the atomic information ******************************///
    
    Init_Atom(ATOM);
    
    ///****************************** Initial the bounds ******************************///

    Init_Bounds(Profile_Data, Wavelength_Data, Num_Wav, Num_Parameter, (*ATOM).Lambda, Parameter_Bounds);
    
  
    ///****************************** GEMC simulation ******************************///
    
RE: Chain_Ini(Num_Chain_Gemc, Num_Wav, Wavelength_Data, Profile_Data, Num_Parameter, Parameter_Bounds, *ATOM, RESULT_GEMC, LIKELIHOOD);
    
    if (ITER>0) {
        
        if (Mean[2]>3.1415926)
            
            Mean[2]-=3.1415926;
        
        if (Mean[2]<0)
            
            Mean[2]+=3.1415926;
        
        for (i=0; i<Num_Parameter; i++)
            
            RESULT_GEMC[0][i]=Mean[i];
        
        LIKELIHOOD[0]=Likelihood_Log(RESULT_GEMC[0], Num_Wav, Wavelength_Data, Profile_Data, *ATOM);
        
    }
    
    Chain_bestfit=GEMC(Num_Chain_Gemc, Num_Wav, Wavelength_Data, Profile_Data, Num_Parameter, Parameter_Bounds, *ATOM, RESULT_GEMC, LIKELIHOOD, Indx_Bounds_Gemc);
    
    
    ///***************** Modify the bounds of the azimuth if the azimuth is close to the bounds *****************///

    if (RESULT_GEMC[Chain_bestfit][2]>2.6||RESULT_GEMC[Chain_bestfit][2]<0.5) {
        
        Parameter_Bounds[2][0]=RESULT_GEMC[Chain_bestfit][2]-3.1415926/2;
        
        Parameter_Bounds[2][1]=RESULT_GEMC[Chain_bestfit][2]+3.1415926/2;
        
    }
    
    
    ///***************************** Sort the likelihoode ***************************///
    ///*************  used as the first genneration of DREAM simulation *************///
    
    for (i=0; i<Num_Chain_Gemc; i++)
        
        LIKELIHOOD_SORT[i]=LIKELIHOOD[i];
    
    HPSORT(Num_Chain_Gemc, LIKELIHOOD_SORT);
    
    i=0;
    
    for (j=0; j<Num_Chain_Dream&&i<Num_Chain_Gemc;) {
        
        if (LIKELIHOOD[i]>=LIKELIHOOD_SORT[Num_Chain_Dream-1]) {
            
            for (k=0; k<Num_Parameter; k++)
                
                RESULT_DREAM[j][k]=RESULT_GEMC[i][k];
            
            LIKELIHOOD_DREAM[0][j]=LIKELIHOOD[i];
                
            Bounds_Enforce(RESULT_DREAM[j], Num_Parameter, Parameter_Bounds, 2);
                
            j++;
                
        }
            
        i++;
            
    }
    
    for (i=0; i<Num_Parameter; i++)
        
        Best_fit[i]=RESULT_GEMC[Chain_bestfit][i];
    
    Best_fit[Num_Parameter]=LIKELIHOOD[Chain_bestfit];
    
    
    ///****************************** DREAM simulation ******************************///
    
    *Accept_Rates=Dream(Num_Wav, Wavelength_Data, Profile_Data, Num_Parameter, Num_Generation, Num_Chain_Dream, Num_Cr, Pair_Number, Parameter_Bounds, Indx_Bounds_Dream, *ATOM, Cr_adapt_flag, RESULT_DREAM, LIKELIHOOD_DREAM, Best_fit, &Convergence);
    
    
    ///************ Calculate the mean values and the standard deviation ***************///

    MCMC_STD(Begin_Generation, Num_Generation, Num_Chain_Dream, Num_Parameter, RESULT_DREAM, STD, Mean);
    
    ///*********** Do the MCMC simulation again if the standard deviation is too large ************///
    
    flag=(Convergence)&&(STD[0]<8);
        
    if (!flag&&ITER<2) {
        
        ITER++;
        
        Parameter_Bounds[2][0]=0;
        
        Parameter_Bounds[2][1]=3.1415926;
        
        goto RE;
        
    }

    
    ///******************** Calculate the correlation coefficients *********************///
    
    if (Flag_Correlation==1) {
        
        double **Correlation;
        
        Correlation=MATRIX_DOUBLE(0, Num_Parameter-1, 0, Num_Parameter-1);
        
        File=fopen(Path_Correlation, "w");
        
        for (m=0; m<Num_Parameter; m++)
            
            for (n=0; n<Num_Parameter; n++)
                
                Correlation[m][n]=0;
            
        
        for (m=0; m<Num_Parameter; m++)
            
            for (n=0; n<Num_Parameter; n++)
                
                for (i=Begin_Generation; i<Num_Generation; i++)
                    
                    for (j=0; j<Num_Chain_Dream; j++)
                        
                        Correlation[m][n] += RESULT_DREAM[j][m+Num_Parameter*i]*RESULT_DREAM[j][n+Num_Parameter*i];

        for (i=0; i<Num_Parameter; i++)
            
            STD[i]=sqrt(STD[i]*STD[i]*((Num_Generation-Begin_Generation)*Num_Chain_Dream-1)/((Num_Generation-Begin_Generation)*Num_Chain_Dream));
        
        for (m=0; m<Num_Parameter; m++) {
            
            for (n=m+1; n<Num_Parameter; n++) {
                
                Correlation[m][n] = (Correlation[m][n]/Num_Chain_Dream/(Num_Generation-Begin_Generation)-Mean[m]*Mean[n])/STD[m]/STD[n];
                
                fprintf(File,"%.5f ",Correlation[m][n]);
                
            }
            
            fprintf(File,"\n");
            
        }
        
        FREE_MATRIX_DOUBLE(Correlation, 0, 0);

    }
    
    
    
    ///**************************** Save the hist data ********************************///
    
    if (Flag_Histdata==1) {
        
        File=fopen(Path_Histdata, "w");
        
        for (i=Begin_Generation; i<Num_Generation; i++) {
            
            for (j=0; j<Num_Chain_Dream; j++) {
                
                for (k=0; k<Num_Parameter; k++)
                    
                    fprintf(File,"%e ",(RESULT_DREAM[j][k+Num_Parameter*i]));
                
                fprintf(File,"\n");
                
            }
            
        }
        
        //fprintf(stderr,"\n Number of Sample %d \n",(Num_Generation-Begin_Generation)*Num_Chain_Dream);
        
        fclose(File);
        
    }
    
    
    ///***************** Modify the azimuth to the range from 0 to pi *********************///

    if (Mean[2]>3.1415926)
        
        Mean[2]-=3.1415926;
    
    if (Mean[2]<0)
        
        Mean[2]+=3.1415926;
    
    if (Best_fit[2]>3.1415926)
        
        Best_fit[2]-=3.1415926;
    
    if (Best_fit[2]<0)
        
        Best_fit[2]+=3.1415926;
    
    for (i=0; i<Num_Parameter; i++) {
        
        Result[i][0]=Mean[i];
        
        Result[i][1]=STD[i];
        
    }
    
    
    ///***************** Calculate the three components of the magnetic field vector *********************///

    for (i=0; i<Num_Chain_Dream; i++) {
        
        for (j=0; j<Num_Generation; j++) {
            
            MAG[i][j*3]=RESULT_DREAM[i][Num_Parameter*j]*cos(RESULT_DREAM[i][Num_Parameter*j+1]);
            
            MAG[i][j*3+1]=RESULT_DREAM[i][Num_Parameter*j]*sin(RESULT_DREAM[i][Num_Parameter*j+1])*cos(2*RESULT_DREAM[i][Num_Parameter*j+2]);
            
            MAG[i][j*3+2]=RESULT_DREAM[i][Num_Parameter*j]*sin(RESULT_DREAM[i][Num_Parameter*j+1])*sin(2*RESULT_DREAM[i][Num_Parameter*j+2]);
            
        }
        
    }

    MCMC_STD(Begin_Generation, Num_Generation, Num_Chain_Dream, 3, MAG, STD, Mean);
    
    for (i=0; i<3; i++) {
        
        Result[i+Num_Parameter][0]=Mean[i];
        
        Result[i+Num_Parameter][1]=STD[i];
        
    }
 
    ///**************************** Free the memory ********************************///
    
    FREE_MATRIX_DOUBLE(Parameter_Bounds, 0, 0);
        
    FREE_MATRIX_DOUBLE(RESULT_GEMC, 0, 0);
        
    FREE_MATRIX_DOUBLE(RESULT_DREAM, 0, 0);
        
    FREE_MATRIX_DOUBLE(LIKELIHOOD_DREAM, 0, 0);
        
    FREE_VECTOR_DOUBLE(LIKELIHOOD, 0);
        
    FREE_VECTOR_DOUBLE(LIKELIHOOD_SORT, 0);
        
    FREE_VECTOR_DOUBLE(Best_fit, 0);
    
    FREE_VECTOR_DOUBLE(STD, 0);
        
    FREE_VECTOR_DOUBLE(Mean, 0);
        
    FREE_MATRIX_DOUBLE(MAG, 0, 0);
    
    free(ATOM);
    
    return Convergence;
    
}

