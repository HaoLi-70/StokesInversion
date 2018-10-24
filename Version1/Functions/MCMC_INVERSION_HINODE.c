#include <stdio.h>

#include <math.h>

#include <time.h>

#include <fitsio.h>

#include "ALLOCATION.h"

#include "INIT.h"

#include "PATH.h"

#include "INVERSION.h"

#include "FITS_FILE.h"

int main(int argc, const char * argv[]) {
    
    
    ///****************** Begin to calculate the running time **********///
    clock_t timea=clock();
    
    
    ///****************** Initialize the Files Information **********///
    
    int Wav_begin, Wav_end, Pixel_begin, Pixel_end, Num_files;
    
    Init_file(&Wav_begin, &Wav_end, &Pixel_begin, &Pixel_end, &Num_files);
    
    int Num_Wav=Wav_end-Wav_begin;
    
    int Num_Pixel=Pixel_end-Pixel_begin;

    
    ///****************************** Memory allocation ******************************///
    
    double **Fits_Data;
    
    Fits_Data=MATRIX_DOUBLE(0, Pixel_end-Pixel_begin-1, 0, Num_Wav*4-1);
    
    double **Profile_Data;
    
    Profile_Data=MATRIX_DOUBLE(0, Num_Wav-1, 0, 3);
    
    double *Wavelength_Data;
    
    Wavelength_Data=VECTOR_DOUBLE(0, Num_Wav-1);
    
    double **Result;
    
    Result=MATRIX_DOUBLE(0, 11, 0, 1);
    
    char **Name_fits;
    
    Name_fits=MATRIX_CHAR(0, Num_files-1, 0, 180);
    
    char **Res_fits;
    
    Res_fits=MATRIX_CHAR(0, Num_files-1, 0, 180);
    
    double **INV_MAP;
    
    INV_MAP=MATRIX_DOUBLE(0, Num_Pixel-1, 0, 23);
    
    
    ///****************** Initialize the Wavelength **********///
    
    Init_Wavelength(Wavelength_Data, Num_Wav);

    
    ///****************** Parameters for the Inversion Procedure **********///
    /* If Flag_Correlation=1, save the correlation coefficients.
     If Flag_Histdata=1, save the hist data.  */
    
    int Flag_Correction=0, Flag_Histdata=0;

    double Accept_Rates;

    int i, j, k, m;
    
    long naxes[2];
    
    naxes[0]=24;
    
    naxes[1]=Num_Pixel;
    
    
    ///****************************** Determines the names of each file ******************************///
    
    filelist(Num_files, Name_fits,Res_fits);
    
    
    for (k=0; k<Num_files; k++) {
        
        ///****************************** Input the profile data ******************************///
        Fitsread(Name_fits[k], Wav_begin, Wav_end, Pixel_begin, Pixel_end, Fits_Data);
        
        fprintf(stderr, "\n%s \n",Name_fits[k]);
        
        for (m=0; m<Num_Pixel; m++) {
            
            printf("\n Pixle = %d \n",m);
            
            for (i=0; i<Num_Wav; i++) {
                
                for (j=0; j<4; j++)
                    
                    Profile_Data[i][j]=Fits_Data[m][i+Num_Wav*j];
                
                Profile_Data[i][0] *=2;
            
            }
            
            Inversion(Profile_Data, Wavelength_Data, Num_Wav, Flag_Correction, Flag_Histdata, Result, &Accept_Rates);
            
            for (i=0; i<12; i++) {
                
                INV_MAP[m][i*2]=Result[i][0];
                
                INV_MAP[m][i*2+1]=Result[i][1];
                
            }
            
        }
        
        Fitswrite(Res_fits[k], naxes, INV_MAP);
        
    }
    
    
    
    ///**************************** Free the memory ********************************///
   
    FREE_MATRIX_CHAR(Name_fits, 0, 0);
    
    FREE_MATRIX_DOUBLE(Profile_Data, 0, 0);
    
    FREE_MATRIX_DOUBLE(Result, 0, 0);
    
    FREE_VECTOR_DOUBLE(Wavelength_Data, 0);
    
    FREE_MATRIX_DOUBLE(INV_MAP, 0, 0);

    
    ///**************************** Print the running time ****************************///
    
    clock_t timeb=clock();
    
    long hours=(timeb-timea)/CLOCKS_PER_SEC/3600;
    
    long minutes=((timeb-timea)/CLOCKS_PER_SEC-hours*3600)/60;
    
    double seconds=(timeb-timea)*1.0/CLOCKS_PER_SEC-hours*3600-minutes*60;
    
    printf("\nrunning time= %lu h %lu min %.2lf sec\n",hours,minutes,seconds);
    
    return 0;
}

