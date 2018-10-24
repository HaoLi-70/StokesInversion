#include <stdio.h>

#include <math.h>

#include <time.h>

#include <mpi.h>

#include "INIT.h"

#include "PATH.h"

#include "ALLOCATION.h"

#include "INVERSION.h"

#include "FITS_FILE.h"


int main(int argc, char *argv[]) {
    
    
    ///****************** Begin to calculate the Running Time **********///

    clock_t timea=clock();
    

    ///****************** Initialize the MPI Execution Environment **********///

    MPI_Init(&argc,&argv);
    
    
    ///****************** Initialize the Files Information **********///

    int Wav_begin, Wav_end, Pixel_begin, Pixel_end, Num_files;
    
    Init_file(&Wav_begin, &Wav_end, &Pixel_begin, &Pixel_end, &Num_files);
    
    int Num_Wav=Wav_end-Wav_begin;
    
    int Num_Pixel=Pixel_end-Pixel_begin;
    
    
    ///****************** Initialize the Wavelength **********///

    double *Wavelength_Data;
    
    Wavelength_Data=VECTOR_DOUBLE(0, Num_Wav-1);
    
    Init_Wavelength(Wavelength_Data, Num_Wav);

    
    ///****************** Memory Allocation **********///

    double **Result;
    
    Result=MATRIX_DOUBLE(0, 11, 0, 1);
    
    double **Profile_Data;
    
    Profile_Data=MATRIX_DOUBLE(0, Num_Wav-1, 0, 3);
    
    
    ///****************** Determines the rank of the calling process **********///

    int myid, numprocs;
    
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    
    
    ///****************** Parameters for the Inversion Procedure **********///
    /* If Flag_Correlation=1, save the correlation coefficients.
       If Flag_Histdata=1, save the hist data.  */
    
    int Flag_Correction=0, Flag_Histdata=0;

    double Accept_Rates;
    
    int i, j, k, m, idnumb;

    int LL=Num_Pixel/numprocs+1;
    
    if (LL*(numprocs-1)>Num_Pixel)
        
        LL--;
    
    if (numprocs<2) {
        
        fprintf(stderr,"Hello World! Process %d \n", numprocs);
        
        exit(0);
        
    }
    
    
    ///****************************** Master Process ******************************///
    
    if (myid==0) {
        
        ///****************************** Memory allocation for the Master Process ******************************///

        double **Fits_Data;
        
        Fits_Data=MATRIX_DOUBLE(0, Num_Pixel-1, 0, Num_Wav*4-1);
        
        double **INV_MAP;
        
        INV_MAP=MATRIX_DOUBLE(0, Num_Pixel-1, 0, 23);
        
        char **Name_fits;
        
        Name_fits=MATRIX_CHAR(0, Num_files-1, 0, 180);
        
        char **Res_fits;
        
        Res_fits=MATRIX_CHAR(0, Num_files-1, 0, 180);
        

        ///****************************** Determines the names of each file ******************************///

        filelist(Num_files, Name_fits, Res_fits);

        FILE *File;
        
        File=fopen(Path_RESULT, "w");
        
        MPI_Status status;

        long naxes[2];
        
        naxes[0]=24;
        
        naxes[1]=Num_Pixel;
        
        
        for (k=0; k<Num_files; k++) {
            
            ///****************************** Read the fits files ******************************///

            Fitsread(Name_fits[k], Wav_begin, Wav_end, Pixel_begin, Pixel_end, Fits_Data);
            
            fprintf(stderr, "%s \n\n",Name_fits[k]);

            
            ///****************************** Send the fits data to the slave process ******************************///

            for (idnumb=1; idnumb<numprocs; idnumb++)
                
                MPI_Send(Fits_Data[LL*(idnumb-1)], LL*Num_Wav*4, MPI_DOUBLE, idnumb, k, MPI_COMM_WORLD);
            
            for (m=LL*(numprocs-1); m<Num_Pixel; m++) {
                
                for (i=0; i<Num_Wav; i++) {
                    
                    for (j=0; j<4; j++)
                        
                        Profile_Data[i][j]=Fits_Data[m][i+Num_Wav*j];
                    
                    Profile_Data[i][0] *=2;
                    
                }
                
                ///****************************** Inversion of the profile ******************************///

                Inversion(Profile_Data, Wavelength_Data, Num_Wav, Flag_Correction, Flag_Histdata, Result, &Accept_Rates);
                
                fprintf(stderr, " File Number = %d Pixel Number = %d %e\n",k,m,Accept_Rates);

                for (i=0; i<12; i++) {
                    
                    INV_MAP[m][i*2]=Result[i][0];
                    
                    INV_MAP[m][i*2+1]=Result[i][1];
                    
                }
                
            }
            
            ///****************************** Receive the results form the slave process ******************************///

            for (idnumb=1; idnumb<numprocs; idnumb++)
                
                MPI_Recv(INV_MAP[LL*(idnumb-1)], LL*24, MPI_DOUBLE, idnumb, k+Num_files, MPI_COMM_WORLD, &status);

            
            ///****************************** print the results ******************************///
            
            Fitswrite(Res_fits[k], naxes, INV_MAP);
            
        }
        
        fclose(File);
        
        ///**************************** Free the memory ********************************///
        
        FREE_MATRIX_CHAR(Name_fits, 0, 0);
        
        FREE_MATRIX_DOUBLE(Fits_Data, 0, 0);
        
        FREE_MATRIX_DOUBLE(INV_MAP, 0, 0);
        
        
        ///**************************** Print the running time ****************************///
        
        clock_t timeb=clock();
        
        long hours=(timeb-timea)/CLOCKS_PER_SEC/3600;
        
        long minutes=((timeb-timea)/CLOCKS_PER_SEC-hours*3600)/60;
        
        double seconds=(timeb-timea)*1.0/CLOCKS_PER_SEC-hours*3600-minutes*60;
        
        printf("\nrunning time= %lu h %lu min %.2lf sec\n",hours,minutes,seconds);
        
    }
    ///****************************** Slave Process ******************************///
    else{
        
        ///****************************** Memory allocation for the Slave Process ******************************///
        
        double **FITS_DATA;
        
        FITS_DATA=MATRIX_DOUBLE(0, LL-1, 0, Num_Wav*4-1);
        
        double **INV_MAP;
        
        INV_MAP=MATRIX_DOUBLE(0, LL-1, 0, 23);
        
        MPI_Status status;
        
        for (k=0; k<Num_files; k++) {
            
            ///****************************** Receive the profiles form the master process ******************************///

            MPI_Recv(FITS_DATA[0], LL*Num_Wav*4, MPI_DOUBLE, 0, k, MPI_COMM_WORLD, &status);
            
            for (m=0; m<LL; m++) {
         
                for (i=0; i<Num_Wav; i++) {
                    
                    for (j=0; j<4; j++)
                        
                        Profile_Data[i][j]=FITS_DATA[m][i+Num_Wav*j];
                    
                    Profile_Data[i][0] *=2;
                    
                }
                
                ///****************************** Inversion of the profile ******************************///

                Inversion(Profile_Data, Wavelength_Data, Num_Wav, Flag_Correction, Flag_Histdata, Result, &Accept_Rates);
                
                fprintf(stderr, " File Number = %d Pixel Number = %d %e\n",k,m+LL*(myid-1),Accept_Rates);

         
                for (i=0; i<12; i++) {
                    
                    INV_MAP[m][i*2]=Result[i][0];
                    
                    INV_MAP[m][i*2+1]=Result[i][1];
                    
                }
                
            }
            
            MPI_Send(INV_MAP[0], LL*24, MPI_DOUBLE, 0, k+Num_files, MPI_COMM_WORLD);
            
        }
        
        ///**************************** Free the memory ********************************///

        FREE_MATRIX_DOUBLE(FITS_DATA, 0, 0);
        
        FREE_MATRIX_DOUBLE(INV_MAP, 0, 0);

    }
    
    ///**************************** Free the memory ********************************///

    FREE_VECTOR_DOUBLE(Wavelength_Data, 0);
    
    FREE_MATRIX_DOUBLE(Profile_Data, 0, 0);
    
    FREE_MATRIX_DOUBLE(Result, 0, 0);

    
    MPI_Finalize();

    return 0;
}

