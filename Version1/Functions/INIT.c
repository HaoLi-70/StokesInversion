#include "INIT.h"


void Init_Wav(int *Wav_Num){
    
    /***********************************************************************************
     Purpose:
     Initial the number of wavelength points. (only used for .dat file read)
     Modified:
     22 Jun. 2018, Hao Li
     Output parameters:
     Wav_Num, the number of wavelength points.
     ***********************************************************************************/
    
    *Wav_Num=51;
    
    return;
}

void Init_Atom(TRANSITION_INFORMATION *Atom){
    
    /***********************************************************************************
     Purpose:
     Initial the atomic information.
     Modified:
     22 Jun. 2018, Hao Li
     Output parameters:
     Atom, return the atomic information
     ***********************************************************************************/
    
    (*Atom).Lambda=6302.4932e-10;
    
    (*Atom).Ju=0.0;
    
    (*Atom).Lu=2.0;
    
    (*Atom).Su=2.0;
    
    (*Atom).Jl=1.0;
    
    (*Atom).Ll=1.0;
    
    (*Atom).Sl=2.0;
    
    (*Atom).Gu=Gfactor((*Atom).Ju, (*Atom).Lu, (*Atom).Su);
    
    (*Atom).Gl=Gfactor((*Atom).Jl, (*Atom).Ll, (*Atom).Sl);
    
    (*Atom).Geffect=Geffect((*Atom).Gu, (*Atom).Gl, (*Atom).Ju, (*Atom).Jl);
    
    return;
    
}

void Init_file(int *Wav_begin, int *Wav_end, int *Pixel_begin, int *Pixel_end, int *Num_file){
    
    /***********************************************************************************
     Purpose:
     Initial the inversion regions of the files (HINODE/SP).
     Modified:
     17 Oct. 2018, Hao Li
     Output parameters:
     Num_file, total number of the files.
     Wav_begin, the first pixel in the wavelength axis used in the inversion.
     Wav_end, the last pixel (not include) in the wavelength axis used in the inversion process.
     Pixel_begin, the first pixel along the slit used in the inversion.
     Pixel_end, the last pixel (not include) along the slit used in the inversion process.
     ***********************************************************************************/
    
    *Num_file=3;
    
    *Wav_begin=55;
    
    *Wav_end=95;
    
    *Pixel_begin=95;
    
    *Pixel_end=226;
    
    return;
    
}

void filelist(int Num_file, char **Name_fits, char **Res_fits){
    
    /***********************************************************************************
     Purpose:
     Initial the fits name used in the inversion process (HINODE/SP).
     Modified:
     17 Oct. 2018, Hao Li
     Input parameters:
     Num_file, total number of the files.
     Output parameters:
     Name_fits, the fits name used in the inversion process (read from a list file).
     Res_fits, the fits name used to save the inversion result.
     ***********************************************************************************/
    
    FILE *File;
    
    int i;
    
    char str[30],sst[180]=Path_fits, st[4]="Re_",stt[4]="!";
    
    File=fopen(Path_List, "r");
    
    for (i=0; i<Num_file; i++) {
        
        fscanf(File, "%s ",str);
        
        strcat(Name_fits[i],sst);
        
        strcat(Name_fits[i],str);
        
        strcat(Res_fits[i],stt);
        
        strcat(Res_fits[i],sst);
        
        strcat(Res_fits[i],st);
        
        strcat(Res_fits[i],str);
    
    }
    
    fclose(File);

    return;
}

void Init_Wavelength(double *Wavelength_Data, int Num_Wav){
    
    /***********************************************************************************
     Purpose:
     Initial the wavelength of each point.
     Modified:
     22 Jun. 2018, Hao Li
     Input parameters:
     Num_Wav, number of profile points.
     Output parameters:
     Wavelength_Data, return the wavelength of each point
     ***********************************************************************************/
    
    int i=0;
    
    for (i=0; i<Num_Wav; i++)
    
        Wavelength_Data[i]=6302e-10+(i-0.5)*0.021549e-10;
    
    return;
    
}

void Init_Parameter(int *Num_Parameter, int *Num_Chain_Gemc, int *Num_Chain_Dream, int *Begin_Generation, int *Num_Generation, int *Num_Cr, int *Pair_Number, int *Indx_Bounds_Gemc, int *Indx_Bounds_Dream, _Bool *Cr_adapt_flag){
    
    /***********************************************************************************
     Purpose:
     Initial the parameters of MCMC simulation.
     Modified:
     22 Jun. 2018, Hao Li
     Output parameters:
     Number_Parameter, the total number of parameters. (now only 9 is accepted)
     Num_Chain_Gemc, the total number of chains for GEMC simulation.
     Num_Chain_Dream, the total number of chains for DREAM simulation.
     Begin_Generation, end of the burn in period.
     Num_Generation, the total number of gererations.
     Num_Cr, the total number of CR values.
     Pair_Number, the max number of chain pairs used to generate the jump.
     Indx_Bounds_Gemc, the index of bounds condition for GEMC simulation. 0, the fold bounds; 1, the reflect bounds; 2, set to the bound value
     Indx_Bounds_Dream, the index of bounds condition for DREAM simulation. 0, the fold bounds; 1, the reflect bounds; 2, set to the bound value
     Cr_adapt_flag, the Cr probability adapt flag.
     ***********************************************************************************/
    
    *Num_Parameter=9;
    
    *Num_Chain_Gemc=100;
    
    *Num_Chain_Dream=36;
    
    *Begin_Generation=2000;
    
    *Num_Generation=7000;
    
    *Num_Cr=3;
    
    *Pair_Number=3;
    
    *Indx_Bounds_Gemc=1;
    
    *Indx_Bounds_Dream=0;
    
    *Cr_adapt_flag=1;
    
    return;
    
}


void Init_Bounds(double **Profile_Data, double *Wavelength_Data, int Num_Wav, int Num_Parameter, float Lambda, double **Parameter_Bounds){
    
    /***********************************************************************************
     Purpose:
     Initial the boud values of each parameter.
     Modified:
     22 Jun. 2018, Hao Li
     Input parameters:
     Profile_Data[][], the Stokes profiles.
     Wavelength_Data[Num_Wav], wavelength of each profile points.
     Num_Wav, number of profile points.
     Num_Parameter, the total number of parameters. (now only 9 is accepted)
     Lambda, wavelength of the line center.
     Output parameters:
     Parameter_Bounds[][], the boud values of each parameter.
     ***********************************************************************************/
    
    const double Con_C=299792458.0;
    
    int i;
    
    double MAX=Profile_Data[0][0];
    
    double sum=0, weight=0;
    
    //S0+S1: MAX value in Stokes I
    for (i=1; i<Num_Wav; i++) {
        
        //MAX=Profile_Data[i][1]>MAX?Profile_Data[i][1]:MAX;
        if (Profile_Data[i][0]>MAX)
            
            MAX=Profile_Data[i][0];
        
    }
    
    Parameter_Bounds[5][0]=MAX*0.6;
    
    Parameter_Bounds[5][1]=MAX*1.4;
    
    //Velocity: we use the COG method in Stokes I
    
    for (i=0; i<Num_Wav; i++) {
        
        sum += Wavelength_Data[i]*(MAX-Profile_Data[i][0]);
        
        weight += MAX-Profile_Data[i][0];
        
    }
    
    double Delta_Velocity;
    
    Delta_Velocity=1026;//=(Profile_Data[1][0]-Profile_Data[0][0])/Lambda*Con_C;
    
    double Velocity=(sum/weight-Lambda)/Lambda*Con_C;
    
    Parameter_Bounds[3][0]=Velocity-Delta_Velocity*2;
    
    Parameter_Bounds[3][1]=Velocity+Delta_Velocity*2;
    
    Parameter_Bounds[0][0]=0;
    
    Parameter_Bounds[0][1]=4000;
    
    Parameter_Bounds[1][0]=0;
    
    Parameter_Bounds[1][1]=3.1415926;
    
    Parameter_Bounds[2][0]=0;
    
    Parameter_Bounds[2][1]=3.1415926;
    
    Parameter_Bounds[4][0]=1e-12;
    
    Parameter_Bounds[4][1]=6e-12;
    
    Parameter_Bounds[6][0]=0;
    
    Parameter_Bounds[6][1]=70;//90

    Parameter_Bounds[7][0]=0;
    
    Parameter_Bounds[7][1]=0.8;
    
    Parameter_Bounds[8][0]=0;
    
    Parameter_Bounds[8][1]=1;//1.2
    
    return;
    
}
