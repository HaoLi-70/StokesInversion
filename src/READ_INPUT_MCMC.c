
#include "READ_INPUT_MCMC.h"

/*--------------------------------------------------------------------------------*/

  /*######################################################################
   
    revision log:
    10 Sept. 2021.
   
   ######################################################################*/

/*--------------------------------------------------------------------------------*/

extern void Keywords_Conversion(STRUCT_KEYS Keywords[], STRUCT_INPUT *Input, \
    STRUCT_MPI *Mpi, STRUCT_ATOM *Atom){
  
    /*######################################################################
      Purpose:
        convert the keywords to the configuration.
      Record of revisions:
        10 Oct. 2022.
      Input parameters:
        Keywords, a structure saved input configuration.
        Mpi, a structure saved the Mpi information.
      Output parameters:
        Input, a structure saved the input information.
     ######################################################################*/
    
      //const char *routine_name = "Keywords_Conversion";
    
    
    int indx, nread;
    /*
    int num_read, tmp_array[20];
    char parameter[Key_Length];
    long len_tot;
    */
    
    indx = 0;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Verbose = true;
    }else{
      Input->Verbose = false;
    }
   
    indx = 1;
    sscanf(Keywords[indx].line, "%s", Input->Path2Profile);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n Path to the profile file : %s \n", Input->Path2Profile);
    }

    indx = 2;
    nread = sscanf(Keywords[indx].line, "%lf", &Atom->Lambda);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n The wavelength (A) : %.4f \n", Atom->Lambda);
    }

    indx = 3;
    nread = sscanf(Keywords[indx].line, "%lf", &Atom->Geffect);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n The effective lande factor : %.2f \n", Atom->Geffect);
    }

    indx = 4;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->GEMC_Simulation = true;
    }else{
      Input->GEMC_Simulation = false;
    }
    if(Input->Verbose && Mpi->Rank == 0){
      if(Input->GEMC_Simulation){
        fprintf(stderr, "\n GEMC simulation : YES \n");
      }else{
        fprintf(stderr, "\n GEMC simulation : No \n");
      }
    }

    indx = 5;
    nread = sscanf(Keywords[indx].line, "%d", &Input->Gemc_Nchains);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n Chains number for GEMC : %d \n", Input->Gemc_Nchains);
    }

    indx = 6;
    nread = sscanf(Keywords[indx].line, "%d", &Input->Dream_Nchains);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n Chains number for DREAM : %d \n", Input->Dream_Nchains);
    }

    indx = 7;
    nread = sscanf(Keywords[indx].line, "%d", &Input->Num_Gener);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n Generation number for DREAM : %d \n", Input->Num_Gener);
    }

    indx = 8;
    nread = sscanf(Keywords[indx].line, "%d", &Input->MaxNum_Pair);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n Max number for the chain pairs : %d \n", Input->MaxNum_Pair);
    }

    indx = 9;
    nread = sscanf(Keywords[indx].line, "%d", &Input->Num_Cr);
    if(Input->Verbose && Mpi->Rank == 0){
      fprintf(stderr, "\n Max number for the crossover probability : %d \n", Input->Num_Cr);
    }

    indx = 10;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Cr_Update = true;
    }else{
      Input->Cr_Update = false;
    }
    if(Input->Verbose && Mpi->Rank == 0){
      if(Input->Cr_Update){
        fprintf(stderr, "\n Update the crossover probability : YES \n");
      }else{
        fprintf(stderr, "\n Update the crossover probability : No \n");
      }
    }

    indx = 11;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Distribution_Update = true;
    }else{
      Input->Distribution_Update = false;
    }
    if(Input->Verbose && Mpi->Rank == 0){
      if(Input->Distribution_Update){
        fprintf(stderr, "\n Update the target distribution : YES \n");
      }else{
        fprintf(stderr, "\n Update the target distribution : No \n");
      }
    }

    indx = 12;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Sigma = true;
    }else{
      Input->Sigma = false;
    }
    if(Input->Verbose && Mpi->Rank == 0){
      if(Input->Sigma){
        fprintf(stderr, "\n compute 1 sigma region region : YES \n");
      }else{
        fprintf(stderr, "\n compute 1 sigma region region : No \n");
      }
    }

    indx = 13;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Output_Samples = true;
    }else{
      Input->Output_Samples= false;
    }
    if(Input->Verbose && Mpi->Rank == 0){
      if(Input->Output_Samples){
        fprintf(stderr, "\n Output the samples : YES \n");
      }else{
        fprintf(stderr, "\n Output the samples : No \n");
      }
    }


    indx = 14;
    String_to_Upper(Keywords[indx].line);
    if(strcmp(Keywords[indx].line, "YES") == 0){
      Input->Debug = true;
    }else{
      Input->Debug = false;
    }
    if(Input->Verbose && Mpi->Rank == 0){
      if(Input->Debug){
        fprintf(stderr, "\n Debug mode : YES \n");
      }else{
        fprintf(stderr, "\n Debug mode : No \n");
      }
    }

    return;
}

/*--------------------------------------------------------------------------------*/

extern void RDINPUT(char Filename[], STRUCT_INPUT *Input, STRUCT_MPI *Mpi, \
    STRUCT_ATOM *Atom, STRUCT_OBSERVATION *Observation){
  
    /*######################################################################
      Purpose:
        Read the input file for the forbidden line calculation.
      Record of revisions:
        10 Sept. 2021.
      Input parameters:
        Filename[], the input file.
        Mpi, a structure saved the Mpi information.
      Output parameters:
        Input, a structure saved the input information.
     ######################################################################*/
    
    const char *routine_name = "RDINPUT";
    
    FILE *fa;
    fa = fopen(Filename, "r");    
    int Num_Keywords = 12;

    STRUCT_KEYS Keywords[] ={
      {"verbose", "YES", false, false},  //0
      {"path2profile", "", false, true}, //1
      {"lambda", "", false, true},  //2
      {"geffect", "", false, true}, //3
      {"gemc_simulation", "NO", false, false}, //4
      {"gemc_nchains", "200", false, false}, //5
      {"dream_nchains", "40", false, false}, //6
      {"dream_ngener", "5000", false, false}, //7
      {"max_pair", "3", false, false}, //8
      {"ncr", "3", false, false}, //9
      {"cr_update", "YES", false, false}, //10
      {"distribution_update", "NO", false, false}, //11
      {"1sigma", "NO", false, false}, //12
      {"output_samples", "NO", false, false}, //13
      {"debug", "NO", false, false}, //14

    };
    Num_Keywords = sizeof(Keywords)/sizeof(STRUCT_KEYS);
    
    char lines[Max_Line_Length], key[Key_Length], *p;
    int len, i;
    long len_tot;
    bool neglect = true;
    
    while (Read_line(lines, fa) > 0 ){
      len_tot = strlen(lines);
      len = Indx_Char(lines, '=', 1);
      if(len > 1){
        len_tot -= len;
        p = lines+len;
        String_Copy(key, lines, len-1, 1);
        Trim(key, 3);
        
        neglect = true;
        for (i=0; i<Num_Keywords; i++){
          if(strcmp(key,Keywords[i].keyword)==0){
            if(Keywords[i].Set == true && (Mpi->Rank) == 0){
              Error(enum_warning, routine_name, "aa");
            }
            String_Copy(Keywords[i].line, p, len_tot, true);
            Keywords[i].Set = true;
            if(Mpi->Debug  && (Mpi->Rank) == 0){
              fprintf(stderr," %d %s %s \n",i, Keywords[i].keyword, Keywords[i].line);
            }
            neglect = false;
            break;
          }
        }
        if(Input->Verbose && Mpi->Rank == 0 && neglect){
          fprintf(stderr, "Warning : Neglect key words %s \n",key);
          fprintf(stderr, "The whole line is %s \n",lines);
        }
      }
    }
    
    Keywords_Conversion(Keywords, Input, Mpi, Atom);

    RDOBSERVATION(Input->Path2Profile, Observation);

    //Init_Atom(Atom);

    fclose(fa);
    
    return;
}

/*--------------------------------------------------------------------------------*/

extern int RDOBSERVATION(char Filename[], STRUCT_OBSERVATION *Observation){

    /*######################################################################
      Purpose:
        read the observe profile.
      Record of revisions:
        10 Oct. 2022.
      Input parameters:
        Filename, the path to the file saved the observed profile.
      Output parameters:
        Observation, a structure saved observed profile.
     ######################################################################*/

    FILE *fa;
    fa = fopen(Filename, "r");
    int i;

    fscanf(fa, "%d ", &Observation->Num_Lam);

    Observation->Lambda = (double *)VECTOR(0, Observation->Num_Lam-1, \
        enum_dbl, false);
    Observation->Profile_Data = (double **)MATRIX(0, 3, 0, \
        Observation->Num_Lam-1, enum_dbl, false);
    Observation->Stokes = (double **)MATRIX(0, 3, 0, \
        Observation->Num_Lam-1, enum_dbl, false);
    Observation->Sigma = (double **)MATRIX(0, 3, 0, \
        Observation->Num_Lam-1, enum_dbl, false);

    for(i =0 ; i<Observation->Num_Lam; i++){
      fscanf(fa, "%lf %lf %lf %lf %lf ", &Observation->Lambda[i], \
          &Observation->Profile_Data[0][i], &Observation->Profile_Data[1][i], \
          &Observation->Profile_Data[2][i], &Observation->Profile_Data[3][i]);
      Observation->Sigma[0][i] = Observation->Profile_Data[0][i]*1e-3;
      Observation->Sigma[1][i] = Observation->Profile_Data[0][i]*1e-3;
      Observation->Sigma[2][i] = Observation->Profile_Data[0][i]*1e-3;
      Observation->Sigma[3][i] = Observation->Profile_Data[0][i]*1e-3;
    }
   
    Observation->weights_flag = false;
    
    fclose(fa);
    CONTROL();

    return 0;
}

/*--------------------------------------------------------------------------------*/
