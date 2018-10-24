#include "MCMC.h"


double Dream(int Num_Wav, double *Wavelength_Data, double **Profile_Data, int Num_Parameter, int Num_Generation, int Num_Chain, int Num_Cr, int Pair_Number, double **Parameter_Bounds, int Indx_Bounds, TRANSITION_INFORMATION ATOM, _Bool Cr_adapt_flag, double **RESULT_DREAM, double **LIKELIHOOD_DREAM, double *Best_fit, int *Convergence){
    
    /***********************************************************************************
     Purpose:
     DREAM (DiffeRential Evolution Adaptive Metropolis) algorithm, used to inverse the magnetic field.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     Num_Wav, number of profile points.
     Wavelength_Data[Num_Wav], wavelength of each profile points.
     Profile_Data[][4], the Stokes profiles used to inverse the magnetic field.
     Num_Parameter, the total number of parameters. (now only 9 is accepted)
     Num_Generation, the total number of gererations.
     Num_Chain, the total number of chains.
     Num_Cr, the total number of CR values.
     Pair_Number, the max number of chain pairs used to generate the jump.
     Parameter_Bounds[][2], the bounds for each patameter.
     Indx_Bounds: the index of bounds condition. 0, the fold bounds; 1, the reflect bounds; 2, set to the bound value
     ATOM, a sturcture which saves Atomic information.
     Cr_adapt_flag, the Cr probability adapt flag.
     Output parameters:
     RESULT_DREAM[][], a matrix which saves MCMC data.
     LIKELIHOOD_DREAM[][], a matrix which saves likelihood of the MCMC simulation.
     Best_fit[], return the best fitting chain.
     Convergence, return the convergence status.
     reference:
     Vrugt, J. A., et al. 2009. International Journal of Nonlinear Science & Numerical Simulation, 10(3), 273-290
     ***********************************************************************************/
    
    // Set the random number generator and generate new seed
    long *idum;
    
    idum=(long *)malloc(sizeof(long));
    
    Random_Seed(idum);
    
    double *Cr, *Cr_Prob, *Cr_Delta;
    
    Cr=VECTOR_DOUBLE(1, Num_Cr);
    
    Cr_Prob=VECTOR_DOUBLE(1, Num_Cr);
    
    Cr_Delta=VECTOR_DOUBLE(1, Num_Cr);
    
    int *Cr_L;
    
    Cr_L=VECTOR_INT(1, Num_Cr);
    
    Ini_Cr(Num_Cr, Cr, Cr_Prob, Cr_Delta, Cr_L);
    
    double *DREAM_SAMPLE;
    
    DREAM_SAMPLE=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    int count_sample=0, count_accept=0;
    
    int Indx_Generation=1, Indx_Chain=0, Indx_Parameter, Indx_Cr, Num_Pair, Num_outlierchain=0;
    
    double TMP_LIKELIHOOD, Acceptrates;
    
    *Convergence=0;
    
    double *Statistics_R;
    
    Statistics_R=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    for (Indx_Generation=1; Indx_Generation<Num_Generation; Indx_Generation++) {
        
        for (Indx_Chain=0; Indx_Chain<Num_Chain; Indx_Chain++) {
            
            Num_Pair=Sample_PairNum(Pair_Number, idum);

            Indx_Cr=Dream_Sample(RESULT_DREAM, Num_Chain, Num_Parameter, Num_Pair, Indx_Chain, Indx_Generation, Num_Cr, Cr, Cr_Prob, Parameter_Bounds, Indx_Bounds, idum, DREAM_SAMPLE);
            
                count_sample++;
            
            TMP_LIKELIHOOD=Likelihood_Log(DREAM_SAMPLE, Num_Wav, Wavelength_Data, Profile_Data, ATOM);
            
            if (TMP_LIKELIHOOD>=LIKELIHOOD_DREAM[Indx_Generation-1][Indx_Chain]) {
                
                count_accept++;
                
                LIKELIHOOD_DREAM[Indx_Generation][Indx_Chain]=TMP_LIKELIHOOD;
                
                for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++) {
                    
                    RESULT_DREAM[Indx_Chain][Indx_Parameter+Num_Parameter*Indx_Generation]=DREAM_SAMPLE[Indx_Parameter];
                    
                }
                
                if (TMP_LIKELIHOOD>Best_fit[Num_Parameter]) {
                    
                    Best_fit[Num_Parameter]=TMP_LIKELIHOOD;
                    
                    for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++) {
                        
                        Best_fit[Indx_Parameter]=DREAM_SAMPLE[Indx_Parameter];
                        
                    }
                    
                }
                
            }else if (Ran1(idum)<exp(TMP_LIKELIHOOD-LIKELIHOOD_DREAM[Indx_Generation-1][Indx_Chain])) {
                
                count_accept++;
                
                LIKELIHOOD_DREAM[Indx_Generation][Indx_Chain]=TMP_LIKELIHOOD;
                
                for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++) {
                    
                    RESULT_DREAM[Indx_Chain][Indx_Parameter+Num_Parameter*Indx_Generation]=DREAM_SAMPLE[Indx_Parameter];
                    
                }
                
            }else{
                
                LIKELIHOOD_DREAM[Indx_Generation][Indx_Chain]=LIKELIHOOD_DREAM[Indx_Generation-1][Indx_Chain];
                
                for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++) {
                    
                    RESULT_DREAM[Indx_Chain][Indx_Parameter+Num_Parameter*(Indx_Generation)] =RESULT_DREAM[Indx_Chain][Indx_Parameter+Num_Parameter*(Indx_Generation-1)];
                    
                }
                
            }
            
            if (Indx_Generation<=0.1*Num_Generation&&Cr_adapt_flag&&!(*Convergence))
                
                Cr_dis_Update(RESULT_DREAM, Num_Chain, Num_Parameter, Indx_Chain, Indx_Generation, Indx_Cr, Cr_L, Cr_Delta);
            
        }
        
        if (!(*Convergence)) {
            
            if (Cr_adapt_flag&&Indx_Generation<=0.1*Num_Generation)
                
                Cr_Pro_Update(Num_Cr, Cr_Delta, Cr_L, Cr_Prob);
            
            if (Indx_Generation%50==0) {
                
                if (Indx_Generation<=0.1*Num_Generation) {
                    
                    Num_outlierchain=Rm_Outlierchain(Num_Chain, Indx_Generation, Num_Parameter, RESULT_DREAM, LIKELIHOOD_DREAM);
                    
                }
                
                PSRF(RESULT_DREAM, Indx_Generation, Num_Chain, Num_Parameter, Statistics_R);
                
                *Convergence=1;
                
                for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++) {
                    
                    if (Statistics_R[Indx_Parameter]>=1.2) {
                        
                        *Convergence=0;
                        
                        break;
                        
                    }
                    
                }
                //printf("Index of Generation %d outlierchain = %d convergence = %d \n",Indx_Generation, Num_outlierchain, *Convergence);
            }
        }//else if (Indx_Generation%500==0) printf("Index of Generation %d \n",Indx_Generation);
        
    }
    
    Acceptrates=(double)(count_accept)/count_sample;
    
    printf("Dream sample = %d accept = %d accept rates = %f convergence = %d\n",count_sample,count_accept,Acceptrates, *Convergence);
    
    // free the memory
    
    free(idum);
    
    FREE_VECTOR_DOUBLE(Cr, 1);
    
    FREE_VECTOR_DOUBLE(Cr_Delta, 1);
    
    FREE_VECTOR_DOUBLE(Cr_Prob, 1);
    
    FREE_VECTOR_INT(Cr_L, 1);
    
    FREE_VECTOR_DOUBLE(DREAM_SAMPLE, 0);
    
    FREE_VECTOR_DOUBLE(Statistics_R, 0);
    
    return Acceptrates;
    
}


int Dream_Sample(double **RESULT_DREAM, int Num_Chain, int Num_Parameter, int Num_Pair, int Indx_Chain, int Indx_Generation, int Num_Cr, double *Cr, double *Cr_Prob, double **Parameter_Bounds, int Indx_Bounds, long *idum, double *Dream_Sample){
    
    /***********************************************************************************
     Purpose:
     Initialize the CR (crossover probability) values.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     RESULT_DREAM[][], a matrix which saves the MCMC data.
     Num_Chain, the total number of chains.
     Num_Parameter, the total number of parameters.
     Num_Pair, the number of pairs of crossover chains.
     Indx_Chain, the index of current chain
     Indx_Generation, the index of current generation.
     Num_Cr, the total number of CR values.
     Cr[], a vector which saves the CR values.
     Cr_Prob[], the probability of each individual CR values.
     Parameter_Bounds[][2], the bounds for each patameter.
     Indx_Bounds: the index of bounds condition. 0, the fold bounds; 1, the reflect bounds; 2, set to the bound value
     idum, the random seed.
     Output parameters:
     Indx_Cr, the CR index.
     Dream_Sample[], a vector which saves the candidat parameters
     ***********************************************************************************/
    
    double Jumprate, Small_c=0.1;//, Small_cs=1e-14;
    
    int i, Num_Jump, Indx_Parameter;
    
    int Indx_Cr;
    
    Indx_Cr=Sample_Cr(Num_Cr, Cr_Prob, idum);
        
    int *Jump_dim;
    
    Jump_dim=VECTOR_INT(0, Num_Parameter-1);
    
    double *Diff;
    
    Diff=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    Dream_Dim(Num_Parameter, Cr, Indx_Cr, idum, &Num_Jump, Jump_dim);
    
    Dream_Diff(RESULT_DREAM, Indx_Generation, Num_Chain, Num_Parameter, Num_Pair, Num_Jump, Jump_dim, idum, Diff);
    
    double *Factor_E, *Factor_Epsilon;
    
    Factor_E=VECTOR_DOUBLE(0, Num_Jump-1);
    
    Factor_Epsilon=VECTOR_DOUBLE(0, Num_Jump-1);
    
    for (i=0; i<Num_Jump; i++) {
        
        Factor_E[i]=(Ran1(idum)-0.5)*Small_c;
        
        Factor_Epsilon[i]=0;
        //Factor_Epsilon[i]=GASDEV(idum)*Small_cs;
    }
    
    for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++) {
        
        Dream_Sample[Indx_Parameter] =RESULT_DREAM[Indx_Chain][Indx_Parameter+Num_Parameter*(Indx_Generation-1)];
        
    }
    
    /*
    if (Indx_Generation%10==0){
        Jumprate=1;
    }else{
        Jumprate=2.38/sqrt(2.*Num_pair*Num_Jump);
    }
    */
    
    Jumprate=2.38/sqrt(2.*Num_Pair*Num_Jump);
    
    for (i=0; i<Num_Jump; i++) {
        
        Dream_Sample[Jump_dim[i]] += Diff[i]*(1.0+Factor_E[i])*Jumprate + Factor_Epsilon[i];
        
    }
    
    Bounds_Enforce(Dream_Sample, Num_Parameter, Parameter_Bounds, Indx_Bounds);
    
    FREE_VECTOR_INT(Jump_dim, 0);
    
    FREE_VECTOR_DOUBLE(Diff, 0);
    
    FREE_VECTOR_DOUBLE(Factor_E, 0);
    
    FREE_VECTOR_DOUBLE(Factor_Epsilon, 0);

    return Indx_Cr;
    
}


int Sample_PairNum(int Pair_Number, long *idum){
    
    /***********************************************************************************
     Purpose:
     Sample the number of chain pairs used to generate the jump.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     Pair_Number, the max number of chain pairs used to generate the jump.
     idum, the random seed
     Output parameters:
     Num_Pair, the number of chain pairs used to generate the jump.
     ***********************************************************************************/
    
    int i, Num_Pair=1;
    
    double tmp;
    
    tmp=Ran1(idum);
    
    for (i=1; i<=Pair_Number; i++) {
        
        if(tmp<1.0/Pair_Number){
            
            Num_Pair=i;
            
            break;
            
        }else{
            
            tmp-=1.0/Pair_Number;
            
        }
        
    }
    
    return Num_Pair;
    
}


void Ini_Cr(int Num_Cr, double *Cr, double *Cr_Prob, double *Cr_Delta, int *Cr_L){
    
    /***********************************************************************************
     Purpose:
     Initialize the CR (crossover probability) values.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     Num_Cr, the total number of CR values.
     Output parameters:
     Cr[], a vector which saves the CR values.
     Cr_Prob[], the probability of each individual CR values.
     Cr_Delta[], the squared normalized jumping distance.
     Cr_L[], the number of updates for each CR value.
     ***********************************************************************************/
    
    int i;
    
    for (i=1; i<=Num_Cr; i++) {
        
        Cr[i]=(double)(i)/Num_Cr;
        
        Cr_Delta[i]=1;
        
        Cr_Prob[i]=1.0/Num_Cr;
        
        Cr_L[i]=0;
        
    }
    
    return;
    
}


int Sample_Cr(int Num_Cr, double *Cr_Prob, long *idum){
    
    /***********************************************************************************
     Purpose:
     Sample a CR index from the total number of Cr values according to the probability of each individual CR values.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     Num_Cr, the total number of CR values.
     Cr_Prob[], the probability of each individual CR values.
     idum, the random seed
     Output parameters:
     Sample_CR, the CR index.
     ***********************************************************************************/
    
    int i;
    
    double tmp;
    
    tmp=Ran1(idum);
    
    for (i=1; i<=Num_Cr; i++) {
        
        tmp -= Cr_Prob[i];
        
        if (tmp<0) {
            
            return i;
            
        }
        
    }
    
    return Num_Cr;
    
}


void Cr_dis_Update(double **RESULT_DREAM, int Num_Chain, int Num_Parameter, int Indx_Chain, int Indx_Generation, int Indx_Cr, int *Cr_L, double *Cr_Delta){
    
    /***********************************************************************************
     Purpose:
     Compute the squared normalized jumping distance.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     RESULT_DREAM[][], a matrix saves MCMC data.
     Num_Chain, the total number of chains
     Num_Parameter, the total number of parameters.
     Indx_Chain, the index of current chain.
     Indx_Generation, the index of current generation.
     Indx_Cr, the index of CR value choosed.
     Cr_L[], the number of updates for each CR value.
     Cr_Delta[], the squared normalized jumping distance.
     Output parameters:
     Cr_L[], update the the number of updates for each CR value.
     Cr_Delta[], update the squared normalized jumping distance.
     ***********************************************************************************/
    
    int i;

    double *STD, *MEAN;
    
    STD=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    MEAN=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    MCMC_STD(Indx_Generation, Indx_Generation+1, Num_Chain, Num_Parameter, RESULT_DREAM, STD, MEAN);
    
    for (i=0; i<Num_Parameter; i++) {
        
        Cr_Delta[Indx_Cr] += (RESULT_DREAM[Indx_Chain][i+Num_Parameter*(Indx_Generation)]-RESULT_DREAM[Indx_Chain][i+Num_Parameter*(Indx_Generation-1)])*(RESULT_DREAM[Indx_Chain][i+Num_Parameter*(Indx_Generation)]-RESULT_DREAM[Indx_Chain][i+Num_Parameter*(Indx_Generation-1)])/STD[i]/STD[i];
        
    }
    
    Cr_L[Indx_Cr]++;
    
    FREE_VECTOR_DOUBLE(STD, 0);
    
    FREE_VECTOR_DOUBLE(MEAN, 0);

    return;
    
}


void Cr_Pro_Update(int Num_Cr, double *Cr_Delta, int *Cr_L, double *Cr_Prob){
    
    /***********************************************************************************
     Purpose:
     Update the probability of each individual CR values.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     Num_Cr, the total number of CR values.
     Cr_Delta[], the squared normalized jumping distance.
     Cr_L[], the number of updates for each CR value.
     Output parameters:
     Cr_Prob[], the probability of each individual CR values.
     ***********************************************************************************/
    
    int i;
    
    double Tot_Prob=0, Tot_Dis=0;
    
    for (i=1; i<=Num_Cr; i++) {
        
        Tot_Dis += Cr_Delta[i];
        
    }
    
    for (i=1; i<=Num_Cr; i++) {
        
        Cr_Prob[i]=Cr_Delta[i]/Cr_L[i]/Tot_Dis;
        
    }
    
    for (i=1; i<=Num_Cr; i++) {
        
        Tot_Prob += Cr_Prob[i];
        
    }
    
    for (i=1; i<=Num_Cr; i++) {
        
        Cr_Prob[i] /= Tot_Prob;
        
    }
    
    return;
    
}


void Dream_Diff(double **RESULT_DREAM, int Indx_Generation, int Num_Chain, int Num_Parameter, int Num_Pair, int Num_Jump, int *Jump_dim, long *idum, double *Diff){
    
    /***********************************************************************************
     Purpose:
     Calculate the pairs differences used to sample new candidates.
     Modified:
     27 Mar 2018, Hao Li
     Input parameters:
     RESULT_DREAM[][], a matrix saves MCMC data.
     Indx_Generation, the index of current generation.
     Num_Chain, the total number of chains
     Num_Parameter, the total number of parameters.
     Num_Pair, the number of pairs of crossover chains.
     Num_Jump, the number of dimensions in which a jump will be made.
     Jump_dim, the dimensions in which a jump is to be made.
     idum, random seed.
     Output parameters:
     Diff[], a vector which saves pairs differences.
     ***********************************************************************************/
    
    int i, j;
    
    int *Pairs;
    
    Pairs=VECTOR_INT(0, 2*Num_Pair-1);
    
    int tmp1, tmp2;
    
    for (i=0; i<Num_Pair; i++) {
        
        do {
            
            tmp1=(int)(Ran1(idum)*Num_Chain);
            
            tmp2=(int)(Ran1(idum)*Num_Chain);
            
        } while (tmp1==tmp2||tmp1==Num_Chain||tmp2==Num_Chain);
        
        Pairs[i*2]=tmp1;
        
        Pairs[i*2+1]=tmp2;
        
    }
    
    for (i=0; i<Num_Jump; i++) {
        
        Diff[i] = 0;
        
        for (j=0; j<Num_Pair; j++) {
            
            Diff[i] += RESULT_DREAM[Pairs[j*2]][Jump_dim[i]+(Indx_Generation-1)*Num_Parameter]-RESULT_DREAM[Pairs[j*2+1]][Jump_dim[i]+(Indx_Generation-1)*Num_Parameter];
            
        }
        
    }
    
    FREE_VECTOR_INT(Pairs, 0);

}


void Dream_Dim(int Num_Parameter, double *Cr, int Indx_Cr, long *idum, int *Num_Jump, int *Jump_dim){
    
    /***********************************************************************************
     Purpose:
     Choose the dimensions in which a jump is to be made.
     Modified:
     26 Mar 2018, Hao Li
     Input parameters:
     Num_Parameter, the total number of parameters.
     Cr[], a vector which saves the Cr values.
     Indx_Cr, index of the choosed Cr.
     idum, random seed.
     Output parameters:
     Num_Jump, the total number of jumping dimention.
     Jump_dim, a vector which saves the indexs of jumping parameters.
     ***********************************************************************************/
    
    int i=0;
    
    *Num_Jump=0;
    
    for (i=0; i<Num_Parameter; i++) {
        
        Jump_dim[i]=Num_Parameter;
        
    }
    
    for (i=0; i<Num_Parameter; i++) {
        
        if (Ran1(idum)<Cr[Indx_Cr]) {
            
            Jump_dim[*Num_Jump] = i;
            
            *Num_Jump += 1;
            
        }
        
    }
    
    ////if jump is empty, one randum parameter will sampled to avoid jump vector to have zero length //
    
    if (*Num_Jump==0) {
        
        Jump_dim[0]=(int)(Ran1(idum)*Num_Parameter);
        
        *Num_Jump = 1;
        
    }
    
    return;
    
}


void Bounds_Enforce(double *Parameter_Sample, int Num_Parameter, double **Parameter_Bounds, int Indx_Bounds){
    
    /***********************************************************************************
     Purpose:
     Checks the bounds of the parameters.
     Modified:
     26 Mar 2018, Hao Li
     Input parameters:
     Parameter_Sample[], a vector which saves the sampled parameters.
     Num_Parameter, the total number of parameters.
     Parameter_Bounds[][2], the bounds for each patameter.
     Indx_Bounds: the index of bounds condition. 0, the fold bounds; 1, the reflect bounds; 2, set to the bound value
     Output parameters:
     Parameter_Sample[], return the checked parameters.
     ***********************************************************************************/
    
    int i;
    
    double tmp;
    
    switch (Indx_Bounds) {
            
        case 0://fold bounds
            
            for (i=0; i<Num_Parameter; i++) {
                
                tmp=Parameter_Bounds[i][1]-Parameter_Bounds[i][0];
                
                if (tmp==0.) {
                    
                    Parameter_Sample[i]=Parameter_Bounds[i][0];
                    
                } else {
                    
                    while (Parameter_Sample[i]<Parameter_Bounds[i][0])
                        
                        Parameter_Sample[i]+=tmp;
                    
                    while (Parameter_Sample[i]>Parameter_Bounds[i][1])
                        
                        Parameter_Sample[i]-=tmp;
                    
                }
                
            }
            
            break;
            
        case 1:///reflect bounds
            
            for (i=0; i<Num_Parameter; i++) {
                
                tmp=Parameter_Bounds[i][1]-Parameter_Bounds[i][0];
                
                if (tmp==0.) {
                    
                    Parameter_Sample[i]=Parameter_Bounds[i][0];
                    
                } else {
                    
                    while (Parameter_Sample[i]<Parameter_Bounds[i][0]||Parameter_Sample[i]>Parameter_Bounds[i][1]) {
                        
                        if (Parameter_Sample[i]<Parameter_Bounds[i][0])
                            
                            Parameter_Sample[i] = 2*Parameter_Bounds[i][0]-Parameter_Sample[i];
                        
                        if (Parameter_Sample[i]>Parameter_Bounds[i][1])
                            
                            Parameter_Sample[i] = 2*Parameter_Bounds[i][1]-Parameter_Sample[i];
                        
                    }
                    
                }
                
            }
            
            break;
            
        case 2:///set to bound value
            
            for (i=0; i<Num_Parameter; i++) {
                
                tmp=Parameter_Bounds[i][1]-Parameter_Bounds[i][0];
                
                if (tmp==0.) {
                    
                    Parameter_Sample[i]=Parameter_Bounds[i][0];
                    
                } else {
                    
                    if (Parameter_Sample[i]<Parameter_Bounds[i][0])
                        
                        Parameter_Sample[i] = Parameter_Bounds[i][0];
                    
                    if (Parameter_Sample[i]>Parameter_Bounds[i][1])
                        
                        Parameter_Sample[i] = Parameter_Bounds[i][1];
                    
                }
                
            }
            
            break;
            
        default:
            
            break;
    }
    
    return;
    
}


void PSRF(double **RESULT_MCMC, int Indx_Generation, int Num_Chain, int Num_Parameter, double *Statistics_R){
    
    /***********************************************************************************
     Purpose:
     Computes the potential scale reducton factor (PSRF) used to check multichain MCMC convergence (omit the correction factor due to D tends to be large at convergence).
     Modified:
     12 Mar 2018, Hao Li
     Input parameters:
     RESULT_MCMC[][],a matrix saves MCMC data.
     Indx_Generation, the Index of current Generation.
     Num_Chain, the total number of chains.
     Num_Parameter, the total number of parameters.
     Output parameters:
     Statistics_R[Num_Parameter], a vector return PSRF for each parameter.
     Method:
     Gelman Rubin R statistics.
     References:
     Brooks, S. and Gelman, A. (1998). General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics, 7(4), 434-55.
     Gelman, A. and Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. Statistical Science, 7, 457-72.
     ***********************************************************************************/
    
    int begin=Indx_Generation/2;
    
    int Length = Indx_Generation-begin+1;
    
    int Parameter, Chain, Generation;
    
    double *AVG_Chain;
    
    AVG_Chain=(double *) malloc(Num_Chain*sizeof(double));
    
    double Average, Parameter_B, Parameter_W, Parameter_Sigma;
    
    for (Parameter=0; Parameter<Num_Parameter; Parameter++) {
        
        Average=0;
        
        for (Chain=0; Chain<Num_Chain; Chain++) {
            
            AVG_Chain[Chain]=0;
            
            for (Generation=begin; Generation<=Indx_Generation; Generation++)
                
                AVG_Chain[Chain] += RESULT_MCMC[Chain][Parameter+Generation*Num_Parameter];
        
            AVG_Chain[Chain] /= Length;
            
            Average += AVG_Chain[Chain];
            
        }
        
        Average /= Num_Chain;
        
        Parameter_B=0;
        
        for (Chain=0; Chain<Num_Chain; Chain++)
            
            Parameter_B+=(AVG_Chain[Chain]-Average)*(AVG_Chain[Chain]-Average);
        
        Parameter_B /= (Num_Chain-1.);
        
        Parameter_W=0;
        
        for (Generation=begin; Generation<=Indx_Generation; Generation++)
            
            for (Chain=0; Chain<Num_Chain; Chain++)
                
                Parameter_W += (RESULT_MCMC[Chain][Parameter+Generation*Num_Parameter]-AVG_Chain[Chain])*(RESULT_MCMC[Chain][Parameter+Generation*Num_Parameter]-AVG_Chain[Chain]);
        
        Parameter_W=Parameter_W/(Length-1.0)/Num_Chain;
        
        Parameter_Sigma=(Length-1.)/Length*Parameter_W+Parameter_B;
        
        Statistics_R[Parameter]=sqrt((Num_Chain+1.0)/Num_Chain*Parameter_Sigma/Parameter_W-(Length-1.0)/Length/Num_Chain);
        
    }
    
    return;
    
}


int Rm_Outlierchain(int Num_Chain, int Indx_Generation, int Num_parameter, double **RESULT_MCMC, double **likelihood){
    
    /***********************************************************************************
     Purpose:
     Sorts an array ra[1..n] into ascending numerical order using the Heapsort algorithm.
     Method:
     Modified:
     15 Jun 2018, Hao Li
     Input parameters:
     Num_Chain, the total number of MCMC chains.
     Indx_Generation, the index of current Generation.
     RESULT_MCMC[][], a matrix which saves MCMC data.
     likelihood[][], a matrix which saves likelihood of the MCMC simulation.
     Output parameters:
     Num_Outlierchain, the total number of removed chains.
     reference:
     Vrugt, J. A., et al. 2009. International Journal of Nonlinear Science & Numerical Simulation, 10(3), 273-290
     ***********************************************************************************/
    
    int begin=Indx_Generation/2;
    
    int Length = Indx_Generation-begin+1;
    
    double *AVG_likelihood, *Tmp_likelihood;
    
    AVG_likelihood=(double *)malloc(Num_Chain*sizeof(double));
    
    Tmp_likelihood=(double *)malloc(Num_Chain*sizeof(double));
    
    int i, j, k, BestIndx=0;
    
    double Bestlikelihood=0;
    
    for (j=begin; j<=Indx_Generation; j++)
        
        AVG_likelihood[0]+=likelihood[0][j];
    
    AVG_likelihood[0]/=Length;
    
    Tmp_likelihood[0]=AVG_likelihood[0];
    
    Bestlikelihood=AVG_likelihood[0];
    
    for (i=1; i<Num_Chain; i++) {
        
        AVG_likelihood[i]=0;
        
        for (j=begin; j<=Indx_Generation; j++)
            
            AVG_likelihood[i]+=likelihood[i][j];
        
        AVG_likelihood[i]/=Length;
        
        Tmp_likelihood[i]=AVG_likelihood[i];
        
        if (AVG_likelihood[i]>Bestlikelihood) {
            
            Bestlikelihood=AVG_likelihood[i];
            
            BestIndx=i;
            
        }
        
    }
    
    HPSORT(Num_Chain, Tmp_likelihood);
    
    double Qr, Omega;
    
    Qr=Tmp_likelihood[Num_Chain/4*3]-Tmp_likelihood[Num_Chain/4];
    
    Omega=Tmp_likelihood[Num_Chain/4]-Qr*2;
    
    int Num_Outlierchain=0;
    
    for (i=0; i<Num_Chain; i++) {
        
        if (AVG_likelihood[i]<Omega) {
            
            Num_Outlierchain++;
            
            AVG_likelihood[i]=AVG_likelihood[BestIndx];
            
            for (k=0; k<Num_parameter; k++)
                
                RESULT_MCMC[i][k+Num_parameter*Indx_Generation]=RESULT_MCMC[BestIndx][k+Num_parameter*Indx_Generation];
            
            for (j=0; j<=Indx_Generation; j++)
                
                likelihood[i][j]=likelihood[BestIndx][j];
            
        }
        
    }
    
    free(AVG_likelihood);
    
    free(Tmp_likelihood);
    
    return Num_Outlierchain;
    
}


void HPSORT(long n, double *ra){
    
    /***********************************************************************************
     Purpose:
     Sorts an array ra[1..n] into ascending numerical order using the Heapsort algorithm.
     Modified:
     Numerical recipes in C 2ed.
     Input parameters:
     n, the number of elements the array.
     ra[], the array input.
     Output parameters:
     ra, replaced by its sorted rearrangement
     reference:
     Numerical recipes in C 2ed.
     ***********************************************************************************/
    
    long i,ir,j,l;
    
    double rra;
    
    if (n < 2)
        return;
    
    l=(n >> 1)+1;
    
    ir=n;
    
    //    The index l will be decremented from its initial value down to 1 during the “hiring” (heap creation) phase. Once it reaches 1, the index ir will be decremented from its initial value down to 1 during the “retirement-and-promotion” (heap selection) phase.
    for (;;) {
        
        if(l>1){
            
            rra=ra[--l];
            
        }else{
            
            rra=ra[ir];
            
            ra[ir]=ra[1];
            
            if (--ir == 1) {
                
                ra[1]=rra;
                
                break;
                
            }
            
        }
        
        i=l;
        
        j=l+l;
        
        while (j <= ir) {
            
            if (j < ir && ra[j] < ra[j+1])
                
                j++;
            
            if (rra < ra[j]) {
                
                ra[i]=ra[j];
                
                i=j;
                
                j <<= 1;
                
            }
            
            else break;
            
        }
        
        ra[i]=rra;
        
    }
    
    return;
    
}


void MCMC_STD(int Begin_Generation, int End_Generation, int Num_Chain, int Num_Parameter, double **RESULT_MCMC, double *STD, double *Mean){
    
    /***********************************************************************************
     Purpose:
     Calculate the standard deviations for each parameters
     Modified:
     15 Jun 2018, Hao Li
     Input parameters:
     Begin_Generation, the fist generation (include) used in the calculation.
     End_Generation, the last generation (not include) used in the calculation.
     Num_Chain, the total number of chains.
     Num_Parameter, the total number of parameters.
     RESULT_MCMC[][], a matrix which saves MCMC data.
     Output parameters:
     STD[], return the standard deviation of each parameter.
     Mean[], return the mean value of each patameter.
     Method:
     Sigma=sqrt((sum(Xi^2)-N*Xbar*Xbar)/(N-1)).
     ***********************************************************************************/
    
    int i, j, k;
    
    double AVG, SUM;
    
    for (i=0; i<Num_Parameter; i++) {
        
        AVG=0.;
        
        SUM=0.;
        
        for (j=Begin_Generation; j<End_Generation; j++)
            
            for (k=0; k<Num_Chain; k++)
                
                AVG += RESULT_MCMC[k][i+j*Num_Parameter];
        
        Mean[i]=AVG/(End_Generation-Begin_Generation)/Num_Chain;
        
        for (j=Begin_Generation; j<End_Generation; j++)
            
            for (k=0; k<Num_Chain; k++)
                
                SUM += (RESULT_MCMC[k][i+j*Num_Parameter]-Mean[i])*(RESULT_MCMC[k][i+j*Num_Parameter]-Mean[i]);
        
        STD[i]=sqrt(SUM/((End_Generation-Begin_Generation)*Num_Chain-1));
        
    }
    
    return;
    
}


void Chain_Ini(int Num_Chain, int Num_Wav, double *Wavelength_Data, double **Profile_Data, int Num_Parameter, double **Parameter_Bounds, TRANSITION_INFORMATION ATOM, double **RESULT_MCMC,  double *likelihood){
    
    /***********************************************************************************
     Purpose:
     Initialize the chains, and generate the first generation.
     Modified:
     21 Mar 2018, Hao Li
     Input parameters:
     Num_Chain_Gemc, the total number of GEMC chains.
     Num_profile, the total wavelength number.
     Profile_Data[Num_profile][5], the Stokes profiles used to inverse the magnetic field.
     Num_Parameter, the total number of parameters. (now only 9 is accepted)
     Parameter_Bounds[Num_Parameter][2], the bound for each patameter.
     ATOM, a sturcture which saves Atomic information.
     Output parameters:
     RESULT_MCMC[][], return the first generations of the MCMC chains.
     Likelihood, return the log-likelihood of each chain in the first generation.
     ***********************************************************************************/
    
    int i, j;
    
    long *idum;
    
    idum=(long *)malloc(sizeof(long));
    
    Random_Seed(idum);
    
    for (i=0; i<Num_Chain; i++) {
        
        for (j=0; j<Num_Parameter; j++)
           
            RESULT_MCMC[i][j]=Ran1(idum)*(Parameter_Bounds[j][1]-Parameter_Bounds[j][0])+Parameter_Bounds[j][0];
        
        likelihood[i]=Likelihood_Log(RESULT_MCMC[i], Num_Wav, Wavelength_Data, Profile_Data, ATOM);
    }
    
    free(idum);
    
    return;
    
}


int GEMC(int Num_Chain_Gemc, int Num_Wav, double *Wavelength_Data, double **Profile_Data, int Num_Parameter, double **Parameter_Bounds, TRANSITION_INFORMATION ATOM, double **GEMC_RESULT, double *LIKELIHOOD, int Indx_Bounds){
    
    /***********************************************************************************
     Purpose:
     GEMC (Genetic Evolution Markov Chain) algorithm, used to inverse the magnetic field.
     Modified:
     21 Mar 2018,  Hao Li
     Input parameters:
     Num_Chain_Gemc, the total number of of GEMC chains.
     Num_profile, the total wavelength number.
     Profile_Data[Num_profile][5], the Stokes profiles used to inverse the magnetic field.
     Num_Parameter, the total number of parameters. (now only 9 is accepted)
     Parameter_Bounds[Num_Parameter][2], the bound for each patameter.
     LIKELIHOOD_THRESHOLD, the threshold to stop the GEMC iteritaion.
     ATOM, a sturcture which saves Atomic information.
     Indx_Bounds, Index of bounds condition.
     Output parameters:
     GEMC_RESULT, return the last generation of GEMC simulation
     LIKELIHOOD, return the likelihood of the last generation
     Chain_bestfit, return the index of best fitting chain in the last generation.
     reference:
     Tregloan-Reed, J., et al. 2012 MNRAS
     ***********************************************************************************/
    
    long *idum;
    
    idum=(long *)malloc(sizeof(long));
    
    Random_Seed(idum);
    
    double *GEMC_SAMPLE;
    
    GEMC_SAMPLE=VECTOR_DOUBLE(0, Num_Parameter-1);
    
    int Indx_Chain=0, Indx_Parameter, Chain_bestfit=0, count_sample=0, count_accept=0, ITERATION=0, REIT=0;
    
    double TMP_LIKELIHOOD, Varance=0, AVG_LIKELIHOOD=0, SUM_SQUARE=0;
    
    for (Indx_Chain=1; Indx_Chain<Num_Chain_Gemc; Indx_Chain++)
        
        if (LIKELIHOOD[Indx_Chain]>LIKELIHOOD[Chain_bestfit])
            
            Chain_bestfit=Indx_Chain;
    
    double Delta=0, Delta_chi=0, Delta_variance=0, chi=0;
    
    for (REIT=0; ; REIT++) {
       
        for (ITERATION=0; ; ITERATION++) {
      
            for (Indx_Chain=0; Indx_Chain<Num_Chain_Gemc; Indx_Chain++) {
                
                if (Indx_Chain!=Chain_bestfit) {
                    
                    for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++) {
                        
                        Delta=GEMC_RESULT[Chain_bestfit][Indx_Parameter] - GEMC_RESULT[Indx_Chain][Indx_Parameter];
                        
                        if (fabs(Delta/GEMC_RESULT[Chain_bestfit][Indx_Parameter])>2e-6)
                            
                            GEMC_SAMPLE[Indx_Parameter]=GEMC_RESULT[Indx_Chain][Indx_Parameter] + Delta*Ran1(idum)*2.;
                        
                        else
                            
                            GEMC_SAMPLE[Indx_Parameter]=GEMC_RESULT[Chain_bestfit][Indx_Parameter] + GEMC_RESULT[Chain_bestfit][Indx_Parameter]*(Ran1(idum)-0.5)*5e-6;
                        
                    }
                    
                    count_sample++;
                    
                    Bounds_Enforce(GEMC_SAMPLE, Num_Parameter, Parameter_Bounds, Indx_Bounds);
                    
                    TMP_LIKELIHOOD=Likelihood_Log(GEMC_SAMPLE, Num_Wav, Wavelength_Data, Profile_Data, ATOM);
                    
                    if (TMP_LIKELIHOOD>LIKELIHOOD[Indx_Chain]) {
                        
                        count_accept++;
                        
                        LIKELIHOOD[Indx_Chain]=TMP_LIKELIHOOD;
                        
                        for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++)
                            
                            GEMC_RESULT[Indx_Chain][Indx_Parameter]=GEMC_SAMPLE[Indx_Parameter];
                        
                    }else if (Ran1(idum)<exp(TMP_LIKELIHOOD-LIKELIHOOD[Indx_Chain])) {
                        
                        count_accept++;
                        
                        LIKELIHOOD[Indx_Chain]=TMP_LIKELIHOOD;
                        
                        for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++)
                            
                            GEMC_RESULT[Indx_Chain][Indx_Parameter]=GEMC_SAMPLE[Indx_Parameter];
                        
                    }
                    
                } else {
                    
                    count_sample++;

                    for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++)
                        
                        GEMC_SAMPLE[Indx_Parameter]=GEMC_RESULT[Chain_bestfit][Indx_Parameter] + GEMC_RESULT[Chain_bestfit][Indx_Parameter]*(Ran1(idum)-0.5)*2e-6;
                    
                    Bounds_Enforce(GEMC_SAMPLE, Num_Parameter, Parameter_Bounds, Indx_Bounds);
                    
                    TMP_LIKELIHOOD=Likelihood_Log(GEMC_SAMPLE, Num_Wav, Wavelength_Data, Profile_Data, ATOM);

                    if (TMP_LIKELIHOOD>LIKELIHOOD[Indx_Chain]) {
                        
                        count_accept++;
                        
                        LIKELIHOOD[Indx_Chain]=TMP_LIKELIHOOD;
                        
                        for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++)
                            
                            GEMC_RESULT[Indx_Chain][Indx_Parameter]=GEMC_SAMPLE[Indx_Parameter];
                        
                    }else if (Ran1(idum)<exp(TMP_LIKELIHOOD-LIKELIHOOD[Indx_Chain])) {
                        
                        count_accept++;
                        
                        LIKELIHOOD[Indx_Chain]=TMP_LIKELIHOOD;
                        
                        for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++)
                            
                            GEMC_RESULT[Indx_Chain][Indx_Parameter]=GEMC_SAMPLE[Indx_Parameter];
                        
                    }

                }
                
            }
            
            Chain_bestfit=0;
            
            for (Indx_Chain=1; Indx_Chain<Num_Chain_Gemc; Indx_Chain++)
                
                if (LIKELIHOOD[Indx_Chain]>LIKELIHOOD[Chain_bestfit])
                    
                    Chain_bestfit=Indx_Chain;
            
            AVG_LIKELIHOOD=0;
            
            SUM_SQUARE=0;
            
            for (Indx_Chain=0; Indx_Chain<Num_Chain_Gemc; Indx_Chain++)
                
                AVG_LIKELIHOOD+=LIKELIHOOD[Indx_Chain];
            
            AVG_LIKELIHOOD/=Num_Chain_Gemc;
            
            for (Indx_Chain=0; Indx_Chain<Num_Chain_Gemc; Indx_Chain++)
                
                SUM_SQUARE+=(LIKELIHOOD[Indx_Chain]-AVG_LIKELIHOOD)*(LIKELIHOOD[Indx_Chain]-AVG_LIKELIHOOD);
            
            Delta_variance=Varance;
            
            Varance=sqrt(SUM_SQUARE/(Num_Chain_Gemc-1));
            
            Delta_variance -= Varance;
            
            if (ITERATION>=50||fabs(Delta_variance/Varance)<1e-4)
                
                break;
            
        }
        
        Delta_chi=LIKELIHOOD[Chain_bestfit]-chi;
        
        chi=LIKELIHOOD[Chain_bestfit];
        
        if (REIT>=15||fabs(Delta_chi/chi)<2e-4) break;

        for (Indx_Chain=0; Indx_Chain<Num_Chain_Gemc; Indx_Chain++) {
            
            if (Indx_Chain!=Chain_bestfit) {
                
                for (Indx_Parameter=0; Indx_Parameter<Num_Parameter; Indx_Parameter++)
                    
                    GEMC_RESULT[Indx_Chain][Indx_Parameter]=Ran1(idum)*(Parameter_Bounds[Indx_Parameter][1]-Parameter_Bounds[Indx_Parameter][0])+Parameter_Bounds[Indx_Parameter][0];
                
                LIKELIHOOD[Indx_Chain]=Likelihood_Log(GEMC_RESULT[Indx_Chain], Num_Wav, Wavelength_Data, Profile_Data, ATOM);
            
            }
            
        }
        
    }
    
    printf("GEMC sample = %d Best log likelihood = %e\n",count_sample, LIKELIHOOD[Chain_bestfit]);
    
    free(idum);
    
    FREE_VECTOR_DOUBLE(GEMC_SAMPLE, 0);
    
    return Chain_bestfit;
    
}
