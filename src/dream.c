
#include "dream.h"
#include "chain_statistics.h"
#include "likelihood.h"
#include "logger.h"
#include "profile_io.h"
#include "random.h"
#include "sample_analysis.h"
#include "sorting.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

#define DREAM_VERBOSE(mpi, level, fmt, ...)                                   \
    do{                                                                       \
      if((mpi)->island_rank == 0 && (mpi)->verbose_level >= (level)){         \
        snprintf(message_buffer, sizeof(message_buffer), fmt, ##__VA_ARGS__); \
        LOG_WRITE(message_buffer, true, true);                                \
      }                                                                       \
    }while(0)

#define rhat_limit 1.2

/*--------------------------------------------------------------------------------*/

static int Chain_Init_Dream(STRUCT_MPI *mpi, STRUCT_DREAM *dream,
    STRUCT_PARA *params, STRUCT_STK *stokes);
static int Dream_Sample(STRUCT_MPI *mpi, STRUCT_DREAM *dream,
    STRUCT_PARA *params, int npairs, double **chains, int chain_idx,
    int generation_idx, double *sample);
static int Sample_PairNum(int max_pairs, STRUCT_MPI *mpi);
static int Init_Cr(STRUCT_DREAM *dream);
static int Cr_Prob(STRUCT_DREAM *dream, STRUCT_MPI *mpi);
static int Cr_distance(STRUCT_DREAM *dream, int nchains, int chain_idx,
    int generation_idx, int cr_idx);
static int Sample_Cr(STRUCT_DREAM *dream, STRUCT_MPI *mpi);
static int Dream_Dim(STRUCT_MPI *mpi, STRUCT_DREAM *dream, int cr_idx);
static int Dream_Diff(STRUCT_MPI *mpi, double **chains, int chain_idx,
    int npairs, int njump, int *jump_dims, double *diff);
static int Rm_Outlierchain(STRUCT_MPI *mpi, STRUCT_DREAM *dream,
    int generation);
static int Bounds_Enforce(double *sample, STRUCT_PARA *params);

/*--------------------------------------------------------------------------------*/

static int DREAM_Allgather(double **chains, double *likelihood, STRUCT_MPI *mpi, 
    int nparams){

    /*######################################################################
      Purpose:
        Gather the current chain states and log-likelihoods within an island.
      Input parameters:
        chains, current parameter values for all chains.
        likelihood, current log-likelihood values for all chains.
        mpi, MPI island layout and local chain distribution.
        nparams, number of parameters in each chain.
      Output parameters:
        chains, populated with the states from every island rank.
        likelihood, populated with the values from every island rank.
      Return:
        0 on success.
     ######################################################################*/

    if(!chains || !chains[0] || !likelihood || !mpi || nparams < 1
        || mpi->chains_per_rank < 1 || mpi->island_comm == MPI_COMM_NULL
        || mpi->chains_per_rank > INT_MAX/nparams)
        return -1;
    if(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, chains[0],
        mpi->chains_per_rank*nparams, MPI_DOUBLE, mpi->island_comm)
        != MPI_SUCCESS) return -1;
    if(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, likelihood,
        mpi->chains_per_rank, MPI_DOUBLE, mpi->island_comm)
        != MPI_SUCCESS) return -1;

    return 0;
}

static int DREAM_Allgather_Likelihood(double *likelihood, STRUCT_MPI *mpi){

    /*######################################################################
      Purpose:
        Gather the current log-likelihood values within an island.
      Input parameters:
        likelihood, local and remote chain log-likelihood values.
        mpi, MPI island layout and local chain distribution.
      Output parameters:
        likelihood, populated with the values from every island rank.
      Return:
        0 on success.
     ######################################################################*/

    if(!likelihood || !mpi || mpi->chains_per_rank < 1
        || mpi->island_comm == MPI_COMM_NULL) return -1;
    if(MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, likelihood,
        mpi->chains_per_rank, MPI_DOUBLE, mpi->island_comm)
        != MPI_SUCCESS) return -1;

    return 0;
}

/*--------------------------------------------------------------------------------*/

static int DREAM_CENTER_PHI(double ***chains, int last_generation,
    STRUCT_MPI *mpi, STRUCT_DREAM *dream, STRUCT_PARA *params,
    double *best_sample){

    /*######################################################################
      Purpose:
        Center the spherical azimuth interval on the best burn-in solution
        and map all retained chains to the equivalent 180-degree branch.
      Input parameters:
        chains, circular DREAM history buffer.
        last_generation, final burn-in generation.
        mpi, MPI island and chain layout.
        dream, DREAM history dimensions.
        params, free-parameter mapping and current bounds.
        best_sample, best sample found during burn-in.
      Output parameters:
        chains, with spherical azimuths mapped around the best solution.
        params, with the azimuth bounds centered on the best solution.
        best_sample, with its azimuth mapped into the new interval.
      Return:
        0 on success or when no azimuth adjustment is required; -1 for
        invalid input.
      Note:
        Stokes profiles are invariant under PhiB -> PhiB+Pi. Burn-in samples
        are discarded, so selecting one equivalent branch before production
        sampling does not change the target distribution.
    ######################################################################*/

    if(!chains || !mpi || !dream || !params || !best_sample
        || last_generation < 0 || dream->history_size < 1
        || mpi->total_chains < 1) return -1;
    if(params->magnetic_mode != MAGNETIC_SPHERICAL) return 0;

    int phi_idx = -1;
    for(int param_idx=0; param_idx<params->npar; param_idx++){
      if(params->free_index[param_idx] == 2){
        phi_idx = param_idx;
        break;
      }
    }
    if(phi_idx < 0) return 0;

    const double old_width = params->bounds[phi_idx][1]
        -params->bounds[phi_idx][0];
    if(!isfinite(old_width) || old_width < M_PI*(1.0-1e-12)) return 0;

    const double center = best_sample[phi_idx];
    if(!isfinite(center)) return -1;
    params->bounds[phi_idx][0] = center-0.5*M_PI;
    params->bounds[phi_idx][1] = center+0.5*M_PI;

    const int history_count = last_generation+1 < dream->history_size
        ? last_generation+1 : dream->history_size;
    const int begin_generation = last_generation+1-history_count;
    for(int generation=begin_generation; generation<=last_generation;
        generation++){
      const int slot = generation%dream->history_size;
      for(int chain_idx=0; chain_idx<mpi->total_chains; chain_idx++){
        const double phi = chains[slot][chain_idx][phi_idx];
        if(!isfinite(phi)) return -1;
        chains[slot][chain_idx][phi_idx] = center+remainder(phi-center, M_PI);
      }
    }
    best_sample[phi_idx] = center;

    DREAM_VERBOSE(mpi, 4,
        "PhiB bounds centered after burn-in: [%.6f, %.6f] Pi.",
        params->bounds[phi_idx][0]/M_PI,
        params->bounds[phi_idx][1]/M_PI);

    return 0;
}

/*--------------------------------------------------------------------------------*/

int INIT_DREAM(STRUCT_MPI *mpi, STRUCT_DREAM *dream, STRUCT_PARA *params, 
    STRUCT_STK *stokes){

    /*######################################################################
      Purpose:
        Initialize all per-profile state required by the DREAM sampler.
      Input parameters:
        input, DREAM controls, parameter limits, and current profile 
          buffer.
        mpi, initialized MPI island layout and random-number-generator 
          states.
        dream, zero-initialized or previously used DREAM state.
        params, initialized spectral-line data and model workspace.
        stokes, wavelength grid and current observed profile.
      Output parameters:
        dream, with the first chain generation and likelihoods 
          initialized.
        params, with DREAM bounds, proposal scales, and policies 
          initialized.
        stokes, with the current profile, initial guess, and noise 
          initialized.
      Return:
        0 on success; 1 for a rejected low-intensity profile; -1 for invalid
          configuration or allocation failure.
     ######################################################################*/

    if(!mpi || !dream || !params || !stokes) return -1;
    const int nmodel = params->nmodel;
    if(nmodel < 1 || nmodel > MAX_MODEL_PARAMS) return -1;
    if(mpi->is_master || !mpi->rank_rng || mpi->total_chains < 3) return -1;
    if(dream->burnin_generations < 2 || dream->sampling_generations < 2
        || dream->ncr < 1 || dream->max_pairs < 1) return -1;
    if(stokes->nw < 1 || !stokes->wavelength || !stokes->synthetic 
        || !params->lines) return -1;

    int nfree = 0;
    for(int imodel=0; imodel<nmodel; imodel++){
      if(params->inv[imodel]) nfree++;
    }
    if(nfree < 1) return -1;

    if(!stokes->profile) return -1;

    if(!params->bounds){
      params->bounds = (double **)malloc(
          (size_t)nfree*sizeof(*params->bounds));
      if(!params->bounds) return -1;
      params->bounds[0] = (double *)malloc(
          (size_t)nfree*2*sizeof(**params->bounds));
      if(!params->bounds[0]){
        free(params->bounds);
        params->bounds = NULL;
        return -1;
      }
      for(int ipar=1; ipar<nfree; ipar++){
        params->bounds[ipar] = params->bounds[0]+ipar*2;
      }
    }

    if(!params->proposal_scale){
      params->proposal_scale = (double *)malloc(
          (size_t)nfree*sizeof(*params->proposal_scale));
    }
    if(!params->type){
      params->type = (enum bounds_type *)malloc(
          (size_t)nfree*sizeof(*params->type));
    }
    if(!params->narrower_guess){
      params->narrower_guess = (bool *)malloc(
          (size_t)nfree*sizeof(*params->narrower_guess));
    }
    if(!params->step){
      params->step = (double *)malloc(
          (size_t)nmodel*sizeof(*params->step));
    }
    if(!params->free_index){
      params->free_index = (int *)malloc(
          (size_t)nfree*sizeof(*params->free_index));
    }
    if(!stokes->inv_noise){
      stokes->inv_noise = (double *)malloc(
          (size_t)stokes->nw*4*sizeof(*stokes->inv_noise));
    }
    if(!params->proposal_scale || !params->type 
        || !params->narrower_guess || !params->step || !params->free_index 
        || !stokes->inv_noise) return -1;

    params->npar = nfree;
    dream->nparams = nfree;
    params->bx_free_index = -1;
    params->by_free_index = -1;
    int ifree = 0;
    for(int imodel=0; imodel<nmodel; imodel++){
      if(params->limits[imodel][0] > params->limits[imodel][1]) return -1;
      if(params->kind[imodel] != PARAM_CONTINUUM && !params->inv[imodel] 
          && (params->value_const[imodel] 
          < params->limits[imodel][0] || params->value_const[imodel] 
          > params->limits[imodel][1])) return -1;

      params->step[imodel] = 1e-3*(params->limits[imodel][1] 
          -params->limits[imodel][0]);
      if(params->inv[imodel]){
        params->free_index[ifree] = imodel;
        if(imodel == 1) params->bx_free_index = ifree;
        if(imodel == 2) params->by_free_index = ifree;
        params->type[ifree] = params->magnetic_mode == MAGNETIC_SPHERICAL 
            && params->kind[imodel] == PARAM_B2 ? enum_fold : enum_reflect;
        params->narrower_guess[ifree] = 
            params->kind[imodel] == PARAM_CONTINUUM;
        ifree++;
      }
    }

    int status = Profile_Prepare(stokes, params);
    if(status != 0) return status;

    if(Noise_Init(stokes) != 0) return -1;
    for(int ipar=0; ipar<nfree; ipar++){
      int imodel = params->free_index[ipar];
      const double lower = params->limits[imodel][0];
      const double upper = params->limits[imodel][1];
      const double width = upper-lower;
      if(!isfinite(lower) || !isfinite(upper) || !isfinite(width)
          || width < 0.0
          || (params->type[ipar] == enum_reflect
            && !isfinite(2.0*width))) return -1;
      params->bounds[ipar][0] = lower;
      params->bounds[ipar][1] = upper;
      params->proposal_scale[ipar] = 1e-3*width;
    }

    /* Reuse the structure safely when initializing another profile. */
    FREE_TENSOR(dream->chains);
    FREE_MATRIX(dream->likelihood);

    return Chain_Init_Dream(mpi, dream, params, stokes);
}

/*--------------------------------------------------------------------------------*/

int DREAM(STRUCT_MPI *mpi, STRUCT_DREAM *dream, STRUCT_PARA *params,
    STRUCT_STK *stokes, STRUCT_PROFILE_IO *input,
    const STRUCT_SUBSET *subset){

    /*######################################################################
      Purpose:
        Run the DREAM Markov-chain Monte Carlo sampler.
      Input parameters:
        mpi, MPI island layout and random-number-generator states.
        dream, DREAM state, history buffers, and crossover statistics.
        params, model parameters, bounds, and proposal scales.
        stokes, observed Stokes profiles used by the likelihood function.
        input, sample-file layout and parallel file handle.
        subset, coordinates of the profile being inverted.
      Output parameters:
        dream, updated samples, likelihoods, and crossover statistics.
        params, proposal scales updated after burn-in when enabled.
      Return:
        Best-chain index on success; -1 for invalid state, allocation, or MPI
        failure; -2 when a non-finite numerical result is detected.
     ######################################################################*/

    if(!mpi || !dream || !params || !stokes || !input || !subset
        || params->npar < 1
        || params->npar != dream->nparams || !params->free_index
        || !params->bounds || !params->proposal_scale
        || !dream->chains.data_ptr || !dream->likelihood.data_ptr
        || dream->history_size < 2 || dream->burnin_generations < 2
        || dream->sampling_generations < 2
        || mpi->total_chains < 3 || mpi->is_master
        || mpi->island_comm == MPI_COMM_NULL) return -1;

    double ***chains = TENSOR_DBL(dream->chains);
    double **likelihood = MAT_DBL(dream->likelihood);
    if(!chains || !likelihood) return -1;

    DREAM_VERBOSE(mpi, 2, "DREAM: burn-in up to %d generations, "
        "%d sampling generations, %d chains.", dream->burnin_generations,
        dream->sampling_generations, mpi->total_chains);

    bool Converged = false;
    int npairs, cr_idx;
    int count_sample=0, count_accept=0, tot_sample=0, tot_accept=0;
    int chain_idx=0, generation_idx, param_idx;
    int burnin_phase, best_chain = 0, last_generation = 0;
    double proposed_likelihood, best_likelihood;
    double sample[params->npar], rhat[params->npar];
    double best_sample[params->npar];
    double stddev[params->npar], mean[params->npar];
    int psrf_indices[4], npsrf = 0;
    int psrf_generation = -1;
    bool psrf_available = false;
    bool phi_centered = false;
    dream->retained_generations = 0;
    best_likelihood = -INFINITY;

    /* Only magnetic field parameters and line-of-sight velocity participate
       in the convergence diagnostic.  Fixed parameters have no chain index. */
    for(param_idx=0; param_idx<params->npar; param_idx++){
      if(params->free_index[param_idx] <= 3){
        psrf_indices[npsrf++] = param_idx;
      }
    }

    int init_status = Init_Cr(dream);
    if(mpi->island_size > 1){
      if(MPI_Allreduce(MPI_IN_PLACE, &init_status, 1, MPI_INT, MPI_MIN,
          mpi->island_comm) != MPI_SUCCESS) init_status = -1;
    }
    if(init_status != 0){
      Free_Dream(dream);
      return -1;
    }

    /* Chain initialization contains no collective operation, so the caller
       can first synchronize allocation failures across the island. */
    if(mpi->island_size > 1){
      if(DREAM_Allgather(chains[0], likelihood[0], mpi, params->npar) != 0){
        Free_Dream(dream);
        return -1;
      }
    }

    for(burnin_phase = 0; burnin_phase<2; burnin_phase++){
      int phase_generations = burnin_phase == 0
          ? dream->burnin_generations : dream->sampling_generations;
      count_sample=0;
      count_accept=0;
      last_generation = 0;
      for(generation_idx=1; generation_idx<phase_generations;
          generation_idx++){
        int curr_slot = generation_idx%dream->history_size;
        int prev_slot = (generation_idx-1)%dream->history_size;
        int numerical_error = 0;
        for(chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end; 
            chain_idx++){
         
          npairs = Sample_PairNum(dream->max_pairs, mpi);
          cr_idx = Dream_Sample(mpi, dream, params, npairs, 
              chains[prev_slot], chain_idx, generation_idx, sample);
          count_sample++;
          if(cr_idx < 1){
            numerical_error = -1;
            likelihood[curr_slot][chain_idx] =
                likelihood[prev_slot][chain_idx];
            for(param_idx=0; param_idx<params->npar; param_idx++){
              chains[curr_slot][chain_idx][param_idx] =
                  chains[prev_slot][chain_idx][param_idx];
            }
            continue;
          }
          int proposal_error = 0;
          proposed_likelihood = Likelihood_Log(sample, stokes, params,
              &proposal_error);
          if(proposal_error != 0) numerical_error = -1;

          
          if(log(RNG_UNIFORM(mpi->rank_rng))< 
              (proposed_likelihood-likelihood[prev_slot][chain_idx])){    
            count_accept++;
            likelihood[curr_slot][chain_idx] = proposed_likelihood;
            for(param_idx=0; param_idx<params->npar; param_idx++){
              chains[curr_slot][chain_idx][param_idx] = sample[param_idx];
            }
                  
          }else{
            likelihood[curr_slot][chain_idx] = 
                likelihood[prev_slot][chain_idx];
            for(param_idx=0; param_idx<params->npar; param_idx++){
              chains[curr_slot][chain_idx][param_idx] = 
                  chains[prev_slot][chain_idx][param_idx];
            }
          }

          if(dream->ncr>1 && burnin_phase == 0 && dream->update_crossover
              && generation_idx<1500){
            if(Cr_distance(dream, mpi->total_chains, chain_idx,
                generation_idx, cr_idx) != 0) numerical_error = -1;
          }
        }

        if(mpi->island_size > 1){
          if(MPI_Allreduce(MPI_IN_PLACE, &numerical_error, 1, MPI_INT,
              MPI_MIN, mpi->island_comm) != MPI_SUCCESS) numerical_error = -1;
        }
        if(numerical_error != 0){
          if(mpi->island_rank == 0){
            LOG_ERROR(ERR_LVL_ERROR, "DREAM",
                "Non-finite value detected in likelihood or chain statistics.\n");
          }
          Free_Dream(dream);
          return -2;
        }

        if(mpi->island_size > 1){
          if(DREAM_Allgather(chains[curr_slot], likelihood[curr_slot], mpi,
              params->npar) != 0){
            Free_Dream(dream);
            return -1;
          }
        }
        last_generation = generation_idx;

        if(burnin_phase == 1
            && (generation_idx+1)%dream->history_size == 0){
          int first_generation = generation_idx+1-dream->history_size;
          if(SAMPLE_WRITE_BLOCK(input, dream, params, mpi, subset,
              first_generation, dream->history_size) != 0){
            Free_Dream(dream);
            return -1;
          }
        }

        if(burnin_phase == 0){
          for(chain_idx=0; chain_idx<mpi->total_chains; chain_idx++){
            if(likelihood[curr_slot][chain_idx] > best_likelihood){
              best_likelihood = likelihood[curr_slot][chain_idx];
              best_chain = chain_idx;
              for(param_idx=0; param_idx<params->npar; param_idx++){
                best_sample[param_idx] = chains[curr_slot][chain_idx][param_idx];
              }
            }
          }
        }

        if(generation_idx%500 == 0){
          DREAM_VERBOSE(mpi, 2, "DREAM generation: %d", generation_idx);
        }

        if(generation_idx%30 == 0 && burnin_phase == 0){

          if(params->magnetic_mode == MAGNETIC_SPHERICAL && !phi_centered){
            if(DREAM_CENTER_PHI(chains, generation_idx, mpi, dream, params,
                best_sample) != 0){
              Free_Dream(dream);
              return -1;
            }
            phi_centered = true;
          }
        
          int psrf_status = PSRF(mpi, chains, generation_idx, psrf_indices, 
              npsrf, params->npar, dream->history_size, rhat);

          if(psrf_status == 0){
            psrf_available = true;
            psrf_generation = generation_idx;
          }
          if(psrf_status == 0 && !Converged){
            Converged = true;   
            for(param_idx=0; param_idx<npsrf; param_idx++){
              if(!(rhat[param_idx] < rhat_limit)){
                Converged = false;          
                break;
              }          
            }
          }

          if(Converged){
            DREAM_VERBOSE(mpi, 3,
                "Best burn-in chain: chain=%d, log likelihood=%e",
                best_chain, best_likelihood);
            if(mpi->island_rank == 0 && mpi->verbose_level >= 3){
              LOG_WRITE("Burn-in finished. Best-fit parameters:", true, true);
              for(param_idx = 0; param_idx<params->npar;param_idx++){
                snprintf(message_buffer, sizeof(message_buffer), 
                    "parameter %d: %e", param_idx, best_sample[param_idx]);
                LOG_WRITE(message_buffer, true, true);
              }
              snprintf(message_buffer, sizeof(message_buffer), 
                  "Log likelihood: %e", best_likelihood);
              LOG_WRITE(message_buffer, true, true);
            }
            break;
          }

          if(dream->ncr > 1 && dream->update_crossover &&generation_idx<1500
              && Cr_Prob(dream, mpi) != 0){
            Free_Dream(dream);
            return -1;
          }

          if(generation_idx < 2000 
              && Rm_Outlierchain(mpi, dream, generation_idx) < 0){
            Free_Dream(dream);
            return -1;
          }
	
        }
      }

      if(burnin_phase == 0){
        if(!isfinite(best_likelihood)){
          Free_Dream(dream);
          return -1;
        }

        if(psrf_available){
          DREAM_VERBOSE(mpi, 2, "Burn-in PSRF at generation %d:", 
              psrf_generation);
          for(int selected_idx=0; selected_idx<npsrf; selected_idx++){
            int free_idx = psrf_indices[selected_idx];
            int model_idx = params->free_index[free_idx];
            char name[64];
            MODEL_PARAMETER_NAME(params, model_idx, name, sizeof(name));
            DREAM_VERBOSE(mpi, 2, "  %s: R-hat = %.6f", name, 
                rhat[selected_idx]);
          }
        }else{
          DREAM_VERBOSE(mpi, 2, "%s", 
              "Burn-in PSRF unavailable: no successful PSRF evaluation " 
              "was completed.");
        }
        if(!Converged){
          DREAM_VERBOSE(mpi, 2, "Burn-in reached its maximum of %d "
              "generations without satisfying the PSRF threshold; "
              "continuing to sampling.", dream->burnin_generations);
        }
      }

      if(burnin_phase == 0){
        int final_slot = last_generation%dream->history_size;
        for(chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end; 
            chain_idx++){
          for(param_idx = 0; param_idx<params->npar;param_idx++){
            chains[0][chain_idx][param_idx] = 
                chains[final_slot][chain_idx][param_idx];
          }
          likelihood[0][chain_idx] = 
              likelihood[final_slot][chain_idx];
        }

        if(mpi->island_size > 1){
          if(DREAM_Allgather(chains[0], likelihood[0], mpi,
              params->npar) != 0){
            Free_Dream(dream);
            return -1;
          }
        }
      
      }else{
        dream->retained_generations = last_generation+1;
        int remaining_generations = dream->retained_generations
            %dream->history_size;
        if(remaining_generations > 0){
          int first_generation = dream->retained_generations
              -remaining_generations;
          if(SAMPLE_WRITE_BLOCK(input, dream, params, mpi, subset,
              first_generation, remaining_generations) != 0){
            Free_Dream(dream);
            return -1;
          }
        }
        int history_count = dream->retained_generations<dream->history_size 
            ? dream->retained_generations : dream->history_size;
        int begin_generation = dream->retained_generations-history_count;
        int best_generation = begin_generation;
        best_chain = 0;
        if(mpi->island_rank == 0){
          for(int retained_generation=begin_generation; 
              retained_generation<=last_generation; retained_generation++){
            int slot = retained_generation%dream->history_size;
            for(chain_idx=0; chain_idx<mpi->total_chains; chain_idx++){
              int best_slot = best_generation%dream->history_size;
              if(likelihood[slot][chain_idx] 
                  > likelihood[best_slot][best_chain]){
                best_generation = retained_generation;
                best_chain = chain_idx;
              }
            }
          }
        }
        if(mpi->island_size > 1){
          if(MPI_Bcast(&best_generation, 1, MPI_INT, 0, mpi->island_comm)
              != MPI_SUCCESS
              || MPI_Bcast(&best_chain, 1, MPI_INT, 0, mpi->island_comm)
              != MPI_SUCCESS){
            Free_Dream(dream);
            return -1;
          }
        }
        int best_slot = best_generation%dream->history_size;
        best_likelihood = likelihood[best_slot][best_chain];
        for(param_idx=0; param_idx<params->npar; param_idx++){
          best_sample[param_idx] = chains[best_slot][best_chain][param_idx];
        }
        DREAM_VERBOSE(mpi, 2, 
            "Sampling finished: generation=%d, best generation=%d, " 
            "best chain=%d, log likelihood=%e", last_generation, 
            best_generation, best_chain, best_likelihood);
        if(mpi->island_rank == 0 && mpi->verbose_level >= 2){
          double best_model[MAX_MODEL_PARAMS];
          for(int model_idx=0; model_idx<params->nmodel; model_idx++){
            best_model[model_idx] = params->value_const[model_idx];
          }
          for(param_idx=0; param_idx<params->npar; param_idx++){
            best_model[params->free_index[param_idx]] = best_sample[param_idx];
          }
          LOG_WRITE("Best-fit model parameters:", true, true);
          for(int model_idx=0; model_idx<params->nmodel; model_idx++){
            char name[64];
            MODEL_PARAMETER_NAME(params, model_idx, name, sizeof(name));
            DREAM_VERBOSE(mpi, 2, "  %s = %.10e%s", name, 
                best_model[model_idx], 
                params->inv[model_idx] ? "" : " (fixed)");
          }
        }
      }

      if(MPI_Allreduce(&count_sample, &tot_sample, 1, MPI_INT, MPI_SUM,
          mpi->island_comm) != MPI_SUCCESS
          || MPI_Allreduce(&count_accept, &tot_accept, 1, MPI_INT, MPI_SUM,
          mpi->island_comm) != MPI_SUCCESS){
        Free_Dream(dream);
        return -1;
      }

      DREAM_VERBOSE(mpi, 2, 
          "%s acceptance probability: %.2f%% (%d accepted / %d proposed).", 
          burnin_phase == 0 ? "Burn-in" : "Sampling", 
          tot_sample > 0 ? tot_accept*100.0/tot_sample : 0.0, 
          tot_accept, tot_sample);

      if(burnin_phase == 0 && params->magnetic_mode == MAGNETIC_SPHERICAL){
        if(DREAM_CENTER_PHI(chains, last_generation, mpi, dream, params,
            best_sample) != 0){
          Free_Dream(dream);
          return -1;
        }
      }

      if(dream->update_proposal_noise && burnin_phase == 0){
        int begin_stats = last_generation*3/4;
        if(last_generation-begin_stats+1 > dream->history_size){
          begin_stats = last_generation+1-dream->history_size;
        }
        if(Chains_STD(mpi, begin_stats, last_generation, params->npar,
            chains, dream->history_size, stddev, mean) != 0){
          Free_Dream(dream);
          return -1;
        }
        if(mpi->island_rank == 0 && mpi->verbose_level >= 4){
          for(param_idx = 0; param_idx<params->npar; param_idx++){
            snprintf(message_buffer, sizeof(message_buffer), 
                "parameter %d: best=%e, mean=%e, std=%e, normalized std=%e", 
                param_idx, best_sample[param_idx],
                mean[param_idx], stddev[param_idx], stddev[param_idx] 
                /(params->bounds[param_idx][1] 
                -params->bounds[param_idx][0]));
            LOG_WRITE(message_buffer, true, true);
          }
        }
        for(param_idx = 0; param_idx<params->npar; param_idx++){
          params->proposal_scale[param_idx] = stddev[param_idx];
        }
      }      
    }

    Free_Dream(dream);
    return best_chain;
}

/*--------------------------------------------------------------------------------*/

static int Chain_Init_Dream(STRUCT_MPI *mpi, STRUCT_DREAM *dream, STRUCT_PARA *params, 
    STRUCT_STK *stokes){
    
    /*######################################################################
      Purpose:
        Initialize DREAM chain states and their log-likelihoods.
      Input parameters:
        input, runtime controls including the requested generation count.
        mpi, MPI island layout and local chain range.
        params, model parameter bounds and initialization controls.
        stokes, observed Stokes profiles used by the likelihood function.
      Output parameters:
        dream, allocated history buffers containing the initial chain states.
      Return:
        0 on success.
     ######################################################################*/

    int param_idx, chain_idx;

    int max_generations = dream->burnin_generations
        > dream->sampling_generations ? dream->burnin_generations
        : dream->sampling_generations;
    dream->history_size = max_generations < DREAM_HISTORY_MAX
        ? max_generations : DREAM_HISTORY_MAX;
    if(dream->history_size < 2) dream->history_size = 2;

    dream->chains = TENSOR(0, dream->history_size-1, 0, 
        mpi->total_chains-1, 0, params->npar-1, enum_dbl, true);
    dream->likelihood = MATRIX(0, dream->history_size-1, 0, 
        mpi->total_chains-1, enum_dbl, false);

    if(!dream->chains.data_ptr || !dream->likelihood.data_ptr){
      FREE_TENSOR(dream->chains);
      FREE_MATRIX(dream->likelihood);
      return -1;
    }

    double ***chains = TENSOR_DBL(dream->chains);
    double **likelihood = MAT_DBL(dream->likelihood);

    for(chain_idx = mpi->chain_begin; chain_idx <= mpi->chain_end; 
        chain_idx++){
      for(param_idx = 0; param_idx < params->npar; param_idx++){
        if(params->narrower_guess[param_idx]){
          chains[0][chain_idx][param_idx] = 
              (RNG_UNIFORM(mpi->rank_rng)-0.5) 
              *(params->bounds[param_idx][1] 
              -params->bounds[param_idx][0])*0.3 
              +0.5*(params->bounds[param_idx][0] 
              +params->bounds[param_idx][1]);
        }else{
           chains[0][chain_idx][param_idx] = 
              RNG_UNIFORM(mpi->rank_rng) 
              *(params->bounds[param_idx][1] 
              -params->bounds[param_idx][0]) 
              +params->bounds[param_idx][0];
        }
      }
      int numerical_error = 0;
      likelihood[0][chain_idx] = Likelihood_Log(chains[0][chain_idx],
          stokes, params, &numerical_error);
      if(numerical_error != 0) return -2;
    }

    return 0;    
}



/*--------------------------------------------------------------------------------*/

static int Dream_Sample(STRUCT_MPI *mpi, STRUCT_DREAM *dream, STRUCT_PARA *params, 
    int npairs, double **chains, int chain_idx, int generation_idx, 
    double *sample){
    
    /*######################################################################
      Purpose:
        Generate one DREAM proposal for the selected chain.
      Input parameters:
        mpi, MPI state and random-number-generator states.
        dream, DREAM crossover and proposal work arrays.
        params, model parameter bounds and proposal scales.
        npairs, number of chain pairs used for the differential jump.
        chains, current states of all chains.
        chain_idx, index of the chain being updated.
        generation_idx, current generation index.
      Output parameters:
        sample, generated proposal after enforcing parameter bounds.
      Return:
        Index of the crossover value used for the proposal.
     ######################################################################*/
    
    if(!mpi || !dream || !params || !chains || !sample || npairs < 1
        || chain_idx < 0 || chain_idx >= mpi->total_chains
        || generation_idx < 1 || dream->nparams != params->npar
        || dream->nparams < 1) return -1;

    int cr_idx = Sample_Cr(dream, mpi);
    if(cr_idx < 1 || cr_idx > dream->ncr) return -1;

    int njump = Dream_Dim(mpi, dream, cr_idx);
    if(njump < 1 || njump > dream->nparams) return -1;

    if(Dream_Diff(mpi, chains, chain_idx, npairs, njump,
        dream->jump_dims, dream->diff) != 0) return -1;

    int i, param_idx;
    double gamma; 

    for(i=0; i<njump; i++){    
      dream->scale_noise[i]=(RNG_UNIFORM(mpi->rank_rng)-0.5)*0.1;
      dream->additive_noise[i]=RNG_GAUSS(mpi->rank_rng) 
          *params->proposal_scale[dream->jump_dims[i]]*1e-5;
    }
      
    for(param_idx=0; param_idx<params->npar; param_idx++){
      sample[param_idx] = chains[chain_idx][param_idx];
    }
    
    if(generation_idx%5==0){
      gamma = 1.0;
    }else{
      gamma = 2.38/sqrt(2.0*npairs*params->npar);
    }
    
    for(i=0; i<njump; i++){
      sample[dream->jump_dims[i]] += dream->diff[i]*(1.0+dream->scale_noise[i])*gamma 
        +dream->additive_noise[i]; 
    }

    if(Bounds_Enforce(sample, params) != 0) return -1;
    


    return cr_idx;
}

/*--------------------------------------------------------------------------------*/

static int Dream_Dim(STRUCT_MPI *mpi, STRUCT_DREAM *dream, int cr_idx){
    
    /*######################################################################
      Purpose:
        Select the parameter dimensions included in a DREAM jump.
      Input parameters:
        mpi, MPI state and random-number-generator states.
        dream, DREAM crossover values.
        params, model parameter metadata.
        cr_idx, index of the selected crossover value.
      Output parameters:
        jump_dims, indices of the selected parameter dimensions.
      Return:
        Number of selected dimensions.
     ######################################################################*/
    
    if(!mpi || !dream || !mpi->rank_rng || !dream->jump_dims
        || !dream->crossover || dream->nparams < 1 || cr_idx < 1
        || cr_idx > dream->ncr) return -1;

    int i=0, njump=0;
    
    for(i=0; i<dream->nparams; i++){
      dream->jump_dims[i] = -1;
    }
    
    for(i=0; i<dream->nparams; i++){
      if(RNG_UNIFORM(mpi->rank_rng)<=dream->crossover[cr_idx]){
        dream->jump_dims[njump] = i;
        njump++;
      }  
    }
    
    if(njump==0){
        dream->jump_dims[0]=(int)(RNG_UNIFORM(mpi->rank_rng) 
            *dream->nparams);
        njump = 1;
    }
    
    return njump;
    
}

/*--------------------------------------------------------------------------------*/

static int Dream_Diff(STRUCT_MPI *mpi, double **chains, int chain_idx, int npairs, 
    int njump, int *jump_dims, double *diff){
    
    /*######################################################################
      Purpose:
        Calculate differential-evolution jumps from random chain pairs.
      Input parameters:
        mpi, MPI state and random-number-generator states.
        chains, current states of all chains.
        chain_idx, index of the current chain, which cannot be selected.
        npairs, number of random chain pairs.
        njump, number of dimensions in the jump.
        jump_dims, indices of the dimensions in the jump.
      Output parameters:
        diff, summed pair differences for the selected dimensions.
      Return:
        0 on success.
     ######################################################################*/
        
    if(!mpi || !mpi->rank_rng || !chains || !jump_dims || !diff
        || mpi->total_chains < 3 || chain_idx < 0
        || chain_idx >= mpi->total_chains || npairs < 1
        || 2LL*npairs > (long long)mpi->total_chains-1 || njump < 1) return -1;
    for(int i=0; i<njump; i++){
      if(jump_dims[i] < 0) return -1;
    }

    int Pairs[2*npairs];
    int i, j;

    /* Select all differential chains without replacement. */
    for(i=0; i<2*npairs; i++){
      int candidate;
      bool duplicate;
      do{
        candidate = (int)(RNG_UNIFORM(mpi->rank_rng) 
            *mpi->total_chains);
        duplicate = candidate == chain_idx;
        for(j=0; j<i && !duplicate; j++){
          duplicate = candidate == Pairs[j];
        }
      }while(duplicate);
      Pairs[i] = candidate;
    }

    for(i=0; i<njump; i++){   
      diff[i] = 0;
      for(j=0; j<npairs; j++){
        diff[i] += chains[Pairs[j*2]][jump_dims[i]] 
            -chains[Pairs[j*2+1]][jump_dims[i]];
      }     
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Cr_distance(STRUCT_DREAM *dream, int nchains, int chain_idx, 
    int generation_idx, int cr_idx){
    
    /*######################################################################
      Purpose:
        Compute the squared normalized jumping distance.
      Input parameters:
        chains, circular history of sampled model parameters.
        nchains, the number of chains
        nparams, the number of model parameters.
        chain_idx, the index of current chain.
        generation_idx, the index of current generation.
        cr_idx, the index of CR value chosen.
        history_size, number of generations retained in the circular history.
        dream, DREAM crossover-distance accumulators.
      Output parameters:
        dream, with the selected crossover accumulator and count updated.
      Return:
        0 on success.
     ######################################################################*/
    
    int param_idx;
    int nparams = dream->nparams;
    int history_size = dream->history_size;
    double ***chains = TENSOR_DBL(dream->chains);
    double stddev[nparams], mean[nparams];

    if(Chains_STD_Single(nchains, generation_idx-1, generation_idx-1,
        nparams, chains, history_size, stddev, mean) != 0) return -1;

    int curr_slot = generation_idx%history_size;
    int prev_slot = (generation_idx-1)%history_size;
    
    for(param_idx=0; param_idx<nparams; param_idx++){ 
      if(stddev[param_idx] > 0.0){
        dream->delta[cr_idx] += pow( 
            (chains[curr_slot][chain_idx][param_idx] 
            -chains[prev_slot][chain_idx][param_idx])/stddev[param_idx], 2);
      }
    }
    
    dream->counts[cr_idx]++;

    return 0 ; 
}

/*--------------------------------------------------------------------------------*/

static int Cr_Prob(STRUCT_DREAM *dream, STRUCT_MPI *mpi){
    
    /*######################################################################
      Purpose:
        Update the sampling probability of each crossover value.
      Input parameters:
        dream, crossover statistics accumulated by the local rank.
        mpi, MPI island communicator and rank information.
      Output parameters:
        dream, synchronized crossover totals and updated probabilities.
      Return:
        0 on success.
     ######################################################################*/
    
    int i;
    double total_score=0.0;

    if(!dream || !mpi || dream->ncr < 1 || !dream->delta
        || !dream->delta_sum || !dream->delta_total || !dream->counts
        || !dream->counts_sum || !dream->counts_total
        || !dream->probabilities) return -1;
    if(MPI_Allreduce(dream->delta+1, dream->delta_sum+1, dream->ncr,
        MPI_DOUBLE, MPI_SUM, mpi->island_comm) != MPI_SUCCESS
        || MPI_Allreduce(dream->counts+1, dream->counts_sum+1, dream->ncr,
        MPI_INT, MPI_SUM, mpi->island_comm) != MPI_SUCCESS) return -1;

    if(mpi->island_rank == 0){
      for(i=1; i<=dream->ncr; i++){
        dream->delta_total[i] += dream->delta_sum[i];
        dream->counts_total[i] += dream->counts_sum[i];
        if(dream->counts_total[i] > 0){
          dream->probabilities[i] = dream->delta_total[i] 
              /dream->counts_total[i];
          total_score += dream->probabilities[i];
        }else{
          dream->probabilities[i] = 0.0;
        }
      }

      if(total_score > 0.0){
        for(i=1; i<=dream->ncr; i++){
          dream->probabilities[i] /= total_score;
        }
      }else{
        for(i=1; i<=dream->ncr; i++){
          dream->probabilities[i] = 1.0/dream->ncr;
        }
      }
    }

    if(mpi->island_size > 1){
      if(MPI_Bcast(dream->counts_total+1, dream->ncr, MPI_INT, 0,
          mpi->island_comm) != MPI_SUCCESS
          || MPI_Bcast(dream->delta_total+1, dream->ncr, MPI_DOUBLE, 0,
          mpi->island_comm) != MPI_SUCCESS
          || MPI_Bcast(dream->probabilities+1, dream->ncr, MPI_DOUBLE, 0,
          mpi->island_comm) != MPI_SUCCESS) return -1;
    }

    for(i=1; i<=dream->ncr; i++){
      dream->delta[i] = 0;
      dream->counts[i] = 0;
    }
   
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Init_Cr(STRUCT_DREAM *dream){
    
    /*######################################################################
      Purpose:
        Initialize the CR (crossover probability) values.
      Input parameters:
        dream, DREAM state with ncr already configured.
        npar, number of model parameters.
      Output parameters:
        dream, with crossover statistics and proposal work arrays allocated.
      Return:
        0 on success.
     ######################################################################*/
    
    if(!dream || dream->ncr < 1 || dream->nparams < 1) return -1;
    Free_Dream(dream);

    dream->crossover = (double *)calloc((size_t)dream->ncr+1, sizeof(double));
    dream->probabilities = (double *)calloc((size_t)dream->ncr+1,
        sizeof(double));
    dream->delta = (double *)calloc((size_t)dream->ncr+1, sizeof(double));
    dream->delta_sum = (double *)calloc((size_t)dream->ncr+1, sizeof(double));
    dream->delta_total = (double *)calloc((size_t)dream->ncr+1,
        sizeof(double));
    dream->counts = (int *)calloc((size_t)dream->ncr+1, sizeof(int));
    dream->counts_sum = (int *)calloc((size_t)dream->ncr+1, sizeof(int));
    dream->counts_total = (int *)calloc((size_t)dream->ncr+1, sizeof(int));

    if(!dream->crossover || !dream->probabilities || !dream->delta 
        || !dream->delta_sum || !dream->delta_total || !dream->counts 
        || !dream->counts_sum || !dream->counts_total){
      Free_Dream(dream);
      return -1;
    }

    int i;

    for(i=1; i<=dream->ncr; i++){
      dream->crossover[i] = ((double)(i))/dream->ncr;
      dream->delta[i] = 0;
      dream->delta_total[i] = 1;
      dream->probabilities[i] = 1.0/dream->ncr;
      dream->counts[i] = 0;
      dream->counts_total[i] = 1;
    }


    dream->jump_dims = (int *)malloc(
        sizeof(int)*(size_t)dream->nparams);
    dream->diff = (double *)malloc(
        sizeof(double)*(size_t)dream->nparams);
    dream->scale_noise = (double *)malloc(
        sizeof(double)*(size_t)dream->nparams);
    dream->additive_noise = (double *)malloc(
        sizeof(double)*(size_t)dream->nparams);

    if(!dream->jump_dims || !dream->diff || !dream->scale_noise 
        || !dream->additive_noise){
      Free_Dream(dream);
      return -1;
    }
    
    return 0 ;  
}

/*--------------------------------------------------------------------------------*/

int Free_Dream(STRUCT_DREAM *dream){
    if(!dream) return -1;
    /*######################################################################
      Purpose:
        Release all dynamically allocated storage owned by DREAM.
      Input parameters:
        dream, DREAM state whose allocated members will be released.
      Output parameters:
        dream, with owned buffers released and allocation descriptors cleared.
      Return:
        0 on success.
     ######################################################################*/


    free(dream->crossover);
    free(dream->probabilities);
    free(dream->delta);
    free(dream->delta_sum);
    free(dream->delta_total);
    free(dream->counts);
    free(dream->counts_sum);
    free(dream->counts_total);
    free(dream->jump_dims);
    free(dream->diff);
    free(dream->scale_noise);
    free(dream->additive_noise);

    dream->crossover = NULL;
    dream->probabilities = NULL;
    dream->delta = NULL;
    dream->delta_sum = NULL;
    dream->delta_total = NULL;
    dream->counts = NULL;
    dream->counts_sum = NULL;
    dream->counts_total = NULL;
    dream->jump_dims = NULL;
    dream->diff = NULL;
    dream->scale_noise = NULL;
    dream->additive_noise = NULL;
   
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int Sample_Cr(STRUCT_DREAM *dream, STRUCT_MPI *mpi){
    
    /*######################################################################
      Purpose:
        Sample a crossover index from the current categorical distribution.
      Input parameters:
        dream, crossover probabilities.
        mpi, random-number-generator state for the current rank.
      Return:
        Selected crossover index in the range 1 through ncr.
     ######################################################################*/
    
    int i;
    double tmp = RNG_UNIFORM(mpi->rank_rng);
    
    for(i=1; i<=dream->ncr; i++){
      tmp -= dream->probabilities[i];
      if(tmp<0) return i;
    }
    
    return dream->ncr; 
}

/*--------------------------------------------------------------------------------*/

static int Sample_PairNum(int max_pairs, STRUCT_MPI *mpi){
    
    /*######################################################################
      Purpose:
        Sample the number of chain pairs used to generate a jump.
      Input parameters:
        max_pairs, maximum permitted number of chain pairs.
        mpi, random-number-generator state for the current rank.
      Return:
        Number of chain pairs, uniformly sampled from 1 through max_pairs.
     ######################################################################*/
    
    int i, npairs=1;
    double tmp = RNG_UNIFORM(mpi->rank_rng);
    
    for(i=1; i<=max_pairs; i++){
      if(tmp<1.0/max_pairs){
        npairs=i;
        break;
      }else{
        tmp-=1.0/max_pairs;
      }
    }
    
    return npairs; 
}

/*--------------------------------------------------------------------------------*/

static int Rm_Outlierchain(STRUCT_MPI *mpi, STRUCT_DREAM *dream, int generation){

    /*######################################################################
      Purpose:
        Remove the outlier chains.
      Input parameters:
        mpi, MPI island layout and local chain range.
        dream, circular sample and likelihood history.
        nparams, number of model parameters.
        generation, current generation index.
      Output parameters:
        dream, with detected outliers replaced by the best chain.
      Return:
        the number of outlier chains.
      Note:
        Vrugt, j. A., et al. 2009. International Journal of Nonlinear 
          Science & Numerical Simulation, 10(3), 273-290
     ######################################################################*/
    
    if(!mpi || !dream || generation < 1 || dream->history_size < 2
        || dream->nparams < 1 || mpi->total_chains < 4
        || mpi->island_size < 1 || mpi->island_rank < 0
        || mpi->island_rank >= mpi->island_size || mpi->chain_begin < 0
        || mpi->chain_end < mpi->chain_begin
        || mpi->chain_end >= mpi->total_chains
        || !dream->chains.data_ptr || !dream->likelihood.data_ptr) return -1;

    double ***chains = TENSOR_DBL(dream->chains);
    double **likelihood = MAT_DBL(dream->likelihood);
    if(!chains || !likelihood) return -1;
    int nparams = dream->nparams;

    int Begin = generation/2;
    if(generation-Begin+1 > dream->history_size){
      Begin = generation+1-dream->history_size;
    }
    int window_length = generation-Begin+1;
    double *mean_likelihood = calloc((size_t)mpi->total_chains,
        sizeof(*mean_likelihood));
    int *indices = malloc((size_t)mpi->total_chains*sizeof(*indices));
    int allocation_status = mean_likelihood && indices ? 0 : -1;
    if(mpi->island_size > 1){
      if(MPI_Allreduce(MPI_IN_PLACE, &allocation_status, 1, MPI_INT, MPI_MIN,
          mpi->island_comm) != MPI_SUCCESS) allocation_status = -1;
    }
    if(allocation_status != 0){
      free(mean_likelihood);
      free(indices);
      return -1;
    }
    if(!mean_likelihood || !indices){
      free(mean_likelihood);
      free(indices);
      return -1;
    }

    int chain_idx, generation_idx, param_idx;

    for(chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end; chain_idx++){
	    mean_likelihood[chain_idx] = 0;
      for(generation_idx=Begin; generation_idx<=generation; generation_idx++){
        mean_likelihood[chain_idx] += 
            likelihood[generation_idx%dream->history_size][chain_idx];
      }
      mean_likelihood[chain_idx] /= window_length;
    }

    if(mpi->island_size > 1){
      if(DREAM_Allgather_Likelihood(mean_likelihood, mpi) != 0){
        free(mean_likelihood);
        free(indices);
        return -1;
      }
    }

    for(chain_idx=0; chain_idx<mpi->total_chains; chain_idx++){
      indices[chain_idx] = chain_idx;
    }
    int sort_status = 0;
    if(mpi->island_rank == 0){
      sort_status = SORT_INDICES(mean_likelihood, 
          (size_t)mpi->total_chains, indices);
    }
    if(mpi->island_size > 1){
      if(MPI_Bcast(&sort_status, 1, MPI_INT, 0, mpi->island_comm)
          != MPI_SUCCESS){
        free(mean_likelihood);
        free(indices);
        return -1;
      }
    }
    if(sort_status != 0){
      free(mean_likelihood);
      free(indices);
      return -1;
    }
    if(mpi->island_size > 1){
      if(MPI_Bcast(indices, mpi->total_chains, MPI_INT, 0, mpi->island_comm)
          != MPI_SUCCESS){
        free(mean_likelihood);
        free(indices);
        return -1;
      }
    }
    
    double Qr = mean_likelihood[indices[mpi->total_chains*3/4]] 
        -mean_likelihood[indices[mpi->total_chains/4]];
    double Omega = mean_likelihood[indices[mpi->total_chains/4]]-Qr*2;
    
    int noutliers=0;


    for(chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end; chain_idx++){
      if(mean_likelihood[chain_idx]<Omega){
        noutliers++;
               
        for(param_idx=0; param_idx<nparams; param_idx++){
          int slot = generation%dream->history_size;
          chains[slot][chain_idx][param_idx] = 
              chains[slot][indices[mpi->total_chains-1]][param_idx];
        }
        int slot = generation%dream->history_size;
        likelihood[slot][chain_idx] =
            likelihood[slot][indices[mpi->total_chains-1]];
      }        
    }

    if(mpi->island_size > 1){
      int gener_slot = generation%dream->history_size;
      if(DREAM_Allgather(chains[gener_slot], likelihood[gener_slot], mpi,
          nparams) != 0){
        free(mean_likelihood);
        free(indices);
        return -1;
      }
    }

    int total_outliers = noutliers;
    if(mpi->island_size > 1
        && MPI_Allreduce(&noutliers, &total_outliers, 1, MPI_INT, MPI_SUM,
        mpi->island_comm) != MPI_SUCCESS){
      free(mean_likelihood);
      free(indices);
      return -1;
    }

    if(total_outliers > 0){
      DREAM_VERBOSE(mpi, 3, "Generation %d: replaced %d outlier chains.", 
          generation, total_outliers);
    }
 
    free(mean_likelihood);
    free(indices);
    return total_outliers;
}

/*--------------------------------------------------------------------------------*/

static int Bounds_Enforce(double *sample, STRUCT_PARA *params){
    
    /*######################################################################
      Purpose:
        Map sampled parameters back into their configured bounds.
      Input parameters:
        sample, proposed model parameters.
        params, parameter bounds and parameter count.
      Output parameters:
        sample, bounded model parameters.
      Return:
        0 on success.
     ######################################################################*/
    
    /* Bounds, policies, and free-parameter indices are initialized and
       validated by INIT_DREAM(); only proposal values vary here. */
    if(!sample || !params) return -1;

    int param_idx;

    if(params->magnetic_mode == MAGNETIC_CARTESIAN){
      const int bx_idx = params->bx_free_index;
      const int by_idx = params->by_free_index;

      /* Resolve the 180-degree transverse-field ambiguity continuously at
         the By=0 boundary: (Bx, By) and (-Bx, -By) are equivalent. */
      if(bx_idx >= 0 && by_idx >= 0
          && fabs(params->bounds[by_idx][0]) <= 1e-12
          && params->bounds[by_idx][1] > 0.0
          && sample[by_idx] < 0.0){
        sample[by_idx] = -sample[by_idx];
        sample[bx_idx] = -sample[bx_idx];
      }
    }

    for(param_idx=0; param_idx<params->npar; param_idx++){
      const double lower = params->bounds[param_idx][0];
      const double upper = params->bounds[param_idx][1];
      const double width = upper-lower;

      if(!isfinite(sample[param_idx])) return -1;

      if(width == 0.0){
        sample[param_idx] = lower;
      }else if(sample[param_idx] < lower || sample[param_idx] > upper){
        switch(params->type[param_idx]){
          case enum_fold:{
            double offset = fmod(sample[param_idx]-lower, width);
            if(offset < 0.0) offset += width;
            sample[param_idx] = lower+offset;
            break;
          }

          case enum_reflect:{
            const double period = 2.0*width;
            double offset = fmod(sample[param_idx]-lower, period);
            if(offset < 0.0) offset += period;
            if(offset > width) offset = period-offset;
            sample[param_idx] = lower+offset;
            break;
          }

          case enum_set:
            sample[param_idx] = sample[param_idx] < lower ? lower : upper;
            break;

          default:
            return -1;
        }
      }
    }
    
    return 0;   
}

/*--------------------------------------------------------------------------------*/
