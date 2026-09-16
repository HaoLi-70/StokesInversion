#include "chain_statistics.h"
#include <math.h>
#include <stdlib.h>

/*--------------------------------------------------------------------------------*/

static void SUM_ADD(double value, double *sum, double *correction){

    /*######################################################################
      Purpose:
        Add one value using Kahan compensated summation.
      Output parameters:
        sum, accumulated sum.
        correction, accumulated rounding correction.
    ######################################################################*/

    double corrected = value-*correction;
    double updated = *sum+corrected;
    *correction = (updated-*sum)-corrected;
    *sum = updated;
}

/*--------------------------------------------------------------------------------*/

static int MPI_SUM_DOUBLE(const STRUCT_MPI *mpi, double local, double *global){

    if(!mpi || !global) return -1;
    if(mpi->island_size == 1){
      *global = local;
      return 0;
    }
    return MPI_Allreduce(&local, global, 1, MPI_DOUBLE, MPI_SUM,
        mpi->island_comm) == MPI_SUCCESS ? 0 : -1;
}

/*--------------------------------------------------------------------------------*/

static int MPI_MIN_INT(const STRUCT_MPI *mpi, int local, int *global){

    if(!mpi || !global) return -1;
    if(mpi->island_size == 1){
      *global = local;
      return 0;
    }
    return MPI_Allreduce(&local, global, 1, MPI_INT, MPI_MIN,
        mpi->island_comm) == MPI_SUCCESS ? 0 : -1;
}

/*--------------------------------------------------------------------------------*/

static bool MPI_STATE_VALID(const STRUCT_MPI *mpi){

    return mpi && !mpi->is_master && mpi->island_comm != MPI_COMM_NULL
        && mpi->island_size >= 1 && mpi->island_rank >= 0
        && mpi->island_rank < mpi->island_size && mpi->total_chains >= 2
        && mpi->chains_per_rank >= 1 && mpi->chain_begin >= 0
        && mpi->chain_end >= mpi->chain_begin
        && mpi->chain_end-mpi->chain_begin+1 == mpi->chains_per_rank
        && mpi->chain_end < mpi->total_chains;
}

/*--------------------------------------------------------------------------------*/

int PSRF(STRUCT_MPI *mpi, double ***chains, int generation,
    const int *param_indices, int nparams, int chain_nparams,
    int history_size, double *rhat){

    /*######################################################################
      Purpose:
        Compute the Gelman-Rubin potential scale reduction factor.
      Input parameters:
        mpi, MPI island and local-chain layout.
        chains, circular chain-history buffer.
        generation, current absolute generation index.
        param_indices, selected chain-parameter indices.
        nparams, number of selected parameters.
        chain_nparams, total parameters stored in each chain.
        history_size, number of retained circular-buffer slots.
      Output parameters:
        rhat, PSRF value for every selected parameter.
      Return:
        0 on success; -1 for invalid input, non-finite samples, allocation
          failure, or an MPI failure.
      Note:
        Brooks and Gelman (1998); Gelman and Rubin (1992).
    ######################################################################*/

    if(!MPI_STATE_VALID(mpi) || !chains || !param_indices || !rhat
        || generation < 1 || nparams < 1 || chain_nparams < 1
        || history_size < 2) return -1;
    for(int selected=0; selected<nparams; selected++){
      if(param_indices[selected] < 0
          || param_indices[selected] >= chain_nparams) return -1;
      rhat[selected] = INFINITY;
    }

    int begin = generation/2;
    int window_length = generation-begin+1;
    if(window_length > 1500){
      window_length = 1500;
      begin = generation+1-window_length;
    }
    if(window_length < 2 || window_length > history_size) return -1;

    double *chain_mean = malloc((size_t)mpi->chains_per_rank
        *sizeof(*chain_mean));
    int local_status = chain_mean ? 0 : -1;
    int global_status = 0;
    if(MPI_MIN_INT(mpi, local_status, &global_status) != 0
        || global_status != 0){
      free(chain_mean);
      return -1;
    }

    for(int selected=0; selected<nparams; selected++){
      int param_idx = param_indices[selected];
      double local_mean_sum = 0.0, local_mean_correction = 0.0;
      local_status = 0;
      for(int chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end;
          chain_idx++){
        int local_chain = chain_idx-mpi->chain_begin;
        double sum = 0.0, correction = 0.0;
        for(int gener=begin; gener<=generation; gener++){
          double value = chains[gener%history_size][chain_idx][param_idx];
          if(!isfinite(value)){
            local_status = -1;
            value = 0.0;
          }
          SUM_ADD(value, &sum, &correction);
        }
        chain_mean[local_chain] = sum/(double)window_length;
        SUM_ADD(chain_mean[local_chain], &local_mean_sum,
            &local_mean_correction);
      }
      if(MPI_MIN_INT(mpi, local_status, &global_status) != 0
          || global_status != 0){
        free(chain_mean);
        return -1;
      }

      double total_mean_sum = 0.0;
      if(MPI_SUM_DOUBLE(mpi, local_mean_sum, &total_mean_sum) != 0
          || !isfinite(total_mean_sum)){
        free(chain_mean);
        return -1;
      }
      double mean = total_mean_sum/(double)mpi->total_chains;

      double local_between = 0.0, correction = 0.0;
      for(int local_chain=0; local_chain<mpi->chains_per_rank; local_chain++){
        double diff = chain_mean[local_chain]-mean;
        SUM_ADD(diff*diff, &local_between, &correction);
      }
      double between = 0.0;
      if(MPI_SUM_DOUBLE(mpi, local_between, &between) != 0
          || !isfinite(between) || between < 0.0){
        free(chain_mean);
        return -1;
      }
      between /= (double)(mpi->total_chains-1);

      double local_within = 0.0;
      correction = 0.0;
      for(int chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end;
          chain_idx++){
        int local_chain = chain_idx-mpi->chain_begin;
        for(int gener=begin; gener<=generation; gener++){
          double diff = chains[gener%history_size][chain_idx][param_idx]
              -chain_mean[local_chain];
          SUM_ADD(diff*diff, &local_within, &correction);
        }
      }
      double within = 0.0;
      if(MPI_SUM_DOUBLE(mpi, local_within, &within) != 0
          || !isfinite(within) || within < 0.0){
        free(chain_mean);
        return -1;
      }
      within /= (double)(window_length-1)*(double)mpi->total_chains;
      if(within == 0.0) continue;

      double variance = (double)(window_length-1)/(double)window_length
          *within+between;
      double rhat_squared = ((double)mpi->total_chains+1.0)
          /(double)mpi->total_chains*variance/within
          -(double)(window_length-1)/(double)window_length
          /(double)mpi->total_chains;
      if(!isfinite(rhat_squared) || rhat_squared < 0.0){
        free(chain_mean);
        return -1;
      }
      rhat[selected] = fmax(1.0, sqrt(rhat_squared));
    }

    free(chain_mean);
    return 0;
}

/*--------------------------------------------------------------------------------*/

int Chains_STD(STRUCT_MPI *mpi, int begin_generation, int end_generation,
    int nparams, double ***chains, int history_size, double *stddev,
    double *mean){

    /*######################################################################
      Purpose:
        Compute distributed means and sample standard deviations over a
          retained circular-history window.
      Input parameters:
        mpi, MPI island and local-chain layout.
        begin_generation, first absolute generation index, inclusive.
        end_generation, last absolute generation index, inclusive.
        nparams, number of chain parameters.
        chains, circular chain-history buffer.
        history_size, number of retained circular-buffer slots.
      Output parameters:
        stddev, sample standard deviation of every parameter.
        mean, mean value of every parameter.
      Return:
        0 on success; -1 for invalid input, non-finite samples, insufficient
          history, or an MPI failure.
    ######################################################################*/

    if(!MPI_STATE_VALID(mpi) || !chains || !stddev || !mean || nparams < 1
        || history_size < 1 || begin_generation < 0
        || end_generation < begin_generation) return -1;

    long long window_length = (long long)end_generation-begin_generation+1;
    if(window_length > history_size) return -1;
    double sample_count = (double)window_length*(double)mpi->total_chains;
    if(!isfinite(sample_count) || sample_count <= 1.0) return -1;

    for(int param_idx=0; param_idx<nparams; param_idx++){
      double local_sum = 0.0, correction = 0.0;
      int local_status = 0;
      for(int gener=begin_generation; gener<=end_generation; gener++){
        int slot = gener%history_size;
        for(int chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end;
            chain_idx++){
          double value = chains[slot][chain_idx][param_idx];
          if(!isfinite(value)){
            local_status = -1;
            value = 0.0;
          }
          SUM_ADD(value, &local_sum, &correction);
        }
      }
      int global_status = 0;
      if(MPI_MIN_INT(mpi, local_status, &global_status) != 0
          || global_status != 0) return -1;

      double total_sum = 0.0;
      if(MPI_SUM_DOUBLE(mpi, local_sum, &total_sum) != 0
          || !isfinite(total_sum)) return -1;
      mean[param_idx] = total_sum/sample_count;
      if(!isfinite(mean[param_idx])) return -1;

      double local_squared_sum = 0.0;
      correction = 0.0;
      for(int gener=begin_generation; gener<=end_generation; gener++){
        int slot = gener%history_size;
        for(int chain_idx=mpi->chain_begin; chain_idx<=mpi->chain_end;
            chain_idx++){
          double diff = chains[slot][chain_idx][param_idx]-mean[param_idx];
          double squared = diff*diff;
          if(!isfinite(squared)) local_status = -1;
          else SUM_ADD(squared, &local_squared_sum, &correction);
        }
      }
      if(MPI_MIN_INT(mpi, local_status, &global_status) != 0
          || global_status != 0) return -1;

      double total_squared_sum = 0.0;
      if(MPI_SUM_DOUBLE(mpi, local_squared_sum, &total_squared_sum) != 0
          || !isfinite(total_squared_sum) || total_squared_sum < 0.0){
        return -1;
      }
      stddev[param_idx] = sqrt(total_squared_sum/(sample_count-1.0));
      if(!isfinite(stddev[param_idx])) return -1;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

int Chains_STD_Single(int nchains, int begin_generation, int end_generation,
    int nparams, double ***chains, int history_size, double *stddev,
    double *mean){

    /*######################################################################
      Purpose:
        Compute means and sample standard deviations over all chains without
          an MPI reduction.
      Input parameters:
        nchains, total number of chains available in the local buffer.
        begin_generation, first absolute generation index, inclusive.
        end_generation, last absolute generation index, inclusive.
        nparams, number of chain parameters.
        chains, circular chain-history buffer.
        history_size, number of retained circular-buffer slots.
      Output parameters:
        stddev, sample standard deviation of every parameter.
        mean, mean value of every parameter.
      Return:
        0 on success; -1 for invalid input, non-finite samples, or
          insufficient history.
    ######################################################################*/

    if(nchains < 2 || begin_generation < 0
        || end_generation < begin_generation || nparams < 1 || !chains
        || history_size < 1 || !stddev || !mean) return -1;
    long long window_length = (long long)end_generation-begin_generation+1;
    if(window_length > history_size) return -1;
    double sample_count = (double)window_length*(double)nchains;
    if(!isfinite(sample_count) || sample_count <= 1.0) return -1;

    for(int param_idx=0; param_idx<nparams; param_idx++){
      double sum = 0.0, correction = 0.0;
      for(int gener=begin_generation; gener<=end_generation; gener++){
        int slot = gener%history_size;
        for(int chain_idx=0; chain_idx<nchains; chain_idx++){
          double value = chains[slot][chain_idx][param_idx];
          if(!isfinite(value)) return -1;
          SUM_ADD(value, &sum, &correction);
        }
      }
      mean[param_idx] = sum/sample_count;
      if(!isfinite(mean[param_idx])) return -1;

      double squared_sum = 0.0;
      correction = 0.0;
      for(int gener=begin_generation; gener<=end_generation; gener++){
        int slot = gener%history_size;
        for(int chain_idx=0; chain_idx<nchains; chain_idx++){
          double diff = chains[slot][chain_idx][param_idx]-mean[param_idx];
          double squared = diff*diff;
          if(!isfinite(squared)) return -1;
          SUM_ADD(squared, &squared_sum, &correction);
        }
      }
      if(!isfinite(squared_sum) || squared_sum < 0.0) return -1;
      stddev[param_idx] = sqrt(squared_sum/(sample_count-1.0));
      if(!isfinite(stddev[param_idx])) return -1;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/
