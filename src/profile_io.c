
#include "profile_io.h"
#include "logger.h"
#include <ctype.h>
#include <errno.h>
#include <fitsio.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

    /*######################################################################
      Module:
        FITS and DAT profile input, wavelength and noise handling, inversion
        cache management, distributed samples, and FITS result output.
    ######################################################################*/

/*--------------------------------------------------------------------------------*/

#define IO_VERBOSE(lvl, ...)                                           \
    do{                                                                \
      if(mpi->verbose_level>=lvl){                                     \
        snprintf(message_buffer, sizeof(message_buffer), __VA_ARGS__); \
        LOG_WRITE(message_buffer, true, true);                         \
      }                                                                \
    }while(0)

#define IO_ERROR(...) IO_FAIL(__VA_ARGS__)

#define CACHE_MAGIC "mciv"

static FILE *cache_file = NULL;
static fitsfile *obs_file = NULL;
static fitsfile *res_file = NULL;
static long fits_first_pixel[4];

/*--------------------------------------------------------------------------------*/

static _Noreturn void IO_FAIL(const char *format, ...){

    const char prefix[] = "-ERROR in routine IO: ";
    size_t offset = sizeof(prefix)-1U;

    memcpy(message_buffer, prefix, offset);
    va_list args;
    va_start(args, format);
    vsnprintf(message_buffer+offset, sizeof(message_buffer)-offset,
        format, args);
    va_end(args);

    LOG_WRITE(message_buffer, true, true);
    CLOSE_FILES();
    ABORTED();
    abort();
}

/*--------------------------------------------------------------------------------*/

static void MODEL_PARAMETER_KEY(int index, char key[9]){

    /*######################################################################
      Purpose:
        Build a classic FITS keyword for one model parameter.
      Input parameters:
        index, model parameter index in the range 0 through 99999.
      Output parameters:
        key, null-terminated keyword PAR000 through PAR99999.
    ######################################################################*/

    if(index < 0 || index > 99999){
      IO_ERROR("model parameter index %d cannot be stored in a FITS "
          "keyword.\n", index);
    }
    snprintf(key, 9, "PAR%03u", (unsigned int)index);
}

/*--------------------------------------------------------------------------------*/

static uint64_t HASH_BYTES(uint64_t hash, const void *data, size_t size){

    /*######################################################################
      Purpose:
        Add a byte sequence to a 64-bit FNV-1a hash.
      Input parameters:
        hash, current hash value.
        data, bytes to include.
        size, number of bytes in data.
      Return:
        Updated hash value.
    ######################################################################*/

    const unsigned char *bytes = (const unsigned char *)data;
    for(size_t i=0; i<size; i++){
      hash ^= bytes[i];
      hash *= UINT64_C(1099511628211);
    }
    return hash;
}

/*--------------------------------------------------------------------------------*/

static uint64_t HASH_FILE(uint64_t hash, const char *path){

    if(!path || path[0] == '\0') return HASH_BYTES(hash, "", 1);
    FILE *file = fopen(path, "rb");
    if(!file) return HASH_BYTES(hash, "FILE_READ_ERROR", 16);

    unsigned char buffer[65536];
    size_t count;
    while((count = fread(buffer, 1, sizeof(buffer), file)) > 0){
      hash = HASH_BYTES(hash, buffer, count);
    }
    int read_error = ferror(file);
    if(fclose(file) != 0 || read_error){
      hash = HASH_BYTES(hash, "FILE_READ_ERROR", 16);
    }
    return hash;
}

/*--------------------------------------------------------------------------------*/

static uint64_t CONFIG_HASH(const STRUCT_PROFILE_IO *input, 
    const STRUCT_PARA *params, const STRUCT_DREAM *dream, 
    const STRUCT_STK *stokes){

    /*######################################################################
      Purpose:
        Hash all inputs that affect cached inversion results and samples.
      Input parameters:
        input, file paths, dimensions, and selected solution box.
        params, model layout, bounds, and spectral-line metadata.
        dream, DREAM and sample-output controls.
        stokes, noise, wavelength, and spectral-region controls.
      Return:
        A deterministic 64-bit configuration hash.
    ######################################################################*/

    uint64_t hash = UINT64_C(1469598103934665603);
    hash = HASH_BYTES(hash, input->data_path, strlen(input->data_path)+1);
    hash = HASH_BYTES(hash, input->wavelength_path, 
        strlen(input->wavelength_path)+1);
    hash = HASH_BYTES(hash, input->result_path, strlen(input->result_path)+1);
    hash = HASH_BYTES(hash, input->sample_path, strlen(input->sample_path)+1);
    hash = HASH_BYTES(hash, &input->profile_format, 
        sizeof(input->profile_format));
    hash = HASH_FILE(hash, input->data_path);
    hash = HASH_FILE(hash, input->wavelength_path);
    hash = HASH_BYTES(hash, input->sol_box, sizeof(input->sol_box));
    hash = HASH_BYTES(hash, &params->magnetic_mode, 
        sizeof(params->magnetic_mode));
    hash = HASH_BYTES(hash, &params->nmodel, sizeof(params->nmodel));
    for(int i=0; i<params->nmodel; i++){
      hash = HASH_BYTES(hash, &params->inv[i], sizeof(params->inv[i]));
      hash = HASH_BYTES(hash, &params->value_const[i],
          sizeof(params->value_const[i]));
      hash = HASH_BYTES(hash, params->limits[i], sizeof(params->limits[i]));
      hash = HASH_BYTES(hash, &params->kind[i], sizeof(params->kind[i]));
    }
    for(int i=0; i<params->nlines; i++){
      const STRUCT_MELINE *line = &params->lines[i];
      hash = HASH_BYTES(hash, &line->wavelength0, sizeof(line->wavelength0));
      hash = HASH_BYTES(hash, &line->lande_factor, sizeof(line->lande_factor));
      hash = HASH_BYTES(hash, &line->zeeman_shift, sizeof(line->zeeman_shift));
      hash = HASH_BYTES(hash, &line->iw_begin, sizeof(line->iw_begin));
      hash = HASH_BYTES(hash, &line->iw_end, sizeof(line->iw_end));
      hash = HASH_BYTES(hash, &line->dopp_index, sizeof(line->dopp_index));
      hash = HASH_BYTES(hash, &line->damp_index, sizeof(line->damp_index));
      hash = HASH_BYTES(hash, &line->eta_index, sizeof(line->eta_index));
      hash = HASH_BYTES(hash, line->custom_limits, sizeof(line->custom_limits));
      hash = HASH_BYTES(hash, line->custom_value, sizeof(line->custom_value));
      hash = HASH_BYTES(hash, line->custom_bounds, sizeof(line->custom_bounds));
      hash = HASH_BYTES(hash, line->custom_inversion,
          sizeof(line->custom_inversion));
      hash = HASH_BYTES(hash, line->custom_inv, sizeof(line->custom_inv));
    }
    for(int i=0; i<params->nregions; i++){
      const STRUCT_REGION *region = &params->regions[i];
      hash = HASH_BYTES(hash, &region->wavelength_min,
          sizeof(region->wavelength_min));
      hash = HASH_BYTES(hash, &region->wavelength_max,
          sizeof(region->wavelength_max));
      hash = HASH_BYTES(hash, &region->iw_begin, sizeof(region->iw_begin));
      hash = HASH_BYTES(hash, &region->iw_end, sizeof(region->iw_end));
      hash = HASH_BYTES(hash, &region->continuum_index,
          sizeof(region->continuum_index));
      hash = HASH_BYTES(hash, &region->beta_index, sizeof(region->beta_index));
      hash = HASH_BYTES(hash, region->custom_limits,
          sizeof(region->custom_limits));
      hash = HASH_BYTES(hash, region->custom_value,
          sizeof(region->custom_value));
      hash = HASH_BYTES(hash, region->custom_bounds,
          sizeof(region->custom_bounds));
      hash = HASH_BYTES(hash, region->custom_inversion,
          sizeof(region->custom_inversion));
      hash = HASH_BYTES(hash, region->custom_inv, sizeof(region->custom_inv));
    }
    hash = HASH_BYTES(hash, &stokes->noise_mode, sizeof(stokes->noise_mode));
    if(stokes->noise_mode == NOISE_FROM_INTENSITY){
      hash = HASH_BYTES(hash, stokes->noise_level,
          sizeof(stokes->noise_level));
    }
    hash = HASH_BYTES(hash, &stokes->region_gap, 
        sizeof(stokes->region_gap));
    hash = HASH_BYTES(hash, &stokes->line_wing_widths, 
        sizeof(stokes->line_wing_widths));
    hash = HASH_BYTES(hash, &dream->max_memory_gb, 
        sizeof(dream->max_memory_gb));
    hash = HASH_BYTES(hash, &dream->max_sample_file_gb,
        sizeof(dream->max_sample_file_gb));
    hash = HASH_BYTES(hash, &stokes->nw, sizeof(stokes->nw));
    hash = HASH_BYTES(hash, stokes->wavelength, 
        (size_t)stokes->nw*sizeof(*stokes->wavelength));
    hash = HASH_BYTES(hash, &dream->nchains, sizeof(dream->nchains));
    hash = HASH_BYTES(hash, &dream->burnin_generations,
        sizeof(dream->burnin_generations));
    hash = HASH_BYTES(hash, &dream->sampling_generations,
        sizeof(dream->sampling_generations));
    hash = HASH_BYTES(hash, &dream->max_pairs, sizeof(dream->max_pairs));
    hash = HASH_BYTES(hash, &dream->ncr, sizeof(dream->ncr));
    hash = HASH_BYTES(hash, &dream->update_crossover, 
        sizeof(dream->update_crossover));
    hash = HASH_BYTES(hash, &dream->update_proposal_noise,
        sizeof(dream->update_proposal_noise));
    hash = HASH_BYTES(hash, &dream->sample_output, 
        sizeof(dream->sample_output));

    return hash;
}

/*--------------------------------------------------------------------------------*/

void MODEL_PARAMETER_NAME(const STRUCT_PARA *params, int index, 
    char *name, size_t size){

    /*######################################################################
      Purpose:
        Build the FITS column name for one dynamic model parameter.
      Input parameters:
        params, model layout and magnetic parameterization.
        index, model parameter index.
        size, capacity of the destination string.
      Output parameters:
        name, null-terminated parameter name.
    ######################################################################*/

    switch(params->kind[index]){
      case PARAM_B0:
        snprintf(name, size, "%s", params->magnetic_mode == MAGNETIC_CARTESIAN 
            ? "Bz" : "Bmod");
        return;
      case PARAM_B1:
        snprintf(name, size, "%s", params->magnetic_mode == MAGNETIC_CARTESIAN 
            ? "Bx" : "ThetaB");
        return;
      case PARAM_B2:
        snprintf(name, size, "%s", params->magnetic_mode == MAGNETIC_CARTESIAN 
            ? "By" : "PhiB");
        return;
      case PARAM_VLOS:
        snprintf(name, size, "Vlos");
        return;
      default:
        break;
    }
    for(int iline=0; iline<params->nlines; iline++){
      const STRUCT_MELINE *line = &params->lines[iline];
      if(index == line->dopp_index){
        snprintf(name, size, "line%d_Dopp", iline);
        return;
      }
      if(index == line->damp_index){
        snprintf(name, size, "line%d_Damp", iline);
        return;
      }
      if(index == line->eta_index){
        snprintf(name, size, "line%d_Eta", iline);
        return;
      }
    }
    for(int iregion=0; iregion<params->nregions; iregion++){
      const STRUCT_REGION *region = &params->regions[iregion];
      if(index == region->continuum_index){
        snprintf(name, size, "region%d_Continuum", iregion);
        return;
      }
      if(index == region->beta_index){
        snprintf(name, size, "region%d_Beta", iregion);
        return;
      }
    }
    snprintf(name, size, "parameter%d", index);
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_DAT_ROW(const char *text, double values[5]){

    /*######################################################################
      Purpose:
        Parse one DAT row containing wavelength and four Stokes values.
      Input parameters:
        text, row text after comma normalization.
      Output parameters:
        values, wavelength followed by I, Q, U, and V.
      Return:
        true for exactly five finite values; false otherwise.
    ######################################################################*/

    const char *cursor = text;
    for(int i=0; i<5; i++){
      while(isspace((unsigned char)*cursor)) cursor++;
      char *end = NULL;
      errno = 0;
      values[i] = strtod(cursor, &end);
      if(end == cursor || errno == ERANGE || !isfinite(values[i])) return false;
      cursor = end;
      if(i < 4 && !isspace((unsigned char)*cursor)) return false;
    }
    while(isspace((unsigned char)*cursor)) cursor++;
    return *cursor == '\0' || *cursor == '#';
}

/*--------------------------------------------------------------------------------*/

static int READ_DAT_PROFILE(const char *path, STRUCT_STK *stokes){

    /*######################################################################
      Purpose:
        Read one wavelength grid and Stokes profile from a text file.
      Input parameters:
        path, path to a file containing wavelength, I, Q, U, and V columns.
      Output parameters:
        stokes, allocated wavelength and profile arrays.
      Return:
        0 on success; -1 for an invalid file or allocation failure.
    ######################################################################*/

    FILE *file = fopen(path, "r");
    if(!file) return -1;

    char *line = NULL;
    size_t capacity = 0;
    int count = 0;
    ptrdiff_t read_status;
    while((read_status=STR_READ_LINE(&line, &capacity, file)) > 0){
      char *text = line;
      while(isspace((unsigned char)*text)) text++;
      if(*text == '\0' || *text == '#') continue;
      for(char *cursor=text; *cursor; cursor++){
        if(*cursor == ',') *cursor = ' ';
      }
      double values[5];
      if(!PARSE_DAT_ROW(text, values)){
        free(line);
        fclose(file);
        return -1;
      }
      if(count >= INT_MAX/4){
        free(line);
        fclose(file);
        return -1;
      }
      count++;
    }

    if(read_status < -1 || count < 2){
      free(line);
      fclose(file);
      return -1;
    }

    rewind(file);
    stokes->wavelength = malloc((size_t)count*sizeof(*stokes->wavelength));
    stokes->profile = malloc((size_t)count*4*sizeof(*stokes->profile));
    if(!stokes->wavelength || !stokes->profile){
      free(stokes->wavelength);
      free(stokes->profile);
      stokes->wavelength = NULL;
      stokes->profile = NULL;
      free(line);
      fclose(file);
      return -1;
    }

    int iw = 0;
    while((read_status=STR_READ_LINE(&line, &capacity, file)) > 0){
      char *text = line;
      while(isspace((unsigned char)*text)) text++;
      if(*text == '\0' || *text == '#') continue;
      for(char *cursor=text; *cursor; cursor++){
        if(*cursor == ',') *cursor = ' ';
      }
      double values[5];
      if(!PARSE_DAT_ROW(text, values)) break;
      stokes->wavelength[iw] = values[0];
      stokes->profile[iw] = values[1];
      stokes->profile[count+iw] = values[2];
      stokes->profile[2*count+iw] = values[3];
      stokes->profile[3*count+iw] = values[4];
      iw++;
    }
    free(line);
    fclose(file);
    stokes->nw = count;
    return read_status >= -1 && iw == count ? 0 : -1;
}

/*--------------------------------------------------------------------------------*/

static int CLOSE_FITS_FILE(fitsfile **file_ptr){

    /*######################################################################
      Purpose:
        Close a FITS file and clear its handle.
      Input parameters:
        file_ptr, address of the FITS file handle.
      Output parameters:
        file_ptr, set to NULL after a successful close.
      Return:
        CFITSIO status code.
    ######################################################################*/

    if(*file_ptr){
      int status = 0;
      fits_close_file(*file_ptr, &status);
      *file_ptr = NULL;
      if(status){ 
        LOG_ERROR(ERR_LVL_WARNING, "CLOSE_FITS_FILE", 
          "Error closing FITS file \n");
        return status;
      }
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

static int CLOSE_STD_FILE(FILE **file_ptr){

    /*######################################################################
      Purpose:
        Close a standard C file and clear its handle.
      Input parameters:
        file_ptr, address of the file handle.
      Output parameters:
        file_ptr, set to NULL after closing.
      Return:
        0 on success; EOF when closing or flushing fails.
    ######################################################################*/

    if(!file_ptr || !*file_ptr) return 0;
    int status = fclose(*file_ptr);
    *file_ptr = NULL;
    return status;
}

/*--------------------------------------------------------------------------------*/

int READ_WAVELENGTH(STRUCT_PROFILE_IO *input, STRUCT_STK *stokes, 
    STRUCT_MPI *mpi){

    /*######################################################################
      Purpose:
        Read the wavelength grid and initialize profile dimensions.
      Input parameters:
        input, paths and requested solution region.
        stokes, Stokes-profile state to initialize.
        mpi, MPI communicator and rank information.
      Output parameters:
        input, validated data dimensions and solution region.
        stokes, allocated wavelength grid and profile dimensions.
      Return:
        0 on success.
    ######################################################################*/

    // CFITSIO status value MUST be initialized to zero!
    // number of the HDUs, the bit of each pixel, the number of axes
    int status = 0, bitpix, naxis;
    // the length of each axis, the first pixel
    long naxes[4];

    fits_first_pixel[0] = 1;
    fits_first_pixel[1] = 1;

    fitsfile *fptr_wav = NULL;
  
    if(mpi->world_rank==0){

      IO_VERBOSE(2, "\n -- reading the wavelength -- \n");

      if(input->profile_format == PROFILE_FORMAT_DAT){
        status = READ_DAT_PROFILE(input->data_path, stokes);
        if(status) IO_ERROR("invalid DAT profile; expected columns: " 
            "wavelength I Q U V.\n");
        naxis = 3;
        input->nx = 1;
        input->ny = 1;
      }else{

      fits_open_file(&obs_file, input->data_path, 
          READONLY, &status);
      if(status) IO_ERROR("error in opening the data file: " 
          "status = %d \n", status);   

      fits_open_file(&fptr_wav, input->wavelength_path, READONLY, &status);
      if(status) IO_ERROR("error in opening the wavelength file: " 
          "status = %d \n", status);   

      // get the dimension of the data hdu
      fits_get_img_dim(obs_file, &naxis, &status);
      if(status) IO_ERROR("error in getting the dimension of "
          "the data hdu: status = %d \n", status);   

      if(naxis<3 || naxis>4){ 
        IO_ERROR("wrong dimension of the data hdu.\n"); 
      }

      // get the size of each dimension 
      fits_get_img_size(obs_file, naxis, naxes, &status);
      if(status) IO_ERROR("error in getting the size of the data hdu: " 
          "status = %d \n", status);   
      if(naxes[1]!=4) IO_ERROR("wrong size of the data hdu.\n");         

      if(naxes[0] < 1 || naxes[0] > INT_MAX/4 
          || naxes[2] < 1 || naxes[2] > INT_MAX 
          || (naxis == 4 && (naxes[3] < 1 || naxes[3] > INT_MAX))){
        IO_ERROR("profile dimensions exceed the supported integer range " 
            "(nw must be at most INT_MAX/4).\n");
      }
      stokes->nw = (int)naxes[0];
      input->nx = (int)naxes[2];
      if(naxis==3){
        input->ny = 1;        
      }else{
        input->ny = (int)naxes[3];
      }
      if((long long)input->nx*input->ny > INT_MAX){
        IO_ERROR("the number of image pixels exceeds the supported range.\n");
      }

      if(stokes->noise_mode != NOISE_FROM_INTENSITY){
        int nhdus = 0;
        fits_get_num_hdus(obs_file, &nhdus, &status);
        if(status || nhdus < 2) IO_ERROR(
            "noise mode requires an image in HDU 2.\n");
        fits_movabs_hdu(obs_file, 2, NULL, &status);
        fits_get_img_dim(obs_file, &naxis, &status);
        if(status || naxis < 1 || naxis > 4) IO_ERROR(
            "noise HDU must have between one and four dimensions.\n");
        fits_get_img_size(obs_file, naxis, naxes, &status);
        if(status) IO_ERROR("error reading the noise HDU dimensions.\n");

        if(stokes->noise_mode == NOISE_GLOBAL){
          if(naxis != 2 || naxes[0] != stokes->nw || naxes[1] != 4){
            IO_ERROR("GLOBAL noise HDU must have dimensions [nw, 4].\n");
          }
          stokes->noise = malloc((size_t)stokes->nw*4U
              *sizeof(*stokes->noise));
          if(!stokes->noise){
            IO_ERROR("cannot allocate the global noise array.\n");
          }
          long first[2] = {1, 1};
          fits_read_pix(obs_file, TDOUBLE, first, stokes->nw*4, NULL,
              stokes->noise, NULL, &status);
          if(status) IO_ERROR("error reading the global noise HDU.\n");
          for(int index=0; index<stokes->nw*4; index++){
            if(!isfinite(stokes->noise[index])
                || stokes->noise[index] <= 0.0){
              IO_ERROR("GLOBAL noise values must be finite and positive.\n");
            }
          }
        }else{
          bool valid_3d = naxis == 3 && input->ny == 1 
              && naxes[0] == stokes->nw && naxes[1] == 4 
              && naxes[2] == input->nx;
          bool valid_4d = naxis == 4 && naxes[0] == stokes->nw 
              && naxes[1] == 4 && naxes[2] == input->nx 
              && naxes[3] == input->ny;
          if(!valid_3d && !valid_4d) IO_ERROR(
              "PIXEL noise HDU must match the profile dimensions.\n");
        }
        if(status) IO_ERROR("error reading the noise HDU.\n");
        fits_movabs_hdu(obs_file, 1, NULL, &status);
      }
      }
    }

    MPI_Bcast(&input->nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&input->ny, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // check the solution box
    if(input->sol_box[0][0] < 0) input->sol_box[0][0] = 0;
    if(input->sol_box[0][0] >= input->nx) input->sol_box[0][0] = input->nx-1;
    if(input->sol_box[0][1] < 0) input->sol_box[0][1] = input->nx-1;
    if(input->sol_box[0][1] >= input->nx) input->sol_box[0][1] = input->nx-1;
    if(input->sol_box[1][0] < 0) input->sol_box[1][0] = 0;
    if(input->sol_box[1][0] >= input->ny) input->sol_box[1][0] = input->ny-1;
    if(input->sol_box[1][1] < 0) input->sol_box[1][1] = input->ny-1;
    if(input->sol_box[1][1] >= input->ny) input->sol_box[1][1] = input->ny-1;
    if(input->sol_box[0][0] > input->sol_box[0][1] 
        || input->sol_box[1][0] > input->sol_box[1][1]){
      if(mpi->world_rank == 0) IO_ERROR(
          "sol_box start coordinates must not exceed the end coordinates.\n");
      return -1;
    }

    if(mpi->world_rank==0){

      if(input->profile_format == PROFILE_FORMAT_FITS){
      // get the type of the hdu
      fits_get_img_type(obs_file, &bitpix, &status);
      if(status) IO_ERROR("error in getting the data type " 
          "of data hdu: status = %d \n", status);   
      
      if(bitpix != SHORT_IMG && bitpix != FLOAT_IMG 
          && bitpix != DOUBLE_IMG){
        IO_ERROR("error in the data type of data hdu: " 
            "bitpix = %d.\n", bitpix);
      }
      
      // get the dimension of the wavlength hdu
      fits_get_img_dim(fptr_wav, &naxis, &status);
      if(status) IO_ERROR("error in getting the dimension " 
          "of wavelength hdu: status = %d \n", status);   
      if(naxis!=1) IO_ERROR("error in the dimension of " 
          "wavelength hdu.\n");

      // get the size of the wavelength 
      fits_get_img_size(fptr_wav, naxis, naxes, &status);
      if(status) IO_ERROR("error in getting the size of " 
          "the wavelength: status = %d \n", status);   
      if(stokes->nw!=naxes[0]) IO_ERROR("wavelength is not " 
          "consistent with the profile.\n");
      
      // get the type of the hdu
      fits_get_img_type(fptr_wav, &bitpix, &status);
      if(status) IO_ERROR("error in getting the data type of " 
          "the wavelength hdu: status = %d \n", status);   

      stokes->wavelength = malloc((size_t)stokes->nw
          *sizeof(*stokes->wavelength));
      if(!stokes->wavelength) IO_ERROR("cannot allocate wavelength array.\n");

      if(bitpix==-32 || bitpix==-64){
        fits_read_pix(fptr_wav, TDOUBLE, fits_first_pixel, stokes->nw, NULL, 
            stokes->wavelength, NULL, &status);     
      }else{
        free(stokes->wavelength);
        IO_ERROR("error in the data type of the wavelength hdu.\n");
      }
      if(status) IO_ERROR("error in reading the wavelength: " 
          "status = %d \n", status);   

      status = CLOSE_FITS_FILE(&fptr_wav);
      fptr_wav = NULL;
      if(status) IO_ERROR("error in closing the wavelength file: " 
          "status = %d \n", status);   

      status = CLOSE_FITS_FILE(&obs_file);
      if(status) IO_ERROR("error in closing the data file: " 
          "status = %d \n", status);
      }

      if(stokes->nw < 5){
        IO_ERROR("at least five wavelength samples are required.\n");
      }
      for(int iw=0; iw<stokes->nw; iw++){
        if(!isfinite(stokes->wavelength[iw])){
          IO_ERROR("wavelength samples must be finite.\n");
        }
        if(iw > 0 && stokes->wavelength[iw] <= stokes->wavelength[iw-1]){
          IO_ERROR("wavelength samples must be strictly increasing.\n");
        }
      }

      // broadcast to slaves
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&stokes->nw, 1, MPI_INT, 0, MPI_COMM_WORLD);
      // broadcast the wavelengths
      MPI_Bcast(stokes->wavelength, stokes->nw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(stokes->noise_mode == NOISE_FROM_INTENSITY){
        MPI_Bcast(stokes->noise_level, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }else if(stokes->noise_mode == NOISE_GLOBAL){
        MPI_Bcast(stokes->noise, stokes->nw*4, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
      }

    }else{

      // if error in reading wavelength
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(status) ABORTED();

      //  receive the size of the wavelength
      MPI_Bcast(&stokes->nw, 1, MPI_INT, 0, MPI_COMM_WORLD);
      stokes->wavelength = malloc((size_t)stokes->nw
          *sizeof(*stokes->wavelength));
      if(!stokes->wavelength) IO_ERROR("cannot allocate wavelength array.\n");
      if(stokes->noise_mode == NOISE_GLOBAL){
        stokes->noise = malloc((size_t)stokes->nw*4U
            *sizeof(*stokes->noise));
        if(!stokes->noise) IO_ERROR("cannot allocate global noise array.\n");
      }

      // receive the wavelength
      MPI_Bcast(stokes->wavelength, stokes->nw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(stokes->noise_mode == NOISE_FROM_INTENSITY){
        MPI_Bcast(stokes->noise_level, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }else if(stokes->noise_mode == NOISE_GLOBAL){
        MPI_Bcast(stokes->noise, stokes->nw*4, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);
      }
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

int READ_PROFILE(STRUCT_PROFILE_IO *input, STRUCT_STK *stokes, 
    STRUCT_SUBSET *subset){

    /*######################################################################
      Purpose:
        Read the observed Stokes profile at the current pixel.
      Input parameters:
        input, the input configuration.
        stokes, profile dimensions.
        subset, current pixel coordinates.
      Output parameters:
        input, with profile_buffer filled from the observation file.
        subset, advanced to the next pixel with its processed count updated.
      Return:
        CFITSIO status code.
    ######################################################################*/

    int status = 0;

    if(input->profile_format == PROFILE_FORMAT_DAT){
      subset->processed = 1;
      subset->coord[0] = 0;
      subset->coord[1] = 1;
      return 0;
    }

    // get the first pixel
    fits_first_pixel[2] = subset->coord[0]+1;
    fits_first_pixel[3] = subset->coord[1]+1;
    const long noise_x = fits_first_pixel[2];
    const long noise_y = fits_first_pixel[3];

    subset->processed++;
    if(subset->coord[0] < input->sol_box[0][1]){
      subset->coord[0]++;
    }else{
      subset->coord[0] = input->sol_box[0][0];
      subset->coord[1]++;
    }

    fits_open_file(&obs_file, input->data_path, 
        READONLY, &status);
    if(status) IO_ERROR("error in opening the data file: " 
        "status = %d \n", status); 

    fits_read_pix(obs_file, TDOUBLE, fits_first_pixel, stokes->nw*4, 
        NULL, stokes->profile, NULL, &status);

    if(status) IO_ERROR("error in reading profiles: " 
        "status = %d \n", status);  

    if(stokes->noise_mode == NOISE_PER_PIXEL){
      long noise_first[4] = {1, 1, noise_x, noise_y};
      fits_movabs_hdu(obs_file, 2, NULL, &status);
      fits_read_pix(obs_file, TDOUBLE, noise_first, stokes->nw*4, NULL, 
          stokes->noise, NULL, &status);
      if(status) IO_ERROR("error reading the per-pixel noise HDU.\n");
    }

    status = CLOSE_FITS_FILE(&obs_file);
    if(status) IO_ERROR("error closing the profile data file: " 
        "status = %d \n", status);  

    return status;
}

/*--------------------------------------------------------------------------------*/

int PIXEL_ADVANCE(STRUCT_PROFILE_IO *input, STRUCT_SUBSET *subset){

    /*######################################################################
      Purpose:
        Advance by one profile without reading it.
      Input parameters:
        input, the input configuration.
        subset, a structure storing the pixels to read.
      Output parameters:
        subset, a structure storing the pixels to read.
      Return:
        1 when work remains; 0 after the final profile.
    ######################################################################*/

    subset->processed++;
    if(subset->processed >= input->counts){
      subset->coord[0] = input->sol_box[0][0];
      subset->coord[1] = input->sol_box[1][0]+1;
      return 0;
    } 
    if(subset->coord[0] < input->sol_box[0][1]) subset->coord[0]++;
    else{
      subset->coord[0] = input->sol_box[0][0];
      subset->coord[1]++;
    }
        
    return 1;
}

/*--------------------------------------------------------------------------------*/

int CACHE_INIT(STRUCT_PROFILE_IO *input, STRUCT_MPI *mpi, 
    const STRUCT_PARA *params, const STRUCT_DREAM *dream, 
    const STRUCT_STK *stokes){

    /*######################################################################
      Purpose:
        Initialize or validate the cache and result files.
      Input parameters:
        input, the input configuration.
        mpi, structure with mpi configuration.
        stokes, structure with Stokes profiles.
      Output parameters:
        input, with cache state allocated and loaded when available.
      Return:
        0 on success; fatal I/O errors terminate through IO_ERROR.
    ######################################################################*/

    if(mpi->world_rank == 0){
      input->config_hash = CONFIG_HASH(input, params, dream, stokes);
    }
    MPI_Bcast(&input->config_hash, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    memcpy(input->cache_header.magic, CACHE_MAGIC,
        sizeof(input->cache_header.magic));
    input->cache_header.x_begin = input->sol_box[0][0];
    input->cache_header.x_end = input->sol_box[0][1];
    input->cache_header.y_begin = input->sol_box[1][0];
    input->cache_header.y_end = input->sol_box[1][1];
    input->cache_header.config_hash = input->config_hash;
    // size of the cache matrix.
    int nxcache = input->cache_header.nx;
    int nycache = input->cache_header.ny;
    int ncache = input->cache_header.ncache;

    input->counts = nxcache*nycache;

    bool has_cache = false;
    input->cache_reused = false;
    STRUCT_CACHE *header = &(input->cache_header);
    char filename[Max_Line_Length+2];

    // CFITSIO status value MUST be initialized to zero!
    int status = 0, npar = params->nmodel+1, naxis;
    long naxes[4];

    if(mpi->world_rank==0){
      input->cache = calloc((size_t)ncache, sizeof(*input->cache));
      if(!input->cache) IO_ERROR("cannot allocate the inversion cache.\n");
      if(!input->use_cache){
        has_cache = false;
      }else{
        cache_file = fopen(input->cache_path, "rb");

        if(cache_file != NULL){
          STRUCT_CACHE disk_header = {0};
          has_cache = fread(&disk_header, sizeof(disk_header), 1, 
              cache_file) == 1;
          has_cache = has_cache
              && memcmp(disk_header.magic, CACHE_MAGIC,
                  sizeof(disk_header.magic)) == 0;

          has_cache = has_cache && nxcache == disk_header.nx 
              && nycache == disk_header.ny && ncache == disk_header.ncache 
              && disk_header.x_begin == header->x_begin 
              && disk_header.x_end == header->x_end 
              && disk_header.y_begin == header->y_begin 
              && disk_header.y_end == header->y_end 
              && disk_header.config_hash == header->config_hash;
          if(has_cache){
            
            has_cache = fread(input->cache, sizeof(int), 
                (size_t)ncache, cache_file) == (size_t)ncache;
            if(has_cache){
              input->cache_reused = true;
              IO_VERBOSE(2, "\n *** cache file read. ***\n");
            }

          }
          if(!has_cache){
            IO_VERBOSE(2, "\n ** %s is not a correct cache file.", 
                input->cache_path);     
          }
          if(CLOSE_STD_FILE(&cache_file) != 0){
            IO_ERROR("cannot close cache file after reading.\n");
          }
        }
      }

      if(!has_cache){
        memcpy(header->magic, CACHE_MAGIC, sizeof(header->magic));
        header->nx = nxcache;
        header->ny = nycache;
        header->ncache = ncache;
        memset(input->cache, 0, (size_t)ncache*sizeof(*input->cache));
        if(input->use_cache){
          IO_VERBOSE(2, "\n  *** no correct cache file found. ***");     
          IO_VERBOSE(2, "\n   ** creating a cache file.");  

          cache_file = fopen(input->cache_path, "wb");
          if(!cache_file) IO_ERROR("cannot create cache file.\n");
          if(fwrite(header, sizeof(*header), 1, cache_file) != 1 
              || fwrite(input->cache, sizeof(int), (size_t)ncache, 
              cache_file) != (size_t)ncache){
            IO_ERROR("cannot initialize cache file.\n");
          }
          if(CLOSE_STD_FILE(&cache_file) != 0){
            IO_ERROR("cannot flush and close the initialized cache file.\n");
          }
        }

        IO_VERBOSE(2, "\n   ** creating a fits file for the result.");

        snprintf(filename, sizeof(filename), "!%s", 
            input->result_path);
        fits_create_file(&res_file, filename, &status);
        if(status) IO_ERROR("error creating the result file: " 
            "status = %d \n", status);   

        naxis = 3;
        naxes[0] = npar;
        naxes[1] = nxcache;
        naxes[2] = nycache;

        fits_create_img(res_file, DOUBLE_IMG, naxis, naxes, &status);
        if(status) IO_ERROR("error creating the image HDU: " 
            "status = %d \n", status);    

        fits_update_key(res_file, TSTRING, "EXTNAME", "BEST_FIT", 
            "model parameters and log likelihood", &status);
        int nmodel = params->nmodel;
        int nlines = params->nlines;
        int nregions = params->nregions;
        fits_update_key(res_file, TINT, "NMODEL", &nmodel, 
            "number of model parameters", &status);
        fits_update_key(res_file, TINT, "NLINES", &nlines, 
            "number of spectral lines", &status);
        fits_update_key(res_file, TINT, "NREGION", &nregions, 
            "number of spectral regions", &status);
        char config_hash_text[17];
        snprintf(config_hash_text, sizeof(config_hash_text), "%016llx",
            (unsigned long long)input->config_hash);
        fits_update_key(res_file, TSTRING, "CFGHASH", config_hash_text,
            "inversion configuration hash", &status);
        for(int imodel=0; imodel<params->nmodel; imodel++){
          char key[9], name[64];
          MODEL_PARAMETER_KEY(imodel, key);
          MODEL_PARAMETER_NAME(params, imodel, name, sizeof(name));
          fits_update_key(res_file, TSTRING, key, name, 
              "model parameter name", &status);
        }

        long error_axes[4] = {2, params->nmodel, nxcache, nycache};
        fits_create_img(res_file, DOUBLE_IMG, 4, error_axes, &status);
        fits_update_key(res_file, TSTRING, "EXTNAME", "PARAMETER_ERROR",
            "lower and upper 68.27 percent errors", &status);
        fits_update_key(res_file, TSTRING, "ERRTYPE", "LOWER_UPPER",
            "axis 1 stores non-negative lower and upper errors", &status);

        int fits_status = status;
        int close_status = CLOSE_FITS_FILE(&res_file);
        if(fits_status) IO_ERROR("error creating result FITS metadata or "
            "parameter-error HDU: status = %d.\n", fits_status);
        if(close_status) IO_ERROR("error in closing the result file.\n");

      }else{

        IO_VERBOSE(2, "\n   ** opening the corresponding "
            "result files.");
        
        fits_open_file(&res_file, input->result_path, 
            READWRITE, &status);
 
        if(status) IO_ERROR("error in opening the result file: " 
            "status = %d \n",status);

        // get the dimension of the data hdu
        fits_get_img_dim(res_file, &naxis, &status);
        if(status) IO_ERROR("error in getting dimension of the"
            " result file: status = %d \n", status);   
        if(naxis!=3) IO_ERROR( "error in the dimension of the"
            " result file: status = %d \n", status);   
        // get the size of each dimension 
        fits_get_img_size(res_file, naxis, naxes, &status);
        if(status) IO_ERROR("error in getting size of the"
            " result file: status = %d \n", status);   

        if(naxes[0]!=npar || naxes[1]!=nxcache || naxes[2]!=nycache){
          IO_ERROR("wrong size of the result file: status = %d \n", 
              status);   
        } 

        char extname[FLEN_VALUE] = "";
        char disk_hash[FLEN_VALUE] = "";
        char expected_hash[17];
        int disk_nmodel = 0, disk_nlines = 0, disk_nregions = 0;
        snprintf(expected_hash, sizeof(expected_hash), "%016llx",
            (unsigned long long)input->config_hash);
        fits_read_key(res_file, TSTRING, "EXTNAME", extname, NULL, &status);
        fits_read_key(res_file, TINT, "NMODEL", &disk_nmodel, NULL, &status);
        fits_read_key(res_file, TINT, "NLINES", &disk_nlines, NULL, &status);
        fits_read_key(res_file, TINT, "NREGION", &disk_nregions, NULL, &status);
        fits_read_key(res_file, TSTRING, "CFGHASH", disk_hash, NULL, &status);
        if(status || strcmp(extname, "BEST_FIT") != 0
            || disk_nmodel != params->nmodel || disk_nlines != params->nlines
            || disk_nregions != params->nregions
            || strcmp(disk_hash, expected_hash) != 0){
          IO_ERROR("result FITS metadata does not match this inversion.\n");
        }
        for(int imodel=0; imodel<params->nmodel; imodel++){
          char key[9], disk_name[FLEN_VALUE] = "", expected_name[64];
          MODEL_PARAMETER_KEY(imodel, key);
          MODEL_PARAMETER_NAME(params, imodel, expected_name,
              sizeof(expected_name));
          fits_read_key(res_file, TSTRING, key, disk_name, NULL, &status);
          if(status || strcmp(disk_name, expected_name) != 0){
            IO_ERROR("result FITS parameter metadata is inconsistent.\n");
          }
        }

        int nhdus = 0;
        fits_get_num_hdus(res_file, &nhdus, &status);
        if(status || nhdus < 2) IO_ERROR(
            "result file does not contain the parameter-error HDU.\n");
        fits_movabs_hdu(res_file, 2, NULL, &status);
        fits_get_img_dim(res_file, &naxis, &status);
        fits_get_img_size(res_file, 4, naxes, &status);
        extname[0] = '\0';
        fits_read_key(res_file, TSTRING, "EXTNAME", extname, NULL, &status);
        if(status || naxis != 4 || naxes[0] != 2 
            || naxes[1] != params->nmodel 
            || naxes[2] != nxcache || naxes[3] != nycache
            || strcmp(extname, "PARAMETER_ERROR") != 0){
          IO_ERROR("wrong size or type of the parameter-error HDU.\n");
        }
        status = CLOSE_FITS_FILE(&res_file);
        if(status) IO_ERROR("error in closing the result file.\n");

      }

      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);

    }else{

      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(status) ABORTED();

    }

    int cache_reused = input->cache_reused ? 1 : 0;
    MPI_Bcast(&cache_reused, 1, MPI_INT, 0, MPI_COMM_WORLD);
    input->cache_reused = cache_reused != 0;

    return 0;
}

/*--------------------------------------------------------------------------------*/

int WRITE_RESULT(STRUCT_PROFILE_IO *input, STRUCT_SUBSET *subset,
    bool mark_complete){

    /*######################################################################
      Purpose:
        Write inversion results and optionally mark the subset as cached.
      Input parameters:
        input, the input configuration.
        subset, result pixel coordinates.
        mark_complete, whether to mark the cache entry as successfully
          completed.
      Output parameters:
        input, whose result buffer is written to disk.
      Return:
        0 on success.
    ######################################################################*/

    if(!input || !subset || !input->result_buffer || !input->error_buffer){
      IO_ERROR("invalid result-writing state.\n");
    }
    if(subset->coord[0] < input->sol_box[0][0]
        || subset->coord[0] > input->sol_box[0][1]
        || subset->coord[1] < input->sol_box[1][0]
        || subset->coord[1] > input->sol_box[1][1]){
      IO_ERROR("result pixel lies outside sol_box.\n");
    }

    int status = 0;
    long fpixel[4] = {1, 1, 1, 1};

    fpixel[0] = 1;
    fpixel[1] = subset->coord[0]+1-input->sol_box[0][0];
    fpixel[2] = subset->coord[1]+1-input->sol_box[1][0];

    fits_open_file(&res_file, input->result_path, 
        READWRITE, &status);
    if(status) IO_ERROR("error in opening the result file: " 
        "status = %d \n",status);
    fits_write_pix(res_file, TDOUBLE, fpixel, input->result_nparams+1, 
          input->result_buffer, &status);
    if(status) IO_ERROR("error in writing the result file: " 
        "status = %d \n", status);  

    long error_pixel[4] = {1, 1, fpixel[1], fpixel[2]};
    fits_movabs_hdu(res_file, 2, NULL, &status);
    fits_write_pix(res_file, TDOUBLE, error_pixel,
        input->result_nparams*2, input->error_buffer, &status);
    if(status) IO_ERROR("error writing the parameter errors: "
        "status = %d \n", status);

    status = CLOSE_FITS_FILE(&res_file);
    if(status) IO_ERROR("error in closing the result file.\n");

    if(input->use_cache && mark_complete){
      long offset = (long)sizeof(input->cache_header) 
          +(long)((subset->coord[1]-input->sol_box[1][0]) 
          *input->cache_header.nx 
          +(subset->coord[0]-input->sol_box[0][0]))*(long)sizeof(int);
      int one = 1;

      cache_file = fopen(input->cache_path, "r+b");
      if(!cache_file || fseek(cache_file, offset, SEEK_SET) != 0){
        IO_ERROR("cannot update cache file.\n");
      }

      if(fwrite(&one, sizeof(int), 1, cache_file) != 1){
        IO_ERROR("cannot write cache entry.\n");
      }
    
      if(CLOSE_STD_FILE(&cache_file) != 0){
        IO_ERROR("cannot flush and close the updated cache file.\n");
      }
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

int CLOSE_FILES(void){

    /*######################################################################
      Purpose:
        Close all files owned by the I/O and logging modules.
      Input parameters:
        None.
      Return:
        Number of cache, FITS, and log files that failed to close.
    ######################################################################*/

    int status = 0;
    if(CLOSE_STD_FILE(&cache_file) != 0){
      LOG_ERROR(ERR_LVL_WARNING, "CLOSE_FILES", "Error closing cache file.");
      status++;
    }
    if(CLOSE_FITS_FILE(&obs_file)!=0) status++;

    if(CLOSE_FITS_FILE(&res_file)!=0) status++;

    if(LOG_FINALIZE() != 0){
      LOG_WRITE("-WARNING in routine CLOSE_FILES: Error closing log file.",
          true, true);
      status++;
    }

    return status;
}

/*--------------------------------------------------------------------------------*/
