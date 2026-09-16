
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input_reader.h"
#include "logger.h"
#include "me_solver.h"
#include "dream.h"

/*--------------------------------------------------------------------------------*/

#define IS_YES(str) (strcasecmp(str,"YES")==0)
#define IS_NO(str)  (strcasecmp(str,"NO")==0)

/*--------------------------------------------------------------------------------*/

static int PARAMETER_LAYOUT_INIT(STRUCT_PARA *params){
  
    /*######################################################################
      Purpose:
        Build the dynamic model-parameter layout for all lines and regions.
      Input parameters:
        params, global defaults and per-line/per-region overrides.
      Output parameters:
        params, populated parameter kinds, indices, bounds, and controls.
      Return:
        0 on success; -1 when the model exceeds MAX_MODEL_PARAMS.
    ######################################################################*/

    double default_limits[5][2];
    bool default_inv[5];
    double default_value[5];
    for(int i=0; i<5; i++){
      memcpy(default_limits[i], params->limits[4+i], 
          sizeof(default_limits[i]));
      default_inv[i] = params->inv[4+i];
      default_value[i] = params->value_const[4+i];
    }

    params->kind[0] = PARAM_B0;
    params->kind[1] = PARAM_B1;
    params->kind[2] = PARAM_B2;
    params->kind[3] = PARAM_VLOS;
    int index = 4;
    for(int iline=0; iline<params->nlines; iline++){
      STRUCT_MELINE *line = &params->lines[iline];
      line->dopp_index = index++;
      line->damp_index = index++;
      line->eta_index = index++;
      int indices[3] = {line->dopp_index, line->damp_index, line->eta_index};
      MODEL_PARAMETER_KIND kinds[3] = {
          PARAM_DOPPLER, PARAM_DAMPING, PARAM_ETA};
      for(int i=0; i<3; i++){
        int model_index = indices[i];
        memcpy(params->limits[model_index], default_limits[i], 
            sizeof(params->limits[model_index]));
        params->inv[model_index] = default_inv[i];
        params->value_const[model_index] = default_value[i];
        params->kind[model_index] = kinds[i];
        if(line->custom_bounds[i]){
          memcpy(params->limits[model_index], line->custom_limits[i], 
              sizeof(params->limits[model_index]));
        }
        if(line->custom_inversion[i]){
          params->inv[model_index] = line->custom_inv[i];
          params->value_const[model_index] = line->custom_value[i];
        }
      }
    }
    for(int iregion=0; iregion<params->nregions; iregion++){
      STRUCT_REGION *region = &params->regions[iregion];
      region->continuum_index = index++;
      region->beta_index = index++;
      int indices[2] = {region->continuum_index, region->beta_index};
      MODEL_PARAMETER_KIND kinds[2] = {PARAM_CONTINUUM, PARAM_BETA};
      for(int i=0; i<2; i++){
        int model_index = indices[i];
        memcpy(params->limits[model_index], default_limits[3+i], 
            sizeof(params->limits[model_index]));
        params->inv[model_index] = default_inv[3+i];
        params->value_const[model_index] = default_value[3+i];
        params->kind[model_index] = kinds[i];
        if(region->custom_bounds[i]){
          memcpy(params->limits[model_index], region->custom_limits[i], 
              sizeof(params->limits[model_index]));
        }
        if(region->custom_inversion[i]){
          params->inv[model_index] = region->custom_inv[i];
          params->value_const[model_index] = region->custom_value[i];
        }
      }
    }
    if(index > MAX_MODEL_PARAMS) return -1;
    params->nmodel = index;
    params->npar = 0;
    for(int i=0; i<params->nmodel; i++){
      if(params->inv[i]) params->npar++;
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int MODEL_LAYOUT_VALIDATE(const STRUCT_PARA *params){

    /*######################################################################
      Purpose:
        Validate numerical and physical constraints in the model layout.
      Input parameters:
        params, completed model-parameter layout.
      Return:
        0 when all parameters are valid; -1 otherwise.
    ######################################################################*/

    for(int index=0; index<params->nmodel; index++){
      double lower = params->limits[index][0];
      double upper = params->limits[index][1];
      double fixed = params->value_const[index];
      if(!isfinite(lower) || !isfinite(upper) || lower > upper
          || !isfinite(upper-lower)
          || (!params->inv[index] && !isfinite(fixed))) return -1;
      if(params->inv[index] && params->kind[index] != PARAM_CONTINUUM 
          && lower == upper) return -1;
      if(!params->inv[index] && params->kind[index] != PARAM_CONTINUUM 
          && (fixed < lower || fixed > upper)) return -1;
      double physical_lower = params->inv[index] ? lower : fixed;
      switch(params->kind[index]){
        case PARAM_DOPPLER:
          if(physical_lower <= 0.0) return -1;
          break;
        case PARAM_DAMPING:
        case PARAM_ETA:
          if(physical_lower < 0.0) return -1;
          break;
        case PARAM_BETA:
          if(physical_lower < 0.0 
              || (params->inv[index] ? upper : fixed) > 1.0) return -1;
          break;
        case PARAM_CONTINUUM:
          if(!params->inv[index] && fixed <= 0.0) return -1;
          break;
        default:
          break;
      }
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int SPECTRAL_REGIONS_INIT(STRUCT_PARA *params, STRUCT_STK *stokes){

    /*######################################################################
      Purpose:
        Split the sorted wavelength grid at gaps larger than region_gap.
      Input parameters:
        params, storage for spectral-region metadata.
        stokes, wavelength grid and configured separation threshold.
      Output parameters:
        params, populated wavelength regions and index ranges.
      Return:
        0 on success; -1 for an invalid grid; -2 for too many regions.
    ######################################################################*/

    if(stokes->nw < 2 || !isfinite(stokes->region_gap) 
        || stokes->region_gap <= 0.0) return -1;

    int region_count = 1;
    for(int iw=1; iw<stokes->nw; iw++){
      if(stokes->wavelength[iw]-stokes->wavelength[iw-1] 
          > stokes->region_gap) region_count++;
    }
    params->nregions = region_count;
    if(region_count > NUM_REGIONS) return -2;

    int region_id = 0;
    int begin = 0;
    for(int iw=1; iw<stokes->nw; iw++){
      if(stokes->wavelength[iw]-stokes->wavelength[iw-1] 
          <= stokes->region_gap) 
          continue;
      STRUCT_REGION *region = &params->regions[region_id++];
      region->iw_begin = begin;
      region->iw_end = iw-1;
      region->wavelength_min = stokes->wavelength[begin];
      region->wavelength_max = stokes->wavelength[iw-1];
      begin = iw;
    }
    STRUCT_REGION *region = &params->regions[region_id];
    region->iw_begin = begin;
    region->iw_end = stokes->nw-1;
    region->wavelength_min = stokes->wavelength[begin];
    region->wavelength_max = stokes->wavelength[stokes->nw-1];
    params->nregions = region_id+1;
    return 0;
}

/*--------------------------------------------------------------------------------*/

int Model_Layout_Init(STRUCT_PARA *params, STRUCT_STK *stokes, STRUCT_MPI *mpi){

    /*######################################################################
      Purpose:
        Construct and validate the complete multi-line model layout.
      Input parameters:
        params, parsed parameter defaults and spectral-line definitions.
        stokes, observed wavelength grid and region controls.
        mpi, MPI rank information used for reporting errors.
      Output parameters:
        params, finalized model, line, and region parameter mappings.
      Return:
        0 on success; -1 when region construction or validation fails.
    ######################################################################*/

    int region_status = SPECTRAL_REGIONS_INIT(params, stokes);
    if(region_status != 0){
      if(mpi->world_rank == 0){
        if(region_status == -2){
          snprintf(message_buffer, sizeof(message_buffer), 
              "region_gap=%g Angstrom produces %d regions; the maximum is " 
              "%d. Increase region_gap or NUM_REGIONS.\n", 
              stokes->region_gap, params->nregions, NUM_REGIONS);
          LOG_ERROR(ERR_LVL_ERROR, "Model_Layout_Init", message_buffer);
        }else{
          LOG_ERROR(ERR_LVL_ERROR, "Model_Layout_Init", 
              "failed to split the wavelength grid into spectral regions.\n");
        }
      }
      return -1;
    }
    for(int iregion=params->nregions; iregion<NUM_REGIONS; iregion++){
      for(int itype=0; itype<2; itype++){
        if(params->regions[iregion].custom_bounds[itype] 
            || params->regions[iregion].custom_inversion[itype]){
          if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, "Model_Layout_Init", 
              "indexed region parameter refers to an undefined region.\n");
          return -1;
        }
      }
    }
    if(PARAMETER_LAYOUT_INIT(params) != 0 
        || MODEL_LAYOUT_VALIDATE(params) != 0){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, "Model_Layout_Init", 
          "failed to build spectral regions or dynamic model parameters.\n");
      return -1;
    }
    return 0;
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_LONG_TOKEN(const char **cursor, long *value){

    /*######################################################################
      Purpose:
        Parse one base-10 long integer and advance the input cursor.
      Input parameters:
        cursor, address of the current position in the input string.
      Output parameters:
        cursor, advanced to the first unparsed character.
        value, parsed integer value.
      Return:
        true on success; false for a missing or out-of-range value.
    ######################################################################*/

    char *end = NULL;
    errno = 0;
    long parsed = strtol(*cursor, &end, 10);
    if(end == *cursor || errno == ERANGE) return false;
    *cursor = end;
    *value = parsed;
    return true;
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_DOUBLE_TOKEN(const char **cursor, double *value){

    /*######################################################################
      Purpose:
        Parse one finite double value and advance the input cursor.
      Input parameters:
        cursor, address of the current position in the input string.
      Output parameters:
        cursor, advanced to the first unparsed character.
        value, parsed finite value.
      Return:
        true on success; false for an invalid or out-of-range value.
    ######################################################################*/

    char *end = NULL;
    errno = 0;
    double parsed = strtod(*cursor, &end);
    if(end == *cursor || errno == ERANGE || !isfinite(parsed)) return false;
    *cursor = end;
    *value = parsed;
    return true;
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_COMMA(const char **cursor){

    /*######################################################################
      Purpose:
        Consume optional whitespace followed by one comma separator.
      Input parameters:
        cursor, address of the current position in the input string.
      Output parameters:
        cursor, advanced past the comma when present.
      Return:
        true when a comma is consumed; false otherwise.
    ######################################################################*/

    while(isspace((unsigned char)**cursor)) (*cursor)++;
    if(**cursor != ',') return false;
    (*cursor)++;
    return true;
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_INT(const char *text, int *value){

    /*######################################################################
      Purpose:
        Parse a complete integer value with an optional trailing comment.
      Input parameters:
        text, input text to parse.
      Output parameters:
        value, parsed integer.
      Return:
        true on success; false when the text is not a valid int.
    ######################################################################*/

    const char *cursor = text;
    long parsed;
    if(!PARSE_LONG_TOKEN(&cursor, &parsed) || parsed < INT_MIN 
        || parsed > INT_MAX) return false;
    while(isspace((unsigned char)*cursor)) cursor++;
    if(*cursor != '\0' && *cursor != '#') return false;
    *value = (int)parsed;
    return true;
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_DOUBLE(const char *text, double *value){

    /*######################################################################
      Purpose:
        Parse a complete finite double with an optional trailing comment.
      Input parameters:
        text, input text to parse.
      Output parameters:
        value, parsed double value.
      Return:
        true on success; false when the complete value is invalid.
    ######################################################################*/

    const char *cursor = text;
    if(!PARSE_DOUBLE_TOKEN(&cursor, value)) return false;
    while(isspace((unsigned char)*cursor)) cursor++;
    return *cursor == '\0' || *cursor == '#';
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_TWO_DOUBLES(const char *text, double *first, double *second){

    /*######################################################################
      Purpose:
        Parse two finite comma-separated double values.
      Input parameters:
        text, input text to parse.
      Output parameters:
        first, first parsed value.
        second, second parsed value.
      Return:
        true on success; false for invalid syntax or trailing characters.
    ######################################################################*/

    const char *cursor = text;
    if(!PARSE_DOUBLE_TOKEN(&cursor, first) || !PARSE_COMMA(&cursor) 
        || !PARSE_DOUBLE_TOKEN(&cursor, second)) return false;
    while(isspace((unsigned char)*cursor)) cursor++;
    return *cursor == '\0' || *cursor == '#';
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_FOUR_DOUBLES(const char *text, double values[4]){

    /*######################################################################
      Purpose:
        Parse four finite comma-separated double values.
      Input parameters:
        text, input text to parse.
      Output parameters:
        values, the four parsed values.
      Return:
        true on success; false for invalid syntax or trailing characters.
    ######################################################################*/

    const char *cursor = text;
    for(int i=0; i<4; i++){
      if(!PARSE_DOUBLE_TOKEN(&cursor, &values[i])) return false;
      if(i < 3 && !PARSE_COMMA(&cursor)) return false;
    }
    while(isspace((unsigned char)*cursor)) cursor++;
    return *cursor == '\0' || *cursor == '#';
}

/*--------------------------------------------------------------------------------*/

static bool PARSE_FOUR_INTS(const char *text, int values[4]){

    /*######################################################################
      Purpose:
        Parse four comma-separated values in the range of int.
      Input parameters:
        text, input text to parse.
      Output parameters:
        values, the four parsed integers.
      Return:
        true on success; false for invalid syntax or trailing characters.
    ######################################################################*/

    const char *cursor = text;
    for(int i=0; i<4; i++){
      long parsed;
      if(!PARSE_LONG_TOKEN(&cursor, &parsed) || parsed < INT_MIN 
          || parsed > INT_MAX) return false;
      values[i] = (int)parsed;
      if(i < 3 && !PARSE_COMMA(&cursor)) return false;
    }
    while(isspace((unsigned char)*cursor)) cursor++;
    return *cursor == '\0' || *cursor == '#';
}

/*--------------------------------------------------------------------------------*/

static bool INDEXED_KEY(const char *key, const char *prefix, int *index){

    /*######################################################################
      Purpose:
        Match a keyword prefix followed by a complete integer index.
      Input parameters:
        key, keyword to examine.
        prefix, required prefix text.
      Output parameters:
        index, parsed suffix index when the keyword matches.
      Return:
        true for a valid indexed keyword; false otherwise.
    ######################################################################*/

    size_t prefix_length = strlen(prefix);
    if(strncmp(key, prefix, prefix_length) != 0) return false;
    const char *cursor = key+prefix_length;
    long parsed;
    if(!PARSE_LONG_TOKEN(&cursor, &parsed) || *cursor != '\0' 
        || parsed < INT_MIN || parsed > INT_MAX) return false;
    *index = (int)parsed;
    return true;
}

/*--------------------------------------------------------------------------------*/

static int INVERSION_VALUE(const char *text, bool *invert, double *value){

    /*######################################################################
      Purpose:
        Parse an inversion switch or a fixed model-parameter value.
      Input parameters:
        text, YES, NO,value, or a standalone fixed numeric value.
      Output parameters:
        invert, whether the parameter participates in the inversion.
        value, fixed value when inversion is disabled.
      Return:
        0 on success; -1 for invalid input syntax.
    ######################################################################*/

    const char *comma = strchr(text, ',');
    if(!comma){
      char mode[Key_Length];
      size_t length = strlen(text);
      if(length >= sizeof(mode)) return -1;
      memcpy(mode, text, length+1);
      char *comment = strchr(mode, '#');
      if(comment) *comment = '\0';
      STR_TRIM(mode);
      STR_TOUPPER(mode);
      if(strcmp(mode, "YES") == 0){
        *invert = true;
        *value = 0.0;
        return 0;
      }
      double fixed;
      if(!PARSE_DOUBLE(text, &fixed)) return -1;
      *invert = false;
      *value = fixed;
      return 0;
    }

    char mode[Key_Length];
    size_t mode_length = (size_t)(comma-text);
    if(mode_length >= sizeof(mode)) return -1;
    memcpy(mode, text, mode_length);
    mode[mode_length] = '\0';
    STR_TRIM(mode);
    STR_TOUPPER(mode);
    double fixed;
    if(strcmp(mode, "NO") == 0 && PARSE_DOUBLE(comma+1, &fixed)){
      *invert = false;
      *value = fixed;
      return 0;
    }
    return -1;
}

/*--------------------------------------------------------------------------------*/

static int Get_Keys(STRUCT_KEYS keywords[], STRUCT_PROFILE_IO *input, 
    STRUCT_PARA *params, STRUCT_DREAM *dream, STRUCT_STK *stokes, 
    STRUCT_MPI *mpi){
  
    /*######################################################################
      Purpose:
        Convert the input keywords into the runtime configuration.
      Input parameters:
        keywords, the parsed input keywords.
        mpi, the MPI configuration.
      Output parameters:
        input, the populated runtime configuration.
     ######################################################################*/
    
    const char *rname = "Get_Keys";
    const double pi = 3.14159265358979323846;
    
    #define INPUT_VERBOSE(...)                                              \
      do{                                                                   \
        if(mpi->verbose_level >= 3 && mpi->world_rank == 0){                \
          snprintf(message_buffer, sizeof(message_buffer), __VA_ARGS__);    \
          LOG_WRITE(message_buffer, true, true);                            \
        }                                                                   \
      }while(0)

    #define INPUT_FLAG(indx, var, desc)                                     \
      do{                                                                   \
        STR_TRIM(keywords[indx].line);                                      \
        if(!IS_YES(keywords[indx].line) && !IS_NO(keywords[indx].line)){    \
          if(mpi->world_rank == 0){                                         \
            snprintf(message_buffer, sizeof(message_buffer),                \
                "%s must be YES or NO.\n", desc);                           \
            LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);                \
          }                                                                 \
          return -1;                                                        \
        }                                                                   \
        var = IS_YES(keywords[indx].line);                                  \
        INPUT_VERBOSE("\n %s: %s", desc, var ? "Yes" : "No");               \
      }while(0)


    #define PAR_BOUNDS(indx, par_indx, desc, unit)                          \
      do{                                                                   \
        if(!PARSE_TWO_DOUBLES(keywords[indx].line,                          \
            &params->limits[par_indx][0], &params->limits[par_indx][1])){   \
          if(mpi->world_rank==0){                                           \
          snprintf(message_buffer, sizeof(message_buffer),                  \
              "Invalid lower and upper bounds for %s.\n", desc);            \
          LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);                  \
          }                                                                 \
          return -1;                                                        \
        }                                                                   \
        INPUT_VERBOSE("\n Bounds on %s: from %e to %e %s ", desc,           \
            params->limits[par_indx][0], params->limits[par_indx][1], unit);\
      }while(0)

    #define PAR_INVERSION(indx, par_indx, desc)                             \
      do{                                                                   \
        if(INVERSION_VALUE(keywords[indx].line, &params->inv[par_indx],     \
            &params->value_const[par_indx]) != 0){                          \
          snprintf(message_buffer, sizeof(message_buffer),                  \
              "Invalid inversion control for %s. "                          \
              "Use YES, NO,value, or a fixed value.\n", desc);              \
          if(mpi->world_rank == 0){                                         \
            LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);                \
          }                                                                 \
          return -1;                                                        \
        }                                                                   \
        if(params->inv[par_indx]){                                          \
          INPUT_VERBOSE("\n Invert %s: Yes", desc);                         \
        }else{                                                              \
          INPUT_VERBOSE("\n Invert %s: No, fixed value: %e", desc,          \
              params->value_const[par_indx]);                               \
        }                                                                   \
      }while(0)

    #define INPUT_INT_RANGE(indx, var, min_value, max_value)                \
      do{                                                                   \
        if(!PARSE_INT(keywords[indx].line, &(var))                          \
            || (var) < (min_value) || (var) > (max_value)){                 \
          if(mpi->world_rank == 0){                                         \
            snprintf(message_buffer, sizeof(message_buffer),                \
                "%s must be an integer from %d to %d.\n",                   \
                keywords[indx].keyword, (min_value), (max_value));          \
            LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);                \
          }                                                                 \
          return -1;                                                        \
        }                                                                   \
      }while(0)
    
    INPUT_INT_RANGE(KEY_VERBOSE, mpi->verbose_level, 0, 4);
    if(mpi->verbose_level > 0){
      if(mpi->world_rank == 0){
        if(LOG_INIT(input->log_path) != 0){
          LOG_ERROR(ERR_LVL_ERROR, rname, "Cannot open the log file.");
        }
      }else if(mpi->verbose_level >= 4){
        char rank_log[Max_Line_Length+32];
        snprintf(rank_log, sizeof(rank_log), "%s.rank_%05d",
            input->log_path, mpi->world_rank);
        if(LOG_INIT(rank_log) != 0){
          LOG_ERROR(ERR_LVL_ERROR, rname, "Cannot open the rank log file.");
        }
      }
    }
    INPUT_VERBOSE("\n Verbose level: %d", mpi->verbose_level);

    STR_TRIM(keywords[KEY_PROFILE_FORMAT].line);
    STR_TOUPPER(keywords[KEY_PROFILE_FORMAT].line);
    if(strcmp(keywords[KEY_PROFILE_FORMAT].line, "FITS") == 0){
      input->profile_format = PROFILE_FORMAT_FITS;
    }else if(strcmp(keywords[KEY_PROFILE_FORMAT].line, "DAT") == 0){
      input->profile_format = PROFILE_FORMAT_DAT;
    }else{
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "profile_format must be FITS or DAT.\n");
      return -1;
    }
    INPUT_VERBOSE("\n Profile format: %s", 
        keywords[KEY_PROFILE_FORMAT].line);

    for(int i=0; i<params->nlines; i++){
      INPUT_VERBOSE( 
          "\n Line center: %e, Effective lande factor: %e", 
          params->lines[i].wavelength0, params->lines[i].lande_factor);
    }

    STR_COPY(input->data_path, Max_Line_Length, 
      keywords[KEY_DATA_PATH].line, 
      strlen(keywords[KEY_DATA_PATH].line), true);
    STR_TRIM(input->data_path);
    INPUT_VERBOSE("\n Data file path: %s", 
        input->data_path);
    if(!FILE_EXIST(input->data_path) && mpi->world_rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname, "Data file does not exist.\n");
    }

    STR_COPY(input->wavelength_path, Max_Line_Length, 
        keywords[KEY_WAV_PATH].line, 
        strlen(keywords[KEY_WAV_PATH].line), true);
    STR_TRIM(input->wavelength_path);
    INPUT_VERBOSE("\n Wavelength file path: %s", 
        input->wavelength_path);
    if(input->profile_format == PROFILE_FORMAT_FITS 
        && !FILE_EXIST(input->wavelength_path) && mpi->world_rank==0){
      LOG_ERROR(ERR_LVL_ERROR, rname, 
          "Wavelength file does not exist.\n");
    }

    STR_COPY(input->result_path, Max_Line_Length, 
        keywords[KEY_RESULT_PATH].line, 
        strlen(keywords[KEY_RESULT_PATH].line), true);
    STR_TRIM(input->result_path);
    INPUT_VERBOSE("\n Result file path: %s", 
        input->result_path);

    STR_COPY(input->cache_path, Max_Line_Length, 
        keywords[KEY_CACHE_PATH].line, 
        strlen(keywords[KEY_CACHE_PATH].line), true);
    STR_TRIM(input->cache_path);
    INPUT_FLAG(KEY_USE_CACHE, input->use_cache, "Use inversion cache");

    STR_COPY(input->sample_path, Max_Line_Length, 
        keywords[KEY_SAMPLE_PATH].line, 
        strlen(keywords[KEY_SAMPLE_PATH].line), true);
    STR_TRIM(input->sample_path);

    INPUT_INT_RANGE(KEY_NCHAINS, dream->nchains, 11, 1000);
    INPUT_VERBOSE("\n Number of DREAM chains per island: %d",
        dream->nchains);

    INPUT_INT_RANGE(KEY_BURNIN_GENERATIONS, dream->burnin_generations, 2,
        1000000);
    INPUT_VERBOSE("\n Maximum number of burn-in generations: %d",
        dream->burnin_generations);

    INPUT_INT_RANGE(KEY_SAMPLING_GENERATIONS, dream->sampling_generations, 2,
        1000000);
    INPUT_VERBOSE("\n Number of sampling generations: %d",
        dream->sampling_generations);

    INPUT_INT_RANGE(KEY_MAX_PAIR, dream->max_pairs, 1, 1000000);
    INPUT_VERBOSE("\n Maximum number of chain pairs: %d", 
        dream->max_pairs);
    if((long long)dream->nchains <= 2LL*dream->max_pairs){
      if(mpi->world_rank == 0){
        LOG_ERROR(ERR_LVL_ERROR, rname,
            "num_chains must be greater than 2*max_pair.\n");
      }
      return -1;
    }


    INPUT_INT_RANGE(KEY_NCR, dream->ncr, 1, MAX_MODEL_PARAMS);
    INPUT_VERBOSE("\n Number of crossover probabilities: %d", 
        dream->ncr);

    INPUT_FLAG(KEY_CR_UPDATE, dream->update_crossover, 
        "Update crossover probability");


    INPUT_FLAG(KEY_PROPOSAL_NOISE_UPDATE, dream->update_proposal_noise,
        "Update proposal noise scale");


    STR_TRIM(keywords[KEY_SAMPLE_OUTPUT].line);
    STR_TOUPPER(keywords[KEY_SAMPLE_OUTPUT].line);
    if(strcmp(keywords[KEY_SAMPLE_OUTPUT].line, "NONE") == 0){
      dream->sample_output = SAMPLE_OUTPUT_NONE;
    }else if(strcmp(keywords[KEY_SAMPLE_OUTPUT].line, "MAGNETIC") == 0){
      dream->sample_output = SAMPLE_OUTPUT_MAGNETIC;
    }else if(strcmp(keywords[KEY_SAMPLE_OUTPUT].line, "ALL") == 0){
      dream->sample_output = SAMPLE_OUTPUT_ALL;
    }else{
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "sample_output must be NONE, MAGNETIC, or ALL.\n");
      return -1;
    }
    INPUT_VERBOSE("\n Sample output: %s", 
        keywords[KEY_SAMPLE_OUTPUT].line);

    INPUT_INT_RANGE(KEY_ISLAND_SIZE, mpi->requested_island_size, 1, 
        1000000);
    INPUT_VERBOSE("\n MPI processes per island: %d", 
        mpi->requested_island_size);

    if(input->profile_format == PROFILE_FORMAT_DAT){
      mpi->requested_island_size = mpi->world_size;
      INPUT_VERBOSE( 
          "\n DAT input uses one island containing all MPI ranks.");
    }

    int sol_box[4];
    if(!PARSE_FOUR_INTS(keywords[KEY_SOL_BOX].line, sol_box)){
      if(mpi->world_rank==0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "sol_box requires four integers.\n");
      return -1;
    }
    input->sol_box[0][0] = sol_box[0];
    input->sol_box[0][1] = sol_box[1];
    input->sol_box[1][0] = sol_box[2];
    input->sol_box[1][1] = sol_box[3];
    INPUT_VERBOSE( 
        "\n Solution region\n X: %d to %d \n Y: %d to %d", 
        input->sol_box[0][0], input->sol_box[0][1], 
        input->sol_box[1][0], input->sol_box[1][1]);


    INPUT_INT_RANGE(KEY_VOIGT_PRECISION, 
        stokes->faddeeva.precision_digits, 1, 8);
    INPUT_VERBOSE("\n Voigt function accuracy: %d", 
        stokes->faddeeva.precision_digits);
    STR_TRIM(keywords[KEY_NOISE_MODE].line);
    STR_TOUPPER(keywords[KEY_NOISE_MODE].line);
    if(strcmp(keywords[KEY_NOISE_MODE].line, "INTENSITY") == 0){
      stokes->noise_mode = NOISE_FROM_INTENSITY;
    }else if(strcmp(keywords[KEY_NOISE_MODE].line, "PIXEL") == 0){
      stokes->noise_mode = NOISE_PER_PIXEL;
    }else if(strcmp(keywords[KEY_NOISE_MODE].line, "GLOBAL") == 0){
      stokes->noise_mode = NOISE_GLOBAL;
    }else{
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "noise_mode must be INTENSITY, PIXEL, or GLOBAL.\n");
      return -1;
    }
    if(stokes->noise_mode == NOISE_FROM_INTENSITY){
      bool noise_values_valid = PARSE_FOUR_DOUBLES(
          keywords[KEY_NOISE_LEVEL].line, stokes->noise_level);
      if(!noise_values_valid){
        if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
            "noise_level requires four finite numeric values for I, Q, U, V.\n");
        return -1;
      }
      for(int istk=0; istk<4; istk++){
        if(stokes->noise_level[istk] <= 0.0){
          if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname,
              "noise_level values must be positive.\n");
          return -1;
        }
      }
    }
    INPUT_VERBOSE("\n Noise mode: %s", 
        keywords[KEY_NOISE_MODE].line);
    if(!PARSE_DOUBLE(keywords[KEY_REGION_GAP].line, 
        &stokes->region_gap) || stokes->region_gap <= 0.0){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "region_gap must be a finite positive wavelength interval.\n");
      return -1;
    }
    INPUT_VERBOSE("\n Spectral region gap [Angstrom]: %g", 
        stokes->region_gap);
    if(!PARSE_DOUBLE(keywords[KEY_LINE_WINGS].line, 
        &stokes->line_wing_widths) || stokes->line_wing_widths <= 0.0 
        || stokes->line_wing_widths > 10000.0){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "line_wing_widths must be in the range (0, 10000].\n");
      return -1;
    }
    INPUT_VERBOSE("\n Line-wing cutoff: %g Doppler widths", 
        stokes->line_wing_widths);
    if(!PARSE_DOUBLE(keywords[KEY_DREAM_MEMORY].line, 
        &dream->max_memory_gb) || dream->max_memory_gb <= 0.0 
        || dream->max_memory_gb > 1024.0){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "max_dream_memory_gb must be in the range (0, 1024].\n");
      return -1;
    }
    INPUT_VERBOSE("\n DREAM memory limit per rank: %g GiB", 
        dream->max_memory_gb);
    if(!PARSE_DOUBLE(keywords[KEY_SAMPLE_FILE_LIMIT].line,
        &dream->max_sample_file_gb) || dream->max_sample_file_gb <= 0.0
        || dream->max_sample_file_gb > 1048576.0){
      if(mpi->world_rank == 0){
        LOG_ERROR(ERR_LVL_ERROR, rname,
            "max_sample_file_gb must be in the range (0, 1048576].\n");
      }
      return -1;
    }
    INPUT_VERBOSE("\n Sample file size limit: %g GiB",
        dream->max_sample_file_gb);
    if(input->profile_format == PROFILE_FORMAT_DAT 
        && stokes->noise_mode != NOISE_FROM_INTENSITY){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "DAT input currently requires noise_mode = INTENSITY.\n");
      return -1;
    }

    STR_TRIM(keywords[KEY_MAGNETIC_MODE].line);
    STR_TOUPPER(keywords[KEY_MAGNETIC_MODE].line);
    if(strcmp(keywords[KEY_MAGNETIC_MODE].line, "CARTESIAN") == 0){
      params->magnetic_mode = MAGNETIC_CARTESIAN;
      PAR_BOUNDS(KEY_BZ_LIMIT, 0, "Bz", "[G]");
      PAR_BOUNDS(KEY_BX_LIMIT, 1, "Bx", "[G]");
      PAR_BOUNDS(KEY_BY_LIMIT, 2, "By", "[G]");
      PAR_INVERSION(KEY_INV_BZ, 0, "Bz");
      PAR_INVERSION(KEY_INV_BX, 1, "Bx");
      PAR_INVERSION(KEY_INV_BY, 2, "By");
    }else if(strcmp(keywords[KEY_MAGNETIC_MODE].line, "SPHERICAL") == 0){
      params->magnetic_mode = MAGNETIC_SPHERICAL;
      PAR_BOUNDS(KEY_BMOD_LIMIT, 0, "Bmod", "[G]");
      PAR_BOUNDS(KEY_BTHETA_LIMIT, 1, "ThetaB", "[Pi]");
      PAR_BOUNDS(KEY_BPHI_LIMIT, 2, "PhiB", "[Pi]");
      PAR_INVERSION(KEY_INV_BMOD, 0, "Bmod");
      PAR_INVERSION(KEY_INV_BTHETA, 1, "ThetaB");
      PAR_INVERSION(KEY_INV_BPHI, 2, "PhiB");
      params->limits[1][0] *= pi;
      params->limits[1][1] *= pi;
      params->limits[2][0] *= pi;
      params->limits[2][1] *= pi;
      params->value_const[1] *= pi;
      params->value_const[2] *= pi;
    }else{
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "magnetic_mode must be CARTESIAN or SPHERICAL.\n");
      return -1;
    }
    INPUT_VERBOSE("\n Magnetic mode: %s", 
        keywords[KEY_MAGNETIC_MODE].line);

    PAR_BOUNDS(KEY_VLOS_LIMIT, 3, "Vlos", "[km/s]");
    PAR_BOUNDS(KEY_DOPPLER_LIMIT, 4, "Doppler width", "[mA]");
    PAR_BOUNDS(KEY_DAMPING_LIMIT, 5, "Damping width", "[Dopp]");
    PAR_BOUNDS(KEY_ETA_LIMIT, 6, "Eta", "");
    PAR_BOUNDS(KEY_BETA_LIMIT, 8, "Beta", "");

    PAR_INVERSION(KEY_INV_VLOS, 3, "Vlos");
    PAR_INVERSION(KEY_INV_DOPPLER, 4, "Doppler width");
    PAR_INVERSION(KEY_INV_DAMPING, 5, "Damping width");
    PAR_INVERSION(KEY_INV_ETA, 6, "Eta");
    PAR_INVERSION(KEY_INV_CONT, 7, "Continuum");
    PAR_INVERSION(KEY_INV_BETA, 8, "Beta");

    for(int imodel=0; imodel<9; imodel++){
      if(imodel == 7) continue;
      if(!isfinite(params->limits[imodel][0]) 
          || !isfinite(params->limits[imodel][1]) 
          || params->limits[imodel][0] > params->limits[imodel][1]){
        if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
            "parameter bounds must be finite and ordered.\n");
        return -1;
      }
      if(!params->inv[imodel] && !isfinite(params->value_const[imodel])){
        if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
            "fixed parameter values must be finite.\n");
        return -1;
      }
    }
    if((params->magnetic_mode == MAGNETIC_SPHERICAL 
        && params->limits[0][0] < 0.0) 
        || params->limits[4][0] <= 0.0 || params->limits[5][0] < 0.0 
        || params->limits[6][0] < 0.0 || params->limits[8][0] < 0.0 
        || params->limits[8][1] > 1.0){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "physical bounds require Bmod>=0, Dopp>0, Damp>=0, Eta>=0 " 
          "and 0<=Beta<=1.\n");
      return -1;
    }
    if((params->magnetic_mode == MAGNETIC_SPHERICAL && !params->inv[0] 
        && params->value_const[0] < 0.0) 
        || (!params->inv[4] && params->value_const[4] <= 0.0) 
        || (!params->inv[5] && params->value_const[5] < 0.0) 
        || (!params->inv[6] && params->value_const[6] < 0.0) 
        || (!params->inv[8] && (params->value_const[8] < 0.0 
        || params->value_const[8] > 1.0))){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "fixed values require Bmod>=0, Dopp>0, Damp>=0, Eta>=0 " 
          "and 0<=Beta<=1.\n");
      return -1;
    }

    bool has_psrf_parameter = false;
    for(int imodel=0; imodel<4; imodel++){
      if(params->inv[imodel]){
        has_psrf_parameter = true;
        break;
      }
    }
    if(!has_psrf_parameter){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "At least one magnetic parameter or Vlos must be inverted " 
          "for the PSRF convergence test.\n");
      return -1;
    }

    INPUT_FLAG(KEY_HMI_REF, params->HMI_REF, "HMI reference direction");

    return 0;
}

/*--------------------------------------------------------------------------------*/

int RDINPUT(const char filename[], STRUCT_PROFILE_IO *input, STRUCT_PARA *params, 
    STRUCT_DREAM *dream, STRUCT_STK *stokes, STRUCT_MPI *mpi){
  
    /*######################################################################
      Purpose:
        Read and parse the inversion control file.
      Input parameters:
        filename[], the path to the input control file.
        mpi, the MPI configuration.
      Output parameters:
        input, the populated runtime configuration.
     ######################################################################*/
    
    const char *rname = "RDINPUT";
    
    snprintf(input->log_path, sizeof(input->log_path), "./log_%05d", 
        mpi->world_rank);

    if(!FILE_EXIST(filename) && mpi->world_rank==0){
      snprintf(message_buffer, sizeof(message_buffer), 
          "Control file does not exist: %s\n", filename);
      LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);
    }

    FILE *Fa = fopen(filename, "r");
    if(!Fa){
      if(mpi->world_rank==0){
        snprintf(message_buffer, sizeof(message_buffer),       
            "Cannot open control file: %s \n", filename);
        LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);
      }
      return -1;
    }

    STRUCT_KEYS keywords[KEY_TOTAL] ={
      [KEY_LINES] = KEY_REQ("lines"),                                      
      [KEY_REGION_GAP] = KEY_DEF("region_gap","20"),
      [KEY_LINE_WINGS] = KEY_DEF("line_wing_widths","20"),
      [KEY_DREAM_MEMORY] = KEY_DEF("max_dream_memory_gb","8"),
      [KEY_SAMPLE_FILE_LIMIT] = KEY_DEF("max_sample_file_gb","20"),
      [KEY_VERBOSE] = KEY_DEF("verboselv","1"),                            
      [KEY_PROFILE_FORMAT] = KEY_DEF("profile_format","FITS"),
      [KEY_DATA_PATH] = KEY_REQ("data_path"),                              
      [KEY_WAV_PATH] = KEY_DEF("wavelength_path",""),                    
      [KEY_RESULT_PATH] = KEY_DEF("result_path","./inversion_result.fits"),
      [KEY_CACHE_PATH] = KEY_DEF("cache_path","./inversion.cache"),
      [KEY_USE_CACHE] = KEY_DEF("use_cache","NO"),
      [KEY_SAMPLE_PATH] = KEY_DEF("sample_path","./Samples.bin"),
      [KEY_NCHAINS] = KEY_DEF("num_chains","40"),                         
      [KEY_BURNIN_GENERATIONS] = KEY_DEF("burnin_generations","5000"),
      [KEY_SAMPLING_GENERATIONS] = KEY_DEF("sampling_generations","5000"),
      [KEY_MAX_PAIR] = KEY_DEF("max_pair","3"),                            
      [KEY_NCR] = KEY_DEF("ncr","3"),                                      
      [KEY_CR_UPDATE] =  KEY_DEF("cr_update","YES"),                       
      [KEY_PROPOSAL_NOISE_UPDATE] = KEY_DEF("proposal_noise_update","NO"),
      [KEY_SAMPLE_OUTPUT] = KEY_DEF("sample_output","NONE"),
      [KEY_ISLAND_SIZE] = KEY_DEF("island_size","1"),                     
      [KEY_SOL_BOX] = KEY_DEF("sol_box","-1, -1, -1, -1"),                 
      [KEY_VOIGT_PRECISION] = KEY_DEF("voigt_precision", "6"),
      [KEY_NOISE_MODE] = KEY_DEF("noise_mode", "INTENSITY"),
      [KEY_NOISE_LEVEL] = KEY_DEF("noise_level", "1e-3, 1e-3, 1e-3, 1e-3"),
      [KEY_MAGNETIC_MODE] = KEY_DEF("magnetic_mode", "CARTESIAN"),
      [KEY_BZ_LIMIT] = KEY_DEF("Bounds_Bz","-4000, 4000"),
      [KEY_BX_LIMIT] = KEY_DEF("Bounds_Bx","-4000, 4000"),
      [KEY_BY_LIMIT] = KEY_DEF("Bounds_By","0, 4000"),
      [KEY_BMOD_LIMIT] = KEY_DEF("Bounds_Bmod","0, 4500"),
      [KEY_BTHETA_LIMIT] = KEY_DEF("Bounds_ThetaB","0, 1"),
      [KEY_BPHI_LIMIT] = KEY_DEF("Bounds_PhiB","0, 1"),
      [KEY_VLOS_LIMIT] = KEY_DEF("Bounds_Vlos","-20., 20."),                  
      [KEY_DOPPLER_LIMIT] = KEY_DEF("Bounds_Dopp","10, 80"),                  
      [KEY_DAMPING_LIMIT] = KEY_DEF("Bounds_Damp","0.4, 0.8"),                
      [KEY_ETA_LIMIT] = KEY_DEF("Bounds_Eta","2, 90"),                        
      [KEY_BETA_LIMIT] = KEY_DEF("Bounds_Beta","0.1, 0.85"),                  
      [KEY_INV_BZ] = KEY_DEF("Invert_Bz", "YES"),
      [KEY_INV_BX] = KEY_DEF("Invert_Bx", "YES"),
      [KEY_INV_BY] = KEY_DEF("Invert_By", "YES"),
      [KEY_INV_BMOD] = KEY_DEF("Invert_Bmod", "YES"),
      [KEY_INV_BTHETA] = KEY_DEF("Invert_ThetaB", "YES"),
      [KEY_INV_BPHI] = KEY_DEF("Invert_PhiB", "YES"),
      [KEY_INV_VLOS] = KEY_DEF("Invert_Vlos", "YES"),
      [KEY_INV_DOPPLER] = KEY_DEF("Invert_Dopp", "YES"),
      [KEY_INV_DAMPING] = KEY_DEF("Invert_Damp", "YES"),
      [KEY_INV_ETA] = KEY_DEF("Invert_Eta", "YES"),
      [KEY_INV_CONT] = KEY_DEF("Invert_Continuum", "YES"),
      [KEY_INV_BETA] = KEY_DEF("Invert_Beta", "YES"),
      [KEY_HMI_REF] = KEY_DEF("HMI_REF", "NO"),
    };

    char *lines = NULL, key[Key_Length], *value;
    size_t size = 0;
    ptrdiff_t read_status, len_tot, len;
    size_t nspace;
    mpi->verbose_level = 0;
    params->nlines = 0;
    params->nregions = 0;
    params->lines = (STRUCT_MELINE *)calloc(NUM_LINES, 
        sizeof(*params->lines));
    params->regions = (STRUCT_REGION *)calloc(NUM_REGIONS, 
        sizeof(*params->regions));
    if(!params->lines || !params->regions){
      fclose(Fa);
      return -1;
    }
    bool found = false;
    int missing_required = 0;
        
    while((read_status=STR_READ_LINE(&lines, &size, Fa)) > 0){
      len_tot = read_status;
      len = STR_INDEX_CHAR(lines, '=', 1);
      if(len <= 1) continue;

      if(len > Key_Length-1){
        if( mpi->world_rank==0){
          snprintf(message_buffer, sizeof(message_buffer), 
              "\n Keyword is too long.\n%s\nLine skipped.\n", lines);
          LOG_ERROR(ERR_LVL_WARNING, rname, message_buffer);
        }
        continue;
      }

      STR_COPY(key, Key_Length, lines, (size_t)len, true);
      STR_TRIM_RIGHT(key);
        
      value = lines+len+1;
      len_tot -= (len+1);
      found = false;
      nspace = STR_TRIM_LEFT(value);

      if(nspace>0) len_tot -= (ptrdiff_t)nspace;

      if(len_tot>Max_Line_Length-1){
        len_tot = Max_Line_Length-1;

        if(mpi->world_rank==0){
          snprintf(message_buffer, sizeof(message_buffer), 
              "\nValue for keyword %s is too long: %s\n", 
              key, value);
          LOG_ERROR(ERR_LVL_WARNING, rname, message_buffer);
        }
      }

      if(strcmp(key,keywords[KEY_LINES].keyword)==0){
        if(params->nlines >= NUM_LINES){
          if(mpi->world_rank==0){
            snprintf(message_buffer, sizeof(message_buffer), 
                "\nKeyword %s may be specified at most %d times.\n", 
                key, NUM_LINES);
            LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);
          }
        }else{

          STRUCT_MELINE *line = &params->lines[params->nlines];
          if(!PARSE_TWO_DOUBLES(value, &line->wavelength0, 
              &line->lande_factor)){
            if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
                "Invalid spectral-line wavelength or Lande factor.\n");
            free(lines);
            fclose(Fa);
            return -1;
          }
          line->zeeman_shift = (4.6686e-10*line->lande_factor 
              *line->wavelength0*line->wavelength0);
          params->nlines++;
          keywords[KEY_LINES].set = true;
        }
        found = true;
      }else{
        static const char *line_bounds[3] = {
            "Bounds_Dopp_", "Bounds_Damp_", "Bounds_Eta_"};
        static const char *line_invert[3] = {
            "Invert_Dopp_", "Invert_Damp_", "Invert_Eta_"};
        static const char *region_bounds = "Bounds_Beta_";
        static const char *region_invert[2] = {
            "Invert_Continuum_", "Invert_Beta_"};
        int item_id = -1;
        if(INDEXED_KEY(key, "Bounds_Continuum_", &item_id)){
          if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
              "Bounds_Continuum_<region> is not supported because continuum " 
              "bounds are derived from each observed profile.\n");
          free(lines);
          fclose(Fa);
          return -1;
        }
        for(int itype=0; itype<3 && !found; itype++){
          if(INDEXED_KEY(key, line_bounds[itype], &item_id)){
            if(item_id < 0 || item_id >= NUM_LINES 
                || !PARSE_TWO_DOUBLES(value, 
                  &params->lines[item_id].custom_limits[itype][0], 
                  &params->lines[item_id].custom_limits[itype][1])){
              if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
                  "invalid indexed line-parameter bounds.\n");
              free(lines);
              fclose(Fa);
              return -1;
            }
            params->lines[item_id].custom_bounds[itype] = true;
            found = true;
          }else if(INDEXED_KEY(key, line_invert[itype], &item_id)){
            if(item_id < 0 || item_id >= NUM_LINES 
                || INVERSION_VALUE(value, 
                  &params->lines[item_id].custom_inv[itype], 
                  &params->lines[item_id].custom_value[itype]) != 0){
              if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
                  "invalid indexed line inversion control.\n");
              free(lines);
              fclose(Fa);
              return -1;
            }
            params->lines[item_id].custom_inversion[itype] = true;
            found = true;
          }
        }
        if(!found && INDEXED_KEY(key, region_bounds, &item_id)){
            const int itype = 1;
            if(item_id < 0 || item_id >= NUM_REGIONS 
                || !PARSE_TWO_DOUBLES(value, 
                  &params->regions[item_id].custom_limits[itype][0], 
                  &params->regions[item_id].custom_limits[itype][1])){
              if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
                  "invalid indexed region-parameter bounds.\n");
              free(lines);
              fclose(Fa);
              return -1;
            }
            params->regions[item_id].custom_bounds[itype] = true;
            found = true;
        }
        for(int itype=0; itype<2 && !found; itype++){
          if(INDEXED_KEY(key, region_invert[itype], &item_id)){
            if(item_id < 0 || item_id >= NUM_REGIONS 
                || INVERSION_VALUE(value, 
                  &params->regions[item_id].custom_inv[itype], 
                  &params->regions[item_id].custom_value[itype]) != 0){
              if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
                  "invalid indexed region inversion control.\n");
              free(lines);
              fclose(Fa);
              return -1;
            }
            params->regions[item_id].custom_inversion[itype] = true;
            found = true;
          }
        }
        if(found) continue;
        for(int i=1; i<KEY_TOTAL; i++){
          if(strcmp(key,keywords[i].keyword)==0){
            if(keywords[i].set && mpi->world_rank==0){
              snprintf(message_buffer, sizeof(message_buffer), 
                  "\nKeyword %s is redefined; using the last value.\n", key);
              LOG_ERROR(ERR_LVL_WARNING, rname, message_buffer);      
            }
            STR_COPY(keywords[i].line, Max_Line_Length, value, 
                (size_t)len_tot, true);
            keywords[i].set = true;
            found = true;
            break;
          }
        }
      }
      if(!found && mpi->world_rank==0){
        snprintf(message_buffer, sizeof(message_buffer), 
            "\n Unknown keyword: %s\n", key);
        LOG_ERROR(ERR_LVL_WARNING, rname, message_buffer); 
      }
    }

    if(read_status < -1){
      if(mpi->world_rank == 0){
        snprintf(message_buffer, sizeof(message_buffer),
            "Failed to read control file %s (status %td).\n", filename,
            read_status);
        LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);
      }
      free(lines);
      fclose(Fa);
      return -1;
    }

    free(lines);
    fclose(Fa);

    for(int i=0; i<KEY_TOTAL; i++){
      if(keywords[i].required && !keywords[i].set){
        missing_required++;
        if(mpi->world_rank==0){
          snprintf(message_buffer, sizeof(message_buffer), 
              "\nRequired keyword is missing: %.64s\n", 
              keywords[i].keyword);
          LOG_ERROR(ERR_LVL_ERROR, rname, message_buffer);
        }
      }
    }

    if(missing_required > 0) return -1;
    for(int iline=params->nlines; iline<NUM_LINES; iline++){
      for(int itype=0; itype<3; itype++){
        if(params->lines[iline].custom_bounds[itype] 
            || params->lines[iline].custom_inversion[itype]){
          if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
              "indexed line parameter refers to an undefined line.\n");
          return -1;
        }
      }
    }
    if(params->nlines < 1){
      if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
          "At least one spectral line is required.\n");
      return -1;
    }
    for(int iline=0; iline<params->nlines; iline++){
      STRUCT_MELINE *line = &params->lines[iline];
      if(!isfinite(line->wavelength0) || line->wavelength0 <= 0.0 
          || !isfinite(line->lande_factor)){
        if(mpi->world_rank == 0) LOG_ERROR(ERR_LVL_ERROR, rname, 
            "line wavelength must be positive and line values finite.\n");
        return -1;
      }
    }

    return Get_Keys(keywords, input, params, dream, stokes, mpi);
}

/*--------------------------------------------------------------------------------*/
