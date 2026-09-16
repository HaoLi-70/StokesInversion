
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdbool.h>
#include <stdint.h>

#include "string_utils.h"
#include "parallel_runtime.h"

/*--------------------------------------------------------------------------------*/

#define NUM_LINES 20
#define NUM_REGIONS 20
#define MAX_MODEL_PARAMS (4+3*NUM_LINES+2*NUM_REGIONS)
#define KEY_REQ(k)  {k,"",false,true}
#define KEY_DEF(k,v){k,v,false,false}

/*--------------------------------------------------------------------------------*/

typedef enum Key_Index{
    KEY_LINES,
    KEY_REGION_GAP,
    KEY_LINE_WINGS,
    KEY_DREAM_MEMORY,
    KEY_SAMPLE_FILE_LIMIT,
    KEY_VERBOSE,
    KEY_PROFILE_FORMAT,
    KEY_DATA_PATH,
    KEY_WAV_PATH,
    KEY_RESULT_PATH,
    KEY_CACHE_PATH,
    KEY_USE_CACHE,
    KEY_SAMPLE_PATH,
    KEY_NCHAINS,
    KEY_BURNIN_GENERATIONS,
    KEY_SAMPLING_GENERATIONS,
    KEY_MAX_PAIR, 
    KEY_NCR,
    KEY_CR_UPDATE,
    KEY_PROPOSAL_NOISE_UPDATE, 
    KEY_SAMPLE_OUTPUT,
    KEY_ISLAND_SIZE,
    KEY_SOL_BOX,
    KEY_VOIGT_PRECISION,
    KEY_NOISE_MODE,
    KEY_NOISE_LEVEL,
    KEY_MAGNETIC_MODE,
    KEY_BZ_LIMIT,
    KEY_BX_LIMIT,
    KEY_BY_LIMIT,
    KEY_BMOD_LIMIT,
    KEY_BTHETA_LIMIT,
    KEY_BPHI_LIMIT,
    KEY_VLOS_LIMIT,
    KEY_DOPPLER_LIMIT,
    KEY_DAMPING_LIMIT,
    KEY_ETA_LIMIT,
    KEY_BETA_LIMIT,
    KEY_INV_BZ,
    KEY_INV_BX,
    KEY_INV_BY,
    KEY_INV_BMOD,
    KEY_INV_BTHETA,
    KEY_INV_BPHI,
    KEY_INV_VLOS,
    KEY_INV_DOPPLER,
    KEY_INV_DAMPING,
    KEY_INV_ETA,
    KEY_INV_CONT,
    KEY_INV_BETA,
    KEY_HMI_REF,
    KEY_TOTAL
}KEY_INDEX;

typedef enum Noise_Mode{
    NOISE_FROM_INTENSITY,
    NOISE_PER_PIXEL,
    NOISE_GLOBAL
}NOISE_MODE;

typedef enum Profile_Format{
    PROFILE_FORMAT_FITS,
    PROFILE_FORMAT_DAT
}PROFILE_FORMAT;

typedef enum Magnetic_Mode{
    MAGNETIC_CARTESIAN,
    MAGNETIC_SPHERICAL
}MAGNETIC_MODE;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Keywords{

    char keyword[Key_Length];
    char line[Max_Line_Length];
    bool set, required;

}STRUCT_KEYS;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Cache{

    // Cache file identity, dimensions, and inversion configuration.
    char magic[4];
    int nx, ny, ncache;
    int x_begin, x_end, y_begin, y_end;
    uint64_t config_hash;

}STRUCT_CACHE;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Subset{

    int coord[2];
    int processed;

}STRUCT_SUBSET;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Profile_IO{

    PROFILE_FORMAT profile_format;
    char data_path[Max_Line_Length];
    char wavelength_path[Max_Line_Length];
    char log_path[Max_Line_Length];
    char cache_path[Max_Line_Length]; 
    char result_path[Max_Line_Length]; 
    char sample_path[Max_Line_Length];

    int sol_box[2][2];

    // of the data file.
    int counts, nx, ny;

    // buffer for profile reading.
    double *result_buffer;
    double *error_buffer;

    MPI_File sample_file;
    bool sample_file_open;
    int sample_max_count;
    int result_nparams;
    long long sample_header_size, sample_record_size;

    // cache file exists or not.
    bool use_cache, cache_reused;
    // cache file header
    STRUCT_CACHE cache_header;
    uint64_t config_hash;
    int *cache;

}STRUCT_PROFILE_IO;

/*--------------------------------------------------------------------------------*/

struct Struct_Parameter;
struct Struct_Dream;
struct Struct_Stokes;

/*--------------------------------------------------------------------------------*/

extern int RDINPUT(const char filename[], STRUCT_PROFILE_IO *input, 
    struct Struct_Parameter *params, struct Struct_Dream *dream, 
    struct Struct_Stokes *stokes, STRUCT_MPI *mpi);

extern int Model_Layout_Init(struct Struct_Parameter *params, 
    struct Struct_Stokes *stokes, STRUCT_MPI *mpi);

/*--------------------------------------------------------------------------------*/
