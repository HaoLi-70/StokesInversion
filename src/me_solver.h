
#pragma once

/*--------------------------------------------------------------------------------*/

#include <stdbool.h>
#include "faddeeva.h"
#include "input_reader.h"

/*--------------------------------------------------------------------------------*/

enum bounds_type {enum_fold, enum_reflect, enum_set};

typedef enum Model_Parameter_Kind{
    PARAM_B0,
    PARAM_B1,
    PARAM_B2,
    PARAM_VLOS,
    PARAM_DOPPLER,
    PARAM_DAMPING,
    PARAM_ETA,
    PARAM_CONTINUUM,
    PARAM_BETA
} MODEL_PARAMETER_KIND;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_SpectralRegion{

    double wavelength_min, wavelength_max;
    int iw_begin, iw_end;
    int continuum_index, beta_index;
    double custom_limits[2][2], custom_value[2];
    bool custom_bounds[2], custom_inversion[2], custom_inv[2];

} STRUCT_REGION;


typedef struct Struct_MEline{

    // center of the wavelength, effective lande factor, a precomputed 
    // coeffcient related to the Zeeman splitting.
    double wavelength0, lande_factor, zeeman_shift;
    // Inclusive wavelength interval where this transition can contribute.
    // A line may contribute to more than one automatically segmented region.
    int iw_begin, iw_end;
    int dopp_index, damp_index, eta_index;
    double custom_limits[3][2], custom_value[3];
    bool custom_bounds[3], custom_inversion[3], custom_inv[3];

}STRUCT_MELINE;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Parameter{

    STRUCT_MELINE *lines;
    STRUCT_REGION *regions;
    
    int nlines, nregions, nmodel;
    MAGNETIC_MODE magnetic_mode;
    // Dynamic model: shared magnetic field and Vlos, per-line Dopp/Damp/Eta,
    // and per-region Continuum/Beta.

    double *step;
    // limits on the parameter
    double limits[MAX_MODEL_PARAMS][2];

    bool inv[MAX_MODEL_PARAMS];
    double value_const[MAX_MODEL_PARAMS];
    MODEL_PARAMETER_KIND kind[MAX_MODEL_PARAMS];
    int *free_index;
    int bx_free_index, by_free_index;
    // number of free parameter.
    int npar;

    double **bounds;
    double *proposal_scale;
    enum bounds_type *type;
    bool *narrower_guess;
    bool HMI_REF;

}STRUCT_PARA;

/*--------------------------------------------------------------------------------*/

typedef struct Struct_Stokes{

    // Total number of wavelength points in all observed spectral regions.
    int nw;
    // the wavelength, input profiles, synthesized profiles, fitting, 
    // best the fit, noise, and the jacobian used for the inversion.
    double *wavelength, *profile, *synthetic, *inv_noise;

    double intensity_threshold;
    // Absolute wavelength gap [Angstrom] used to split spectral regions.
    double region_gap;
    double line_wing_widths;
    NOISE_MODE noise_mode;
    double noise_level[4];
    double *noise;

    STRUCT_FADDEEVA faddeeva;

}STRUCT_STK;

/*--------------------------------------------------------------------------------*/

extern int Milne_Eddington(const double *model_values, STRUCT_STK *stokes,
    const STRUCT_PARA *params);
  
extern int Noise_Init(STRUCT_STK *stokes);

extern int Profile_Prepare(STRUCT_STK *stokes, STRUCT_PARA *params);

extern int Spectral_Lines_Init(STRUCT_STK *stokes, STRUCT_PARA *params);

/*--------------------------------------------------------------------------------*/
