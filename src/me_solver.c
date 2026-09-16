
#include "me_solver.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*--------------------------------------------------------------------------------*/

// Light speed
static const double L_C = 299792458.0;

/*--------------------------------------------------------------------------------*/

static void SUM_ADD(double value, double *sum, double *correction){

    double corrected = value-*correction;
    double updated = *sum+corrected;
    *correction = (updated-*sum)-corrected;
    *sum = updated;
}

/*--------------------------------------------------------------------------------*/

int Milne_Eddington(const double *model_values, STRUCT_STK *stokes,
    const STRUCT_PARA *params){
  
    /*######################################################################
      Purpose:
        Calculate the Stokes profiles under M-E atmosphere model (normal 
            Zeeman effect)
      Input parameters:
        model_values, the model parameters.
        stokes, a structure storing the Stokes profiles.
        params, a structure with the model parameters.
      Output parameters:
        stokes, with synthesized Stokes profiles populated.
      Return:
        0 on success.
      Note:
        LL04 5.47, LL04 9.4 9.8, del Toro Iniesta 2003
        Cartesian mode: values[0] Bz; values[1] Bx; values[2] By.
        Spherical mode: values[0] Bmod; values[1] ThetaB; values[2] PhiB.
        values[3] is the shared Vlos. Dopp, Damp, and Eta are indexed per
        spectral line; S0+S1 and S0/(S0+S1) are indexed per spectral region.
        Each line may have a different center and Lande factor.
        At every wavelength, propagation-matrix contributions from overlapping
        lines are accumulated before the Stokes vector is solved.
    ######################################################################*/
  
    static const double L_SqrtPi = 1.77245385090551588191;
    static const double L_Pi = 3.14159265358979323846;

    if(!model_values || !stokes || !params) return -1;

    STRUCT_FADDEEVA *faddeeva = &stokes->faddeeva;

    double bmod, btheta, bphi;
    
    if(params->magnetic_mode == MAGNETIC_CARTESIAN){
      const double bz = model_values[0];
      const double bx = model_values[1];
      const double by = model_values[2];
      bmod = hypot(hypot(bx, by), bz);
      double direction_cosine = bmod > 0.0 ? bz/bmod : 1.0;
      direction_cosine = fmax(-1.0, fmin(1.0, direction_cosine));
      btheta = acos(direction_cosine);
      bphi = bmod > 0.0 ? atan2(by, bx) : 0.0;
    }else{
      bmod = model_values[0];
      btheta = model_values[1];
      bphi = model_values[2];
      if(bmod < 0.0) return -1;
    }

    if(params->HMI_REF) bphi += 0.5*L_Pi; 

    if(!isfinite(bmod) || !isfinite(btheta) || !isfinite(bphi)
        || !isfinite(model_values[3])) return -1;
    const double sin_theta = sin(btheta);
    const double cos_theta = cos(btheta);
    const double sin_theta_sq = sin_theta*sin_theta;
    const double COS_THB_SQ_1 = cos_theta*cos_theta+1;
    const double cos_2phi = cos(2.*bphi);
    const double sin_2phi = sin(2.*bphi);
    const int nw = stokes->nw;
    double *synthetic_i = stokes->synthetic;
    double *synthetic_q = synthetic_i+nw;
    double *synthetic_u = synthetic_q+nw;
    double *synthetic_v = synthetic_u+nw;

    const double *wavelength = stokes->wavelength;

    double phi[3], psi[3];

    for(int iregion=0; iregion<params->nregions; iregion++){
      const STRUCT_REGION *region = &params->regions[iregion];
      const double continuum = model_values[region->continuum_index];
      const double beta = model_values[region->beta_index];
      if(continuum <= 0.0 || beta < 0.0 || beta > 1.0) return -1;
      const double B0 = continuum*beta;
      const double B1 = continuum-B0;

      for(int iw=region->iw_begin; iw<=region->iw_end; iw++){
      double eta[4] = {0.0, 0.0, 0.0, 0.0};
      double rho[4] = {0.0, 0.0, 0.0, 0.0};

      for(int iline=0; iline<params->nlines; iline++){
        const STRUCT_MELINE *line = &params->lines[iline];
        if(iw < line->iw_begin || iw > line->iw_end) continue;

        const double dopp = model_values[line->dopp_index];
        const double damp = model_values[line->damp_index];
        const double opacity = 0.5*model_values[line->eta_index];
        if(dopp <= 0.0 || damp < 0.0 || opacity < 0.0) return -1;
        const double field_shift = bmod*line->zeeman_shift/dopp;
        const double shifted_center = line->wavelength0
            *(1.0+model_values[3]*1e3/L_C);
        const double shift_width = 1e3
            *(wavelength[iw]-shifted_center)/dopp;
        if(!isfinite(field_shift) || !isfinite(shift_width)) return -1;

        for(int q=-1; q<=1; q++){
          double h = 0.0, l = 0.0;
          int iq = q+1;
          if(Faddeeva(shift_width+q*field_shift, damp, &h, &l, faddeeva)
              != 0 || !isfinite(h) || !isfinite(l)) return -1;
          phi[iq] = h/L_SqrtPi;
          psi[iq] = l/L_SqrtPi;
        }

        const double phisum = 0.5*(phi[0]+phi[2]);
        const double phitot = phi[1]-phisum;
        const double psitot = psi[1]-0.5*(psi[0]+psi[2]);
        const double polarized_cos = opacity*sin_theta_sq*cos_2phi;
        const double polarized_sin = opacity*sin_theta_sq*sin_2phi;
        const double longitudinal = opacity*cos_theta;

        eta[0] += opacity*(phi[1]*sin_theta_sq+phisum*COS_THB_SQ_1);
        eta[1] += phitot*polarized_cos;
        eta[2] += phitot*polarized_sin;
        eta[3] += (phi[0]-phi[2])*longitudinal;
        rho[1] += psitot*polarized_cos;
        rho[2] += psitot*polarized_sin;
        rho[3] += (psi[0]-psi[2])*longitudinal;
      }

      double temp = 1+eta[0];
      double tempsq = temp*temp;
      double temp1 = eta[1]*eta[1]+eta[2]*eta[2]+eta[3]*eta[3];
      double temp2 = rho[1]*rho[1]+rho[2]*rho[2]+rho[3]*rho[3];
      double temp3 = eta[1]*rho[1]+eta[2]*rho[2]+eta[3]*rho[3];
      double delta = (tempsq*(tempsq-temp1+temp2)-temp3*temp3);

      if(!isfinite(delta) || delta <= 0.0) return -1;

      double B1D = B1/delta;
      double c1 = (eta[3]*rho[2]-eta[2]*rho[3]);
      double c2 = (eta[1]*rho[3]-eta[3]*rho[1]);
      double c3 = (eta[2]*rho[1]-eta[1]*rho[2]); 
      synthetic_i[iw] = B0+B1D*(temp*(tempsq+temp2));
      synthetic_q[iw] = -B1D*(tempsq*eta[1]+temp*c1+rho[1]*temp3);
      synthetic_u[iw] = -B1D*(tempsq*eta[2]+temp*c2+rho[2]*temp3);
      synthetic_v[iw] = -B1D*(tempsq*eta[3]+temp*c3+rho[3]*temp3);

      }
    }
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

static int LOWER_BOUND(const double *values, int count, double target){

    /*######################################################################
      Purpose:
        Locate the first sorted value that is not smaller than target.
      Input parameters:
        values, sorted array to search.
        count, number of values in the array.
        target, lower-bound search value.
      Return:
        Index in the inclusive range 0 through count.
    ######################################################################*/

    if(!values || count < 0 || !isfinite(target)) return -1;

    int left = 0, right = count;
    while(left < right){
      int middle = left+(right-left)/2;
      if(values[middle] < target) left = middle+1;
      else right = middle;
    }
    return left;
}

/*--------------------------------------------------------------------------------*/

static double MAX_ABS_PARAMETER(const STRUCT_PARA *params, int index){

    /*######################################################################
      Purpose:
        Obtain the largest possible absolute value of one model parameter.
      Input parameters:
        params, inversion controls, fixed values, and parameter bounds.
        index, model parameter index.
      Return:
        Absolute fixed value or maximum absolute bound.
    ######################################################################*/

    if(!params || index < 0 || index >= params->nmodel) return NAN;
    if(!params->inv[index]) return fabs(params->value_const[index]);
    return fmax(fabs(params->limits[index][0]), 
        fabs(params->limits[index][1]));
}

/*--------------------------------------------------------------------------------*/

int Spectral_Lines_Init(STRUCT_STK *stokes, STRUCT_PARA *params){

    /*######################################################################
      Purpose:
        Precompute the wavelength interval potentially affected by each line.
      Input parameters:
        stokes, the sorted wavelength grid.
        params, spectral lines and configured model-parameter bounds.
      Output parameters:
        params, with inclusive wavelength-index bounds for every line.
      Return:
        0 on success; -4 when a configured line does not intersect the
        observed wavelength grid.
      Note:
        The interval includes the maximum velocity and Zeeman shifts plus 
        the configured number of maximum Doppler widths, retaining broad 
        Voigt wings while allowing distant lines to be skipped.
    ######################################################################*/

    if(!stokes || !params || !stokes->wavelength || !params->lines
        || stokes->nw < 1 || params->nmodel < 4
        || params->nmodel > MAX_MODEL_PARAMS || params->nlines < 1
        || params->nlines > NUM_LINES || !isfinite(stokes->line_wing_widths)
        || stokes->line_wing_widths <= 0.0) return -1;
    for(int iw=0; iw<stokes->nw; iw++){
      if(!isfinite(stokes->wavelength[iw])
          || (iw > 0 && stokes->wavelength[iw]
          <= stokes->wavelength[iw-1])) return -1;
    }

    const double wing_widths = stokes->line_wing_widths;
    double max_velocity = MAX_ABS_PARAMETER(params, 3);
    double max_field;
    if(params->magnetic_mode == MAGNETIC_CARTESIAN){
      double component[3];
      for(int i=0; i<3; i++){
        component[i] = MAX_ABS_PARAMETER(params, i);
      }
      max_field = hypot(hypot(component[0], component[1]), component[2]);
    }else{
      max_field = MAX_ABS_PARAMETER(params, 0);
    }

    for(int iline=0; iline<params->nlines; iline++){
      STRUCT_MELINE *line = &params->lines[iline];
      if(!isfinite(line->wavelength0) || line->wavelength0 <= 0.0
          || !isfinite(line->zeeman_shift) || line->dopp_index < 0
          || line->dopp_index >= params->nmodel || line->damp_index < 0
          || line->damp_index >= params->nmodel || line->eta_index < 0
          || line->eta_index >= params->nmodel) return -1;
      double max_doppler = MAX_ABS_PARAMETER(params, line->dopp_index);
      double velocity_extent = line->wavelength0*max_velocity*1e3/L_C;
      double zeeman_extent = fabs(line->zeeman_shift)*max_field*1e-3;
      double wing_extent = wing_widths*max_doppler*1e-3;
      double extent = velocity_extent+zeeman_extent+wing_extent;
      if(!isfinite(max_velocity) || !isfinite(max_field)
          || !isfinite(max_doppler) || max_doppler <= 0.0
          || !isfinite(extent) || extent < 0.0) return -1;
      int begin = LOWER_BOUND(stokes->wavelength, stokes->nw, 
          line->wavelength0-extent);
      int end = LOWER_BOUND(stokes->wavelength, stokes->nw, 
          line->wavelength0+extent);
      if(begin < 0 || end < 0) return -1;
      if(end < stokes->nw && stokes->wavelength[end] <= 
          line->wavelength0+extent) end++;
      int global_end = end-1;
      line->iw_begin = begin;
      line->iw_end = global_end;
      if(line->iw_begin >= stokes->nw || line->iw_end < 0 
          || line->iw_begin > line->iw_end) return -4;
    }

    return 0;
}

/*--------------------------------------------------------------------------------*/

int Noise_Init(STRUCT_STK *stokes){
  
    /*######################################################################
      Purpose:
        Initialize inverse-noise weights when no noise array is provided.
      Input parameters:
        stokes, a structure storing the Stokes profiles.
      Output parameters:
        stokes, a structure storing the Stokes profiles.
      Return:
        0 on success.
      Note:
        GLOBAL and PIXEL modes both use wavelength-dependent noise arrays;
        GLOBAL reuses one array for every observed pixel.
    ######################################################################*/

    if(!stokes || stokes->nw < 1 || !stokes->profile
        || !stokes->inv_noise) return -1;
    const int nw = stokes->nw;
    if((size_t)nw > SIZE_MAX/4U) return -1;
    size_t profile_count = (size_t)nw*4U;
    for(size_t i=0; i<profile_count; i++){
      if(!isfinite(stokes->profile[i])) return -1;
    }
    if(stokes->noise_mode == NOISE_FROM_INTENSITY){
      double iref = 0.0;
      for(int iw=0; iw<nw; iw++){
        if(stokes->profile[iw] > iref) iref = stokes->profile[iw];
      }
      if(iref <= 0.0) return -1;
      for(int istk=0; istk<4; istk++){
        double sigma = stokes->noise_level[istk]*iref;
        if(!isfinite(sigma) || sigma <= 0.0) return -1;
        double inv_variance = 1.0/(sigma*sigma);
        if(!isfinite(inv_variance) || inv_variance <= 0.0) return -1;
        for(int iw=0; iw<nw; iw++){
          stokes->inv_noise[istk*nw+iw] = inv_variance;
        }
      }
    }else if(stokes->noise_mode == NOISE_PER_PIXEL
        || stokes->noise_mode == NOISE_GLOBAL){
      if(!stokes->noise) return -1;
      for(size_t i=0; i<profile_count; i++){
        double sigma = stokes->noise[i];
        if(!isfinite(sigma) || sigma <= 0.0) return -1;
        double inv_variance = 1.0/(sigma*sigma);
        if(!isfinite(inv_variance) || inv_variance <= 0.0) return -1;
        stokes->inv_noise[i] = inv_variance;
      }
    }else return -1;
    
    return 0;
}

/*--------------------------------------------------------------------------------*/

int Profile_Prepare(STRUCT_STK *stokes, STRUCT_PARA *params){
  
    /*######################################################################
      Purpose:
        Validate an observed profile and derive each region's continuum range.
      Input parameters:
        stokes, a structure storing the Stokes profiles.
        params, model parameter limits to update for this profile.
      Output parameters:
        params, with per-region continuum limits derived from Stokes I.
      Return:
        0 on success; 1 when the profile is below the intensity threshold.
    ######################################################################*/

    if(!stokes || !params || stokes->nw < 1 || !stokes->profile
        || !params->regions || !params->step || params->nregions < 1
        || params->nregions > NUM_REGIONS || params->nmodel < 1
        || params->nmodel > MAX_MODEL_PARAMS
        || !isfinite(stokes->intensity_threshold)
        || stokes->intensity_threshold < 0.0) return -1;

    const double *profile_i = stokes->profile;
    if((size_t)stokes->nw > SIZE_MAX/4U) return -1;
    size_t profile_count = (size_t)stokes->nw*4U;
    for(size_t i=0; i<profile_count; i++){
      if(!isfinite(stokes->profile[i])) return -1;
    }
    int next_index = 0;
    for(int iregion=0; iregion<params->nregions; iregion++){
      const STRUCT_REGION *region = &params->regions[iregion];
      if(region->iw_begin != next_index || region->iw_end < region->iw_begin
          || region->iw_end >= stokes->nw || region->continuum_index < 0
          || region->continuum_index >= params->nmodel) return -1;
      next_index = region->iw_end+1;
      double i_sum = 0.0, correction = 0.0, i_max = 0.0;
      int count = region->iw_end-region->iw_begin+1;
      for(int iw=region->iw_begin; iw<=region->iw_end; iw++){
        if(i_max < profile_i[iw]) i_max = profile_i[iw];
        SUM_ADD(profile_i[iw], &i_sum, &correction);
      }
      double i_mean = i_sum/count;
      if(!isfinite(i_mean) || !isfinite(i_max)) return -1;
      if(i_mean <= 0.0 || i_mean < stokes->intensity_threshold) return 1;
      int continuum_index = region->continuum_index;
      double continuum_lower = i_max*0.6;
      double continuum_upper = i_max*1.4;
      if(!isfinite(continuum_lower) || !isfinite(continuum_upper) 
          || continuum_lower <= 0.0 || continuum_lower >= continuum_upper){
        return -1;
      }
      params->limits[continuum_index][0] = continuum_lower;
      params->limits[continuum_index][1] = continuum_upper;
      params->step[continuum_index] = i_max*0.1;
    }
    if(next_index != stokes->nw) return -1;

    return 0;
}

/*--------------------------------------------------------------------------------*/
