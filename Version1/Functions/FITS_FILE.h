#ifndef FITS_FILE_h

#define FITS_FILE_h


#include <stdio.h>

#include <math.h>

#include <string.h>

#include "PATH.h"

#include <fitsio.h>


int Fitsread(char *Filename, int Wav_begin, int Wav_end, int Pixel_begin, int Pixel_end, double **Fits_Data);

int Fitswrite(char *Filename, long *naxes, double **image);

int Fitsprint(char *Filename);


#endif /* FITS_FILE_h */
