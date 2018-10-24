#include "FITS_FILE.h"


int Fitsread(char *Filename, int Wav_begin, int Wav_end, int Pixel_begin, int Pixel_end, double **Fits_Data){
    
    /***********************************************************************************
     Purpose:
     Read the fits file of Hinode/SP used in the inversion process.
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     17 Oct. 2018, Hao Li
     Input parameters:
     Filename, the name of the file.
     Wav_begin, the first pixel in the wavelength axis used in the inversion.
     Wav_end, the last pixel (not include) in the wavelength axis used in the inversion process.
     Pixel_begin, the first pixel along the slit used in the inversion.
     Pixel_end, the last pixel (not include) along the slit used in the inversion process.
     Output parameters:
     Fits_Data[][], a matrix which saves Stokes profiles of each pixel.
     ***********************************************************************************/
    
    // Fits file pointer, defined in fitsio.h
    fitsfile *fptr;
    
    // CFITSIO status value MUST be initialized to zero!
    int status = 0;
    
    // Total number of HDUs
    int Total_numb;
    
    // The current HDU position
    int hdupos;
    
    // The type of HDU, IMAGE_HDU (0), ASCII_TBL (1), BINARY_TBL (2)
    int hdutype;
    
    // Naxis, the number of axes; bitpix, the bit of each pixel
    int naxis, bitpix;
    
    // The length of each axis, and the max dimension of the array is 10.
    long naxes[3];
    
    long fpixel[3] = {1,1,1};
    
    fpixel[0]=1;
    
    fpixel[1]=Pixel_begin+1;
    
    int Num_wav=(Wav_end-Wav_begin);
    
    int Num_Pixel=(Pixel_end-Pixel_begin);

    int Stokes=0;

    int i,j;

    if (!fits_open_file(&fptr, Filename, READONLY, &status)){
        
        // Get the total number of HDUs
        fits_get_num_hdus(fptr, &Total_numb, &status);
        
        // Get the current HDU position
        fits_get_hdu_num(fptr, &hdupos);
       
        // Get the HDU type
        fits_get_hdu_type(fptr, &hdutype, &status);
        
        // Get the parameters of the HDU. Bitpix, the bit of pixel; Naxis, the number of axes; naxes[], the number of each axis
        fits_get_img_param(fptr, 10, &bitpix, &naxis, naxes, &status);

        double *pixels;
        
        pixels = (double *) malloc(naxes[0]* Num_Pixel * sizeof(double));
        
        for (Stokes=0; Stokes<4; Stokes++) {
            
            fpixel[2]=Stokes+1;
            
            fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0] * Num_Pixel, NULL, pixels, NULL, &status);
            
            for (i=0; i<Num_Pixel; i++)
                
                for (j=0; j<Num_wav; j++)

                    Fits_Data[i][j+Stokes*Num_wav]=pixels[j+i*naxes[0]+Wav_begin];
            
        }
        
        free(pixels);
        
    }
    
    return status;
    
}


int Fitswrite(char *Filename, long *naxes, double **image){
    
    /***********************************************************************************
     Purpose:
     Create a FITS, containing an image.
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     17 Oct. 2018, Hao Li
     Input parameters:
     Filename, the name of the file.
     naxes[2], the wide and rows of the image.
     map[][], the inversion result.
     ***********************************************************************************/
    
    // Fits file pointer, defined in fitsio.h
    fitsfile *fptr;
    
    int status=0;
    
    long fpixel = 1, nelements;
    
    int naxis = 2;
    
    //create a new file
    fits_create_file(&fptr, Filename, &status);
    
    // Create the primary array image
    fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    
    // Write a keyword; must pass the ADDRESS of the value
    char *str1="B";
    
    char *str2="theta";
    
    char *str3="phi";
    
    char *str4="V_los";
    
    char *str5="Lambda_D";
    
    char *str6="S";
    
    char *str7="eta";
    
    char *str8="beta";
    
    char *str9="a";
    
    char *str10="B_z";
    
    char *str11="B_x";
    
    char *str12="B_y";

    fits_update_key(fptr, TSTRING, "NAXIS1[0]", str1,"(G)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[2]", str2,"(rad)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[4]", str3,"(rad)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[6]", str4,"(m/s)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[8]", str5,"(m)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[10]", str6,"()", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[12]", str7,"()", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[14]", str8,"()", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[16]", str9,"(Lambda_D)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[18]", str10,"(G)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[20]", str11,"(G)", &status);
    
    fits_update_key(fptr, TSTRING, "NAXIS1[22]", str12,"(G)", &status);

    // number of pixels to write
    nelements = naxes[0] * naxes[1];
    
    // Write the data to the image
    fits_write_img(fptr, TDOUBLE, fpixel, nelements, image[0], &status);
    
    // close the file
    fits_close_file(fptr, &status);
    
    // print out any error messages
    fits_report_error(stderr, status);
    
    return status;
    
}


int Fitsprint(char *Filename){
    
    /***********************************************************************************
     Purpose:
     Read the fits file, print some data in the file.
     Licensing:
     This code is distributed under the GNU LGPL license.
     Modified:
     17 Oct. 2018, Hao Li
     Input parameters:
     Filename, the name of the file.
     ***********************************************************************************/

    // Fits file pointer, defined in fitsio.h
    fitsfile *fptr;
    
    // CFITSIO status value MUST be initialized to zero!
    int status = 0;
    
    // Naxis, the number of axes; bitpix, the bit of each pixel
    int naxis, bitpix;
    
    // The length of each axis, and the max dimension of the array is 10.
    long naxes[10];
    
    // The position of the first pixel to print
    long fpixel[10] = {1,1,1,1,1,1,1,1,1,1};
    
    // The data of array returned
    double *pixels;
    
    // The number of keywords
    int nkeys;
    
    // Fits header card, FLEN_CARD (81) is the max length of the card
    char card[FLEN_CARD];
    
    // The format to print
    char format[20], hdformat[20];
    
    // The current HDU position
    int hdupos;
    
    // The type of HDU, IMAGE_HDU (0), ASCII_TBL (1), BINARY_TBL (2)
    int hdutype;
    
    char keyname[FLEN_KEYWORD], colname[FLEN_VALUE], coltype[FLEN_VALUE];
    
    int ncols;
    
    long nrows;
    
    char LPRT[7];
    
    char yes[3]="y";
    
    char oneline[7]="1line";
    
    int ii, jj, Total_numb;
    
    if (!fits_open_file(&fptr, Filename, READONLY, &status)){
        
        // Get the total number of HDUs
        fits_get_num_hdus(fptr, &Total_numb, &status);
        
        printf("\n Total number of HDUs is %d  ", Total_numb);
        
        // Get the current HDU position
        fits_get_hdu_num(fptr, &hdupos);
        
        // Main loop for each HDU
        for (; !status; hdupos++){
            
            // Get the HDU type
            fits_get_hdu_type(fptr, &hdutype, &status);
            
            printf("\n Position of the currently opened HDU #%d  ", hdupos);
            
            //image HDU
            if (hdutype == IMAGE_HDU){
                
                // Get the parameters of the HDU. Bitpix, the bit of pixel; Naxis, the number of axes; naxes[], the number of each axis
                fits_get_img_param(fptr, 10, &bitpix, &naxis, naxes, &status);
                
                // Set the default output format string
                if (bitpix > 0) {
                    
                    strcpy(hdformat, " %7d");
                    
                    strcpy(format,   " %7.0f");
                    
                } else {
                    
                    strcpy(hdformat, " %15d");
                    
                    strcpy(format,   " %15.5f");
                    
                }
                
                printf("Array:  NAXIS = %d,  BITPIX = %d\n", naxis, bitpix);
                
                for (ii = 0; ii < naxis; ii++)
                    
                    printf("   NAXIS%d = %ld\n",ii+1, naxes[ii]);
                
                // Get the number of keywords
                fits_get_hdrspace(fptr, &nkeys, NULL, &status);
                
                printf("Header listing for HDU #%d:\n", hdupos);
                
                // Read and print each keywords
                for (ii = 1; ii <= nkeys; ii++) {
                    
                    if (fits_read_record(fptr, ii, card, &status))
                        
                        break;
                    
                    printf("%s\n", card);
                    
                }
                
                // Terminate listing with END
                printf("END\n\n");
                
                /////////////////////////////////// Print the array ////////////////////////////////////
                switch (naxis) {
                        
                    case 0:
                        
                        break;
                        
                    case 1:
                        
                        // 1 Dimensional image
                        printf("1D image\n");
                        
                        // Allocate the memory
                        pixels = (double *) malloc(naxes[0] * sizeof(double));
                        
                        if (pixels == NULL) {
                            
                            printf("Memory allocation error\n");
                            
                            return(1);
                            
                        }
                        
                        // Read the pixels, and jump out of loop on error
                        if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status) )
                            
                            break;
                        
                        // Print the values
                        for (ii = 0; ii < naxes[0]; ii++)
                            
                            printf(format, pixels[ii]);
                        
                        printf("\n");
                        
                        // Free the memory
                        free(pixels);
                        
                        break;
                        
                    case 2:
                        
                        // 2 Dimensional image
                        printf("2D image\n");
                        
                        // Allocate the memory
                        pixels = (double *) malloc(naxes[0] * sizeof(double));
                        
                        if (pixels == NULL) {
                        
                            printf("Memory allocation error\n");
                            
                            return(1);
                        
                        }
                        
                        for (fpixel[1] = 1; fpixel[1] <= naxes[1] ; fpixel[1]++){
                            
                            // Read the pixels, and jump out of loop on error
                            if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status) )
                            
                                break;
                            
                            // Print the row number and the values of each row
                            printf("row number: %4ld \n",fpixel[1]);
                            
                            for (ii = 0; ii < naxes[0]; ii++)
                                
                                printf(format, pixels[ii]);

                            printf("\n");
                        
                        }
                        
                        free(pixels);
                        
                        break;
                        
                    default:
                        
                        // Multi-demensional image
                        printf("2D image of a multi-demensional-array\n");
                        
                        if (naxes[0] * naxes[1]<=5000) {
                        
                            // Allocate the memory
                            pixels = (double *) malloc(naxes[0] * naxes[1] * sizeof(double));
                            
                            if (pixels == NULL) {
                            
                                printf("Memory allocation error\n");
                                
                                return(1);
                            
                            }
                            
                            printf("%ld %ld \n",fpixel[0],fpixel[1]);
                            
                            // Read the pixels, and jump out of loop on error
                            if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], NULL, pixels, NULL, &status) )
                            
                                break;
                            
                            printf("status = %d \n",status);
                            
                            jj=0;//naxes[1]-1;
                            
                            for (jj=0; jj<naxes[1]; jj++) {
                                
                                for (ii = 0; ii < naxes[0]; ii++)
                                
                                    printf(format, pixels[ii+jj*naxes[0]]);
                                
                                printf("\n");
                            
                            }
                            
                            // terminate line
                            
                            free(pixels);
                        }else {
                         
                            printf("The image is really large, (%ld * %ld )\n",naxes[0], naxes[1]);
                            
                            printf("Do you want to print ? (y/n/1line)\n");
                            
                            scanf("%s",LPRT);
                            
                            printf("%s \n",LPRT);
                            
                            if ( !strcmp (LPRT,yes) ) {
                                
                                // Allocate the memory
                                pixels = (double *) malloc(naxes[0] * naxes[1] * sizeof(double));
                                
                                if (pixels == NULL) {
                                
                                    printf("Memory allocation error\n");
                                    
                                    return(1);
                                
                                }
                                
                                printf("%ld %ld \n",fpixel[0],fpixel[1]);
                                
                                // Read the pixels, and jump out of loop on error
                                if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], NULL, pixels, NULL, &status) )
                                    
                                    break;
                                
                                printf("status = %d \n",status);
                                
                                jj=0;//naxes[1]-1;
                                
                                for (jj=0; jj<naxes[1]; jj++) {
                                
                                    for (ii = 0; ii < naxes[0]; ii++)
                                    
                                        printf(format, pixels[ii+jj*naxes[0]]);
                                    
                                    printf("\n");
                                
                                }
                                
                                // terminate line
                                free(pixels);
                            
                            }else if ( !strcmp (LPRT,oneline) ) {
                                
                                // Allocate the memory
                                pixels = (double *) malloc(naxes[0]  * sizeof(double));
                                
                                if (pixels == NULL) {
                                
                                    printf("Memory allocation error\n");
                                    
                                    return(1);
                                
                                }
                                
                                // Read the pixels, and jump out of loop on error
                                if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status) )
                                    
                                    break;
                                
                                jj=0;//naxes[1]-1;
                                
                                for (ii = 0; ii < naxes[0]; ii++)
                                
                                    printf(format, pixels[ii+jj*naxes[0]]);
                                
                                printf("\n");
                                
                                // terminate line
                                free(pixels);
                            
                            }
                        
                        }
                        
                        break;
                
                }
                
            }else{  // a table HDU

                fits_get_num_rows(fptr, &nrows, &status);
                
                fits_get_num_cols(fptr, &ncols, &status);
                
                if (hdutype == ASCII_TBL)
                    
                    printf("ASCII Table:  ");
                
                else
                
                    printf("Binary Table:  ");
                
                printf("%d columns x %ld rows\n", ncols, nrows);
                
                printf(" COL NAME             FORMAT\n");
                
                for (ii = 1; ii <= ncols; ii++){
                    
                    fits_make_keyn("TTYPE", ii, keyname, &status); // make keyword
                    
                    fits_read_key(fptr, TSTRING, keyname, colname, NULL, &status);
                    
                    fits_make_keyn("TFORM", ii, keyname, &status); // make keyword
                    
                    fits_read_key(fptr, TSTRING, keyname, coltype, NULL, &status);
                    
                    printf(" %3d %-16s %-16s\n", ii, colname, coltype);
                
                }
            
            }
            
            //      if (single)break;  // quit if only listing a single HDU
            fits_movrel_hdu(fptr, 1, NULL, &status);  // try move to next ext
        
        }
        
        if (status == END_OF_FILE) status = 0; // Reset normal error
        
        fits_close_file(fptr, &status);
        
    }
    
    if (status) fits_report_error(stderr, status); // print any error message
    
    return status;

}


