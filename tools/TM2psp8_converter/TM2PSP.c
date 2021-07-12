/**
 * @file    TM2PSP.c
 * @brief   This program converts a Improved Troullier-Martins pseudopotential
 *          to a psp format pseudopotential file.
 *
 * @author  Qimen Xu <qimenxu@gatech.edu>
 *          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
 * 
 * Copyright (c) 2018 Material Physics & Mechanics Group at Georgia Tech.
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M_PI 3.14159265358979323846

// TO PRINT IN COLOR
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

typedef struct _PSD_OBJ {
    char   atomtype[8]; // atom type
    double zion; // number of valence electrons
    double *rVloc;     // stores local part of pseudopotential times radius
    double *UdV; 
    double *rhoIsoAtom;       // stores isolated atom electron density
    double *RadialGrid;   // stores the radial grid
    double *SplinerVlocD;  // stores derivative of rVloc from Spline
    double *SplineFitUdV; // derivative of UdV from Spline
    double *SplineFitIsoAtomDen;
    double *rc;     // component pseudopotential cutoff
    double *Gamma;
    double Vloc_0;
    int lmax;       // maximum pseudopotential component
    int lloc;       // local pseudopotential component
    int size;       // size of the arrays storing the pseudopotentials   
    int *ppl;        // number of nonlocal projectors per l
} PSD_OBJ;

// function declarations
void TM2PSP(char *TM_fname, char *PSP_fname, int lloc);
void read_pseudopotential_TM(char *TM_fname, PSD_OBJ *psd, int lloc);
void getYD_gen(double *X, double *Y, double *YD, int len);
void SplineCoeff(double *X,double *Y,double *YD,int len,double *A3,double *A2);
void tridiag_gen(double *A, double *B, double *C, double *D, int len);
void write_psp(char *PSP_fname, PSD_OBJ psd, int lloc);
void print_usage();



int main(int argc, char *argv[]) {
    char TM_filename[128], PSP_filename[128];
    int lloc;
    
    lloc = 0;
    strncpy(TM_filename, "psd_Ge.pot", sizeof(TM_filename));
    strncpy(PSP_filename, "psd_TM_Ge.pot", sizeof(TM_filename));
    
    if (argc == 3) {
        lloc = atoi(argv[1]);
        strncpy(TM_filename, argv[2], sizeof(TM_filename));
        sprintf(PSP_filename, "%s.psp8", TM_filename);
    } else if (argc == 4) {
        lloc = atoi(argv[1]);
        strncpy(TM_filename, argv[2], sizeof(TM_filename));
        strncpy(PSP_filename, argv[3], sizeof(PSP_filename));
    } else if (argc > 4) {
        print_usage();
    }
    
    printf(GRN "============================\n" RESET);
    printf(GRN "lloc = %d\n" RESET, lloc);
    printf(GRN "input filename: %s\n" RESET, TM_filename);
    printf(GRN "output filename: %s\n" RESET, PSP_filename);
    printf(GRN "============================\n" RESET);
    
    TM2PSP(TM_filename,PSP_filename,lloc);

    printf(GRN "======\n" RESET);
    printf(GRN "Done. \n" RESET);
    printf(GRN "======\n" RESET);
    return 0;
}



/**
 * @brief   Read Troulier-Martins potential and convert to psp format.
 */
void TM2PSP(char *TM_fname, char *PSP_fname, int lloc)
{
    PSD_OBJ psd;
    // read Troulier-Martins potential
    read_pseudopotential_TM(TM_fname, &psd, lloc);
    
    printf(YEL "zion = %f \n" RESET, psd.zion);

    // write psp format potential file
    write_psp(PSP_fname, psd, lloc);    
}



/**
 * @brief   Read pseudopotential files (Troulier-Martins format).
 */
void write_psp(char *PSP_fname, PSD_OBJ psd, int lloc)
{
    int i, l;
    FILE *psp_fp = fopen(PSP_fname,"w");
    if (psp_fp == NULL) {
        printf("\nCannot open file \"%s\"\n",PSP_fname);
        exit(EXIT_FAILURE);
    }
    //B     ONCVPSP-3.3.1  r_core=   1.12502   1.12502
    fprintf(psp_fp,"%s     Troulier-Martins  r_core=",psd.atomtype);
    for (i = 0; i <= psd.lmax; i++) {
        fprintf(psp_fp, "   %.5lf", psd.rc[i]);
    }
    fprintf(psp_fp, "\n");
    fprintf(psp_fp, "0.0000      %.4f      000000    zatom,zion,pspd\n",psd.zion);
    fprintf(psp_fp, "     8       2   %d      %d   %d     0    pspcod,pspxc,lmax,lloc,mmax,r2well\n", psd.lmax,lloc,psd.size);
    fprintf(psp_fp, "  0.00000000  0.00000000  0.00000000    rchrg fchrg qchrg\n");
    
    for (i = 0; i <= psd.lmax; i++) {
        fprintf(psp_fp, "%6d", psd.ppl[i]);
    }
    
    for (i = 0; i < 4-psd.lmax; i++) {
        fprintf(psp_fp, "     0");
    }
    fprintf(psp_fp, "     nproj\n");
    fprintf(psp_fp, "     1     1           extension_switch\n");
    // write projectors
    for (l = 0; l <= psd.lmax; l++) {
        if (l == psd.lloc) {
            // write Vloc
            fprintf(psp_fp, "%5d\n",l); // simply setting to lloc is problematic for reading
            fprintf(psp_fp, "     1  0.0000000000000e+00  %.13e\n", psd.rVloc[1]/psd.RadialGrid[1]);
            for (i = 1; i < psd.size; i++) {
                fprintf(psp_fp, "%6d  %.13e  %.13e\n", i+1, psd.RadialGrid[i], psd.rVloc[i]/psd.RadialGrid[i]);
            }
        } else {
            fprintf(psp_fp, "   %-d                         %.13e\n",l,psd.Gamma[l]);
            for (i = 0; i < psd.size; i++) {
                fprintf(psp_fp, "%6d  %.13e  %.13e\n", i+1, psd.RadialGrid[i], psd.UdV[l*psd.size+i]*psd.RadialGrid[i]);
            }
        }
    }

    if (psd.lloc > psd.lmax) {
        // write Vloc
        fprintf(psp_fp, "%5d\n",4); // simply setting to lloc is problematic for reading
        fprintf(psp_fp, "     1  0.0000000000000e+00  %.13e\n", psd.rVloc[1]/psd.RadialGrid[1]);
        for (i = 1; i < psd.size; i++) {
            fprintf(psp_fp, "%6d  %.13e  %.13e\n", i+1, psd.RadialGrid[i], psd.rVloc[i]/psd.RadialGrid[i]);
        }
    }
    

    // write isolated atom electron density (3rd line of )
    for (i = 0; i < psd.size; i++) {
        fprintf(psp_fp, "%6d  %.13e  %.13e  0.0000000000000e+00  0.0000000000000e+00\n", i+1, psd.RadialGrid[i], psd.rhoIsoAtom[i] * 4 * M_PI);
    }
    fprintf(psp_fp, "<INPUT>\n");
    fprintf(psp_fp, "#\n");
    fprintf(psp_fp, "#This is a Troulier-Martins pseudopotential file coverted to psp format\n");
    fprintf(psp_fp, "#   l,   rc,      ep,       ncon, nbas, qcut\n");
    for (l = 0; l <= psd.lmax; l++) {
        fprintf(psp_fp, "%5d   %.5lf  -0.00000   0    0    0.00000\n", l,psd.rc[l]);
    }
    fprintf(psp_fp, "</INPUT>\n");
    
}

/**
 * @brief   Read pseudopotential files (Troulier-Martins format).
 */
void read_pseudopotential_TM(char *TM_fname, PSD_OBJ *psd, int lloc) {
    printf("Reading Troulier-Martins pseudopotential file ...\n");
    int j, k, l, lmax, ppl, count;
    char str[128];
    //char atomtype[8];
    double *U, *rV, *YD, *YD2, *A3, *A2, *B3, *B2, A0, B0, r, dr, vtemp; //, *r2U2dV, *U2dV;
    struct timespec t1, t2;
    FILE *psd_fp;
    
    //PSD_OBJ psd;
    psd_fp = fopen(TM_fname,"r");    
    
    if (psd_fp == NULL) {
        printf("\nCannot open file \"%s\"\n",TM_fname);
        exit(EXIT_FAILURE);
    }
    
    psd->lloc = lloc;

    /* first check the pseudopotential file and find size of pseudopotential data */
    // check atom type
    fscanf(psd_fp, "%s", str);
    stpcpy((*psd).atomtype,str);
    printf("Atom type: %s\n",(*psd).atomtype);
    
    // find size of pseudopotential data
    j = 0; 
    int len_last = 0;
    do {
        len_last = strlen(str);
        fgets(str, 128, psd_fp);
        j++;
    } while (strcmp(str," Radial grid follows\n") != 0 && j < 1e6);
    
    if (strcmp(str," Radial grid follows\n") != 0) {
        printf("\nCannot find radial grid data in pseudopotential file: %s\n", TM_fname);
        exit(EXIT_FAILURE);
    }
    
    // go back one line to read zion
    fseek (psd_fp, -1*(strlen(str)+len_last), SEEK_CUR );
    fscanf(psd_fp, "%lf", &vtemp);
    fscanf(psd_fp, "%lf", &vtemp);
    fscanf(psd_fp, "%lf", &vtemp);
    fscanf(psd_fp, "%lf", &vtemp);
    fscanf(psd_fp, "%lf", &vtemp);
    fscanf(psd_fp, "%lf", &psd->zion);
    printf("zion = %f\n", psd->zion);

    // find size of pseudopotential data
    j = 0;
    do {
        fgets(str, 128, psd_fp);
        j++;
    } while (strcmp(str," Radial grid follows\n") != 0 && j < 1e6);
    

    count = 0;
    do {
        fscanf(psd_fp,"%s",str);
        count++;
    } while (strcmp(str,"Pseudopotential") != 0 && count < 1e6);
    // note: last data read is "Pseudopotential", so count = length + 1;
    
    (*psd).size = count; // one extra value for r = 0
    
    if (j >= 1e6) {
        printf("\nCannot find Pseudopotential data in pseudopotential file: %s\n", TM_fname);
        exit(EXIT_FAILURE);
    }
    
    printf(YEL "length of pseudopotential data for %s: %d\n" RESET, (*psd).atomtype, (*psd).size);
    
    fscanf(psd_fp, "%*[^\n]\n"); 
    fscanf(psd_fp, "%d", &lmax); // read first l value
    while ( fscanf(psd_fp,"%s",str) != EOF && strcmp(str, "Core") != 0) {
        if (strcmp(str,"Pseudopotential") == 0) {
            fscanf(psd_fp, "%*[^\n]\n"); 
            fscanf(psd_fp, "%d", &lmax);
            fscanf(psd_fp, "%*[^\n]\n"); 
        } else {
            // skip current line
            fscanf(psd_fp, "%*[^\n]\n"); 
            continue;
        }
    }
    printf("lmax = %d for %s\n", lmax, (*psd).atomtype);
    (*psd).lmax = lmax;
    
    // reset file pointer to the start of the file
    fseek(psd_fp, 0L, SEEK_SET);  // returns 0 if succeeded, can use to check status
    
    // number of nonlocal projectors per l
    (*psd).ppl = (int *)malloc((lmax + 1) * sizeof(int));
    for (j = 0; j <= lmax; j++) {
        (*psd).ppl[j] = 1; // for Troulier-Martins, there's only 1 projector per l
    }
    ppl = 1;
    
    // allocate memory
    (*psd).RadialGrid = (double *)malloc(count * sizeof(double));
    (*psd).UdV = (double *)malloc((lmax+1) * count * sizeof(double));
    (*psd).rVloc = (double *)malloc(count * sizeof(double));
    (*psd).rhoIsoAtom = (double *)malloc(count * sizeof(double));
    (*psd).rc = (double *)malloc((lmax+1) * sizeof(double));
    (*psd).Gamma = (double *)malloc((lmax+1) * sizeof(double));
    
    U = (double *)malloc((lmax+1) * count * sizeof(double));
    rV = (double *)malloc((lmax+1) * count * sizeof(double));
    
    if ((*psd).RadialGrid == NULL || (*psd).UdV == NULL || 
        (*psd).rVloc == NULL || (*psd).rhoIsoAtom == NULL || 
        (*psd).rc == NULL) 
    {
        printf("\nCannot allocate memory for psd struct members!\n");
        exit(EXIT_FAILURE);
    }
    
    if (U == NULL || rV == NULL) {
        printf("\nCannot allocate memory for U and rV\n");
        exit(EXIT_FAILURE);
    }
    
    // start reading again
    j = 0;
    do {
        fgets(str, 128, psd_fp);
        j++;
    } while (strcmp(str," Radial grid follows\n") != 0 && j < 1e6);
    
    if (j >= 1e6) {
        printf("\nCannot find radial grid data in pseudopotential file: %s\n", TM_fname);
        exit(EXIT_FAILURE);
    }

    // read radial grid
    (*psd).RadialGrid[0]=0.0; 
    for (j = 1; j < count; j++) {
        fscanf(psd_fp, "%lf", &(*psd).RadialGrid[j]);
    }

    j = 0;
    do {
        fgets(str, 128, psd_fp);
        j++;
    } while (strcmp(str," Pseudopotential follows (l on next line)\n") != 0 && j < 1e6);

    if (j >= 1e6) {
        printf("\nCannot find Pseudopotential data in pseudopotential file: %s\n", TM_fname);
        exit(EXIT_FAILURE);
    }

    // read pseudopotential
    while (strcmp(str," Pseudopotential follows (l on next line)\n") == 0) {
        fscanf(psd_fp, "%d", &l);
        if (l == lloc) {
            for (j = 1; j < count; j++) {
                fscanf(psd_fp, "%lf", &vtemp);                 
                //vtemp = vtemp * 0.5 / (*psd).RadialGrid[j];
                vtemp = vtemp * 0.5;
                rV[l*count + j] = vtemp;
                (*psd).rVloc[j] = vtemp;
            }
            rV[l*count] = rV[l*count + 1];
            //(*psd).rVloc[0] = (*psd).rVloc[1]; // TODO: perhaps it's better to set rVloc[0] = 0
            (*psd).rVloc[0] = 0.0;
            (*psd).Vloc_0 = (*psd).rVloc[1] / (*psd).RadialGrid[1]; 
        } else {
            for (j = 1; j < count; j++) {
                fscanf(psd_fp, "%lf", &vtemp);
                //vtemp = vtemp * 0.5 / (*psd).RadialGrid[j];
                vtemp = vtemp * 0.5;
                rV[l*count + j] = vtemp;
            }
            rV[l*count] = rV[l*count + 1]; 
        }
        fgets(str, 128, psd_fp);
        fgets(str, 128, psd_fp);
    }
    
    j = 0;
    while(strcmp(str," Valence charge follows\n") != 0 && j < 1e6) {       
        fgets(str, 128, psd_fp); 
        j++;                 
    }
    
    if (j >= 1e6) {
        printf("\nCannot find Valence charge data in pseudopotential file: %s\n", TM_fname);
        exit(EXIT_FAILURE);
    }

    // read valence charge
    while (strcmp(str," Valence charge follows\n") == 0) {
        for (j = 1; j < count; j++) {
            fscanf(psd_fp, "%lf", &vtemp);
            vtemp = vtemp / (4 * M_PI * (*psd).RadialGrid[j] * (*psd).RadialGrid[j]);
            (*psd).rhoIsoAtom[j] = vtemp;
        }
        (*psd).rhoIsoAtom[0] = (*psd).rhoIsoAtom[1];
        fgets(str, 128, psd_fp);
        fgets(str, 128, psd_fp);
    }

    j = 0;
    while(strcmp(str," Pseudo-wave-function follows (l, zelect, rc)\n") != 0 && j < 1e6) {       
        fgets(str, 128, psd_fp); 
        j++;                 
    }

    if (j >= 1e6) {
        printf("\nCannot find Pseudo-wave-function data in pseudopotential file: %s\n", TM_fname);
        exit(EXIT_FAILURE);
    }

    // read pseudowavefunction 
    while (strcmp(str," Pseudo-wave-function follows (l, zelect, rc)\n") == 0) {
        fscanf(psd_fp, "%d", &l);
        fscanf(psd_fp, "%lf", &vtemp); // skip this value
        fscanf(psd_fp, "%lf", &(*psd).rc[l]);
        for (j = 1; j < count; j++) {
            fscanf(psd_fp, "%lf", &vtemp);
            vtemp = vtemp / (*psd).RadialGrid[j];
            U[l*count + j] = vtemp;
        }
        U[l*count] = U[l*count + 1];
        fgets(str, 128, psd_fp);
        fgets(str, 128, psd_fp);
        if (feof(psd_fp)) break;
    }
    
    // calculate UdV
    for (l = 0; l <= lmax; l++) {
        for (j = 1; j < count; j++) {
            (*psd).UdV[l*count+j] = U[l*count+j] * (rV[l*count+j] - (*psd).rVloc[j]) / (*psd).RadialGrid[j];
        }
        (*psd).UdV[l*count] = (*psd).UdV[l*count+1];
    }
    
    // deallocate memory
    free(rV); rV = NULL;       

    /******************************************************/
    /*           Express U, UdV as a spline               */
    /******************************************************/
    YD = (double *)malloc(count * sizeof(double));
    YD2 = (double *)malloc(count * sizeof(double));
    A3 = (double *)malloc((count-1) * sizeof(double));
    A2 = (double *)malloc((count-1) * sizeof(double));
    B3 = (double *)malloc((count-1) * sizeof(double));
    B2 = (double *)malloc((count-1) * sizeof(double));

    // calculate nonlocal pseudopotential projector coefficients (Gamma)
    for (l = 0; l < lmax+1; l++) {
        // TODO: only 1 projector per l for T-M, add a loop here for ONCV 
        if (l == lloc) {
            (*psd).Gamma[l] = 1.0;
            continue;
        }

        // start timer
        clock_gettime(CLOCK_MONOTONIC, &t1);
        // find polynomail coefficients
        getYD_gen((*psd).RadialGrid, U+l*count, YD, count);
        getYD_gen((*psd).RadialGrid, (*psd).UdV+l*count, YD2, count);
        SplineCoeff((*psd).RadialGrid, U+l*count, YD, count, A3, A2);
        SplineCoeff((*psd).RadialGrid, (*psd).UdV+l*count, YD2, count, B3, B2);
        // integrate r^2 * U * UdV from 0 to Inf (actually to rc+tol is enough)
        vtemp = 0;
        for (k = 0; k < count-1; k++) {
            if ( (*psd).RadialGrid[k] >  (*psd).rc[l]+1.0)
                break;
            r = (*psd).RadialGrid[k];
            dr = (*psd).RadialGrid[k+1] - (*psd).RadialGrid[k];
            A0 = U[l*count+k];
            B0 = (*psd).UdV[l*count+k];
            // use analytical formula to find the integral
            vtemp += (dr * (42.0*A0 * (20.0*B0 * (3.0*r*r + 3.0*r*dr + dr*dr) + dr * (5.0*YD2[k] * (6.0*r*r + 8.0*r*dr + 3*dr*dr) 
                    + dr * (2.0*B2[k] * (10.0*r*r + 15.0*r*dr + 6.0*dr*dr) + B3[k]*dr * (15.0*r*r + 24.0*r*dr + 10.0*dr*dr)))) 
                    + dr * (6.0*YD[k] * (35.0*B0 * (6.0*r*r + 8.0*r*dr + 3.0*dr*dr) + dr * (14.0*YD2[k] * (10.0*r*r + 15.0*r*dr 
                    + 6.0*dr*dr) + dr * (7.0*B2[k] * (15.0*r*r + 24.0*r*dr + 10.0*dr*dr) + 4.0*B3[k]*dr * (21.0*r*r + 35.0*r*dr 
                    + 15.0*dr*dr))))+ dr * (3.0*A2[k] * (28.0*B0 * (10.0*r*r + 15.0*r*dr + 6.0*dr*dr) 
                    + dr * (14.0*YD2[k] * (15.0*r*r + 24.0*r*dr + 10.0*dr*dr) + dr * (8.0*B2[k] * (21.0*r*r + 35.0*r*dr + 15.0*dr*dr) 
                    + 5.0*B3[k]*dr * (28.0*r*r + 48.0*r*dr + 21.0*dr*dr)))) + A3[k]*dr * (42.0*B0 * (15.0*r*r + 24.0*r*dr + 10.0*dr*dr) 
                    + dr * (24.0*YD2[k] * (21.0*r*r + 35.0*r*dr + 15.0*dr*dr) + 5.0*dr * (3.0*B2[k] * (28.0*r*r + 48.0*r*dr + 21.0*dr*dr) 
                    + 2.0*B3[k]*dr * (36.0*r*r + 63.0*r*dr + 28.0*dr*dr))))))))/2520.0;
        }
        (*psd).Gamma[l] = 1.0/vtemp;
        //(*psd).Gamma[l] = 1.0; // TODO: REMOVE AFTER CHECK!
        printf("\nGamma_Jl[%d] = %19.16f\n",l,1.0/vtemp); 
        // stop timer
        clock_gettime(CLOCK_MONOTONIC, &t2);
         printf(BLU "time for calculating Gamma_Jl is " 
                RED "%5.3lf " 
                BLU "ms , l = %d.\n" RESET,
                1e3 *(t2.tv_sec - t1.tv_sec + (double) (t2.tv_nsec - t1.tv_nsec)*1e-9),l);


    }
    
    free(YD); YD = NULL;
    free(YD2); YD2 = NULL;
    free(A3); A3 = NULL;
    free(A2); A2 = NULL;     
    free(B3); B3 = NULL;
    free(B2); B2 = NULL;
    free(U); U = NULL;
    fclose(psd_fp);

}




/**
 * @brief   Calculates derivatives of a tabulated function required for spline interpolation.
 */
void getYD_gen(double *X, double *Y, double *YD, int len) {
    int i;
    double h0,h1,r0,r1,*A,*B,*C;
          
    A = (double *)malloc(sizeof(double)*len);
    B = (double *)malloc(sizeof(double)*len);
    C = (double *)malloc(sizeof(double)*len);
    if (A == NULL || B == NULL || C == NULL) {
        printf("Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    h0 =  X[1]-X[0]; h1 = X[2]-X[1];
    r0 = (Y[1]-Y[0])/h0; r1=(Y[2]-Y[1])/h1;
    B[0] = h1*(h0+h1);
    C[0] = (h0+h1)*(h0+h1);
    YD[0] = r0*(3*h0*h1 + 2*h1*h1) + r1*h0*h0;
               
    for(i=1;i<len-1;i++) {
        h0 = X[i]-X[i-1]; h1=X[i+1]-X[i];
        r0 = (Y[i]-Y[i-1])/h0;  r1=(Y[i+1]-Y[i])/h1;
        A[i] = h1;
        B[i] = 2*(h0+h1);
        C[i] = h0;
        YD[i] = 3*(r0*h1 + r1*h0);
    }
           
    A[i] = (h0+h1)*(h0+h1);
    B[i] = h0*(h0+h1);
    YD[i] = r0*h1*h1 + r1*(3*h0*h1 + 2*h0*h0);
     
    tridiag_gen(A,B,C,YD,len);
    
    free(A); free(B); free(C);                                     
}



/**
 * @brief   Cubic spline coefficients, but without copying Y and YD into A0 and A1, respectively.
 */
void SplineCoeff(double *X,double *Y,double *YD,int len,double *A3,double *A2) {
    int j;
    double dx,dy;
    for (j = 0; j < len-1; j++) {
        dx = 1.0 / (X[j+1] - X[j]);
        dy = (Y[j+1] - Y[j]) * dx;
        // note: there's no need to copy Y or YD into another array
        // A0[j] = Y[j]; 
        // A1[j] = YD[j];
        A2[j] = dx * (3.0 * dy - 2.0 * YD[j] - YD[j+1]);
        A3[j] = dx * dx * (-2.0*dy + YD[j] + YD[j+1]);
    }
}

 
/**
 * @brief   Solves a tridiagonal system using Gauss Elimination.
 */
void tridiag_gen(double *A, double *B, double *C, double *D, int len) {
    int i;
    double b, *F;
    F = (double *)malloc(sizeof(double)*len);
    if (F == NULL) {
        printf("Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }
    
    // Gauss elimination; forward substitution
    b = B[0];
    if (b == 0) {
        printf("Divide by zero in tridiag_gen\n"); 
        exit(EXIT_FAILURE);
    }
    D[0] = D[0]/b;
    for (i = 1; i<len; i++) {
        F[i] = C[i-1] / b;
        b= B[i] - A[i] * F[i];
        if (b == 0) {
            printf("Divide by zero in tridiag_gen\n"); 
            exit(EXIT_FAILURE);
        }
        D[i] = (D[i] - D[i-1] * A[i])/b;
    }
    // backsubstitution 
    for (i = len-2; i >= 0; i--)
        D[i] -= (D[i+1]*F[i+1]);
        
    free(F);
}
 
 
void print_usage() {
    printf("\n");
    printf("Usage: TM2psp8 [lloc] [input fname] [output fname]\n");
    printf("\n");
}



