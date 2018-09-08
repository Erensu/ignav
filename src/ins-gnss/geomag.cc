/*------------------------------------------------------------------------------
* geomag.cc : international geomagnetic reference field functions
*
* reference :
*    [1] https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
*
* version : $Revision: 1.1 $ $Date: 2008/09/05 01:32:44 $
* history : 2018/09/06 1.0 new
*-----------------------------------------------------------------------------*/

/****************************************************************************/
/*                                                                          */
/*     NGDC's Geomagnetic Field Modeling software for the IGRF and WMM      */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Disclaimer: This program has undergone limited testing. It is        */
/*     being distributed unofficially. The National Geophysical Data        */
/*     Center does not guarantee it's correctness.                          */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Version 7.0:                                                         */
/*     - input file format changed to                                       */
/*            -- accept new DGRF2005 coeffs with 0.01 nT precision          */
/*            -- make sure all values are separated by blanks               */
/*            -- swapped n and m: first is degree, second is order          */
/*     - new my_isnan function improves portablility                        */
/*     - corrected feet to km conversion factor                             */
/*     - fixed date conversion errors for yyyy,mm,dd format                 */
/*     - fixed lon/lat conversion errors for deg,min,sec format             */
/*     - simplified leap year identification                                */
/*     - changed comment: units of ddot and idot are arc-min/yr             */
/*     - added note that this program computes the secular variation as     */
/*            the 1-year difference, rather than the instantaneous change,  */
/*            which can be slightly different                               */
/*     - clarified that height is above ellipsoid, not above mean sea level */
/*            although the difference is negligible for magnetics           */
/*     - changed main(argv,argc) to usual definition main(argc,argv)        */
/*     - corrected rounding of angles close to 60 minutes                   */
/*     Thanks to all who provided bug reports and suggested fixes           */
/*                                                                          */
/*                                          Stefan Maus Jan-25-2010         */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Version 6.1:                                                         */
/*     Included option to read coordinates from a file and output the       */
/*     results to a new file, repeating the input and adding columns        */
/*     for the output                                                       */
/*                                          Stefan Maus Jan-31-2008         */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Version 6.0:                                                         */
/*     Bug fixes for the interpolation between models. Also added warnings  */
/*     for declination at low H and corrected behaviour at geogr. poles.    */
/*     Placed print-out commands into separate routines to facilitate       */
/*     fine-tuning of the tables                                            */
/*                                          Stefan Maus Aug-24-2004         */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      This program calculates the geomagnetic field values from           */
/*      a spherical harmonic model.  Inputs required by the user are:       */
/*      a spherical harmonic model data file, coordinate preference,        */
/*      altitude, date/range-step, latitude, and longitude.                 */
/*                                                                          */
/*         Spherical Harmonic                                               */
/*         Model Data File       :  Name of the data file containing the    */
/*                                  spherical harmonic coefficients of      */
/*                                  the chosen model.  The model and path   */
/*                                  must be less than PATH chars.           */
/*                                                                          */
/*         Coordinate Preference :  Geodetic (WGS84 latitude and altitude   */
/*                                  above ellipsoid (WGS84),                */
/*                                  or geocentric (spherical, altitude      */
/*                                  measured from the center of the Earth). */
/*                                                                          */
/*         Altitude              :  Altitude above ellipsoid (WGS84). The   */
/*                                  program asks for altitude above mean    */
/*                                  sea level, because the altitude above   */
/*                                  ellipsoid is not known to most users.   */
/*                                  The resulting error is very small and   */
/*                                  negligible for most practical purposes. */
/*                                  If geocentric coordinate preference is  */
/*                                  used, then the altitude must be in the  */
/*                                  range of 6370.20 km - 6971.20 km as     */
/*                                  measured from the center of the earth.  */
/*                                  Enter altitude in kilometers, meters,   */
/*                                  or feet                                 */
/*                                                                          */
/*         Date                  :  Date, in decimal years, for which to    */
/*                                  calculate the values of the magnetic    */
/*                                  field.  The date must be within the     */
/*                                  limits of the model chosen.             */
/*                                                                          */
/*         Latitude              :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for northern    */
/*                                  hemisphere, negative for the southern   */
/*                                  hemisphere.                             */
/*                                                                          */
/*         Longitude             :  Entered in decimal degrees in the       */
/*                                  form xxx.xxx.  Positive for eastern     */
/*                                  hemisphere, negative for the western    */
/*                                  hemisphere.                             */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "geomag.h"

/* IEEE: only NaN is not equal to itself-------------------------------------*/
#define EARTH_RADIUS 6371.2
#define NaN log(-1.0)
#define FT2KM (1.0/0.0003048)
#define PI 3.141592654
#define RAD2DEG (180.0/PI)

#ifndef SEEK_SET
#define SEEK_SET 0
#endif

#define IEXT 0
#define RECL 81

#define MAXINBUFF RECL+14        /* max size of in buffer */
#define MAXREAD MAXINBUFF-2      /* max to read 2 less than total size (just to be safe) */
#define PATH MAXREAD             /* max path and filename length */

#define EXT_COEFF1 (double)0
#define EXT_COEFF2 (double)0
#define EXT_COEFF3 (double)0

static int interpsh(const double date, double dte1, int nmax1, double dte2, int nmax2,
                    const int gh,
                    const double gh1[MAXCOEFF], const double gh2[MAXCOEFF],
                    double gha[MAXCOEFF], double ghb[MAXCOEFF]);
static int extrapsh(double date, double dte1, int nmax1, int nmax2, int gh,
                    const double gh1[MAXCOEFF], const double gh2[MAXCOEFF],
                    double gha[MAXCOEFF], double ghb[MAXCOEFF]);
static int shval3(const CoordinateSystem coordSys, double flat, double flon,
                  const double elev, const int nmax, const int gh,
                  int iext, double ext1, double ext2, double ext3,
                  const double gha[MAXCOEFF], const double ghb[MAXCOEFF],
                  double *const x, double *const y, double *const z);
static void dihf(const double x, const double y, const double z,
                 double *const d, double *const i, double *const h, double *const f);
#ifndef TARGET_EMBEDDED
static int getshc(FILE *stream, int iflag, long int strec,
                  int nmax_of_gh, const int gh,
                  double gh1[MAXCOEFF], double gh2[MAXCOEFF]);
#endif

/**
 * @param  alt Altitude, in units specified by altUnits.
 * @param  altUnits
 * @param  latitude North latitude, in degrees.
 * @param  longitude East (?) longitude, in degrees.
 * @param  sdate Start date.
 * @param  mdfile Filename of the model file.
 */
extern int get_field_components(BField *const bfield,
                                BFieldModel const*const model,
                                double alt,
                                const Units altUnits,
                                const CoordinateSystem coordSys,
                                const double latitude,
                                const double longitude,
                                const double sdate)
{
    int warn_H,warn_H_strong,warn_P;

    int modelI;       /* Which model (Index) */
    int nmax;

    double warn_H_val,warn_H_strong_val;

    warn_H=0;
    warn_H_val=99999.0;
    warn_H_strong=0;
    warn_H_strong_val=99999.0;
    warn_P=0;

    double minAlt;         /* Minimum height of selected model. */
    double maxAlt;         /* Maximum height of selected model. */

    double dtemp=0.0,ftemp=0.0,htemp=0.0,itemp=0.0;
    double xtemp=0.0,ytemp=0.0,ztemp=0.0;

    double gha[MAXCOEFF]; /* Coefficients of resulting model. */
    double ghb[MAXCOEFF]; /* Coefficients of rate of change model. */

    /* Warn if the date is past end of validity. */
    if ((sdate>model->maxyr)&&(sdate<model->maxyr+1)) {
        fprintf(stderr,"\nWarning: The date %4.2f is out of range,\n"
                       "         but still within one year of model expiration date.\n"
                       "         An updated model file is available before 1.1.%4.0f\n",
                sdate,model->maxyr);
    }
    if (sdate<model->minyr||sdate>model->maxyr+1) {
        return 0;
    }
    /* pick model */
    for (modelI=0;modelI<model->nmodel;modelI++) {
        if (sdate<model->yrmax[modelI]) break;
    }
    /* if beyond end of last model use last model */
    if (modelI == model->nmodel) modelI--;

    /* Get altitude min and max for selected model. */
    minAlt=model->altmin[modelI];
    maxAlt=model->altmax[modelI];

    if (coordSys==kCoordSysGeocentric) {
        /* Add Earth radius to ranges. */
        minAlt+=EARTH_RADIUS;
        maxAlt+=EARTH_RADIUS;
    }
    if (coordSys==kCoordSysGeocentric&&altUnits!=kUnitsKilometers) {
        fprintf(stderr,"get_field_components: altitude must be specified in "
                "kilometers with geocentric coordinate system\n");
        return 0;
    }
    /* Do unit conversions if necessary */
    switch (altUnits) {
        case kUnitsMeters:
            minAlt*=1000.0;
            maxAlt*=1000.0;
            break;

        case kUnitsFeet:
            minAlt*=FT2KM;
            maxAlt*=FT2KM;
            break;

        default:
            break;
    }
    if (alt<minAlt||alt>maxAlt) {
        return 0;
    }
    /* Convert altitude to km */
    switch (altUnits) {
        case kUnitsMeters:
            alt*=0.001;
            break;

        case kUnitsFeet:
            alt/=FT2KM;
            break;

        default:
            break;
    }
    if (latitude<-90||latitude>90||longitude<-180||longitude>180){
        return 0;
    }
    if (model->max2[modelI]==0) {
        nmax=interpsh(sdate,
                      model->yrmin[modelI],model->max1[modelI],
                      model->yrmin[modelI+1],model->max1[modelI+1],3,
                      model->gh1[modelI],model->gh2[modelI],
                      gha,ghb);
        nmax=interpsh(sdate+1,model->yrmin[modelI],model->max1[modelI],
                      model->yrmin[modelI+1],model->max1[modelI+1],4,
                      model->gh1[modelI],model->gh2[modelI],
                      gha,ghb);
    } else {
        nmax=extrapsh(sdate,
                      model->epoch[modelI],
                      model->max1[modelI],model->max2[modelI],3,
                      model->gh1[modelI],model->gh2[modelI],
                      gha,ghb);
        nmax=extrapsh(sdate+1,model->epoch[modelI],model->max1[modelI],
                      model->max2[modelI],4,
                      model->gh1[modelI],model->gh2[modelI],
                      gha,ghb);
    }
    /* Do the first calculations */
    shval3(coordSys,latitude,longitude,alt,nmax,3,
           IEXT,EXT_COEFF1,EXT_COEFF2,EXT_COEFF3,
           gha, ghb,
           &(bfield->x),&(bfield->y),&(bfield->z));
    dihf(bfield->x,
         bfield->y,
         bfield->z,
         &(bfield->d),
         &(bfield->i),
         &(bfield->h),
         &(bfield->f));
    shval3(coordSys,latitude,longitude,alt,nmax,4,
           IEXT,EXT_COEFF1,EXT_COEFF2,EXT_COEFF3,
           gha,ghb,
           &xtemp,&ytemp,&ztemp);
    dihf(xtemp,ytemp,ztemp,&dtemp,&itemp,&htemp,&ftemp);

    bfield->ddot=((dtemp-bfield->d)*RAD2DEG);
    if (bfield->ddot > 180.0) bfield->ddot-=360.0;
    if (bfield->ddot<=-180.0) bfield->ddot+=360.0;
    bfield->ddot*=60.0;

    bfield->idot=((itemp-bfield->i)*RAD2DEG)*60;
    bfield->d*=RAD2DEG;
    bfield->i*=RAD2DEG;
    bfield->hdot=htemp-bfield->h;
    bfield->xdot=xtemp-bfield->x;
    bfield->ydot=ytemp-bfield->y;
    bfield->zdot=ztemp-bfield->z;
    bfield->fdot=ftemp-bfield->f;

    /* Deal with geographic and magnetic poles */
    if (bfield->h<100.0) /* at magnetic poles */
    {
        bfield->d=NaN;
        bfield->ddot=NaN;
        /* while rest is ok */
    }
    if (bfield->h<1000.0) {
        warn_H=0;
        warn_H_strong=1;
        if (bfield->h<warn_H_strong_val) {
            warn_H_strong_val=bfield->h;
        }
    } else if (bfield->h<5000.0&&!warn_H_strong)  {
        warn_H=1;
        if (bfield->h<warn_H_val) {
            warn_H_val=bfield->h;
        }
    }
    /* at geographic poles */
    if (90.0-fabs(latitude)<=0.001) {
        bfield->x=NaN;
        bfield->y=NaN;
        bfield->d=NaN;
        bfield->xdot=NaN;
        bfield->ydot=NaN;
        bfield->ddot=NaN;
        warn_P=1;
        warn_H=0;
        warn_H_strong=0;
        /* while rest is ok */
    }
    return 1;
}
#ifndef TARGET_EMBEDDED
/**
 * Reads model coefficients from a file.
 *
 * @param  model Structure into which model coefficients will be read.
 * @param  mdfile Model file name.
 * @return Zero on failure; non-zero on success.
 */
extern int read_model(BFieldModel *const model, const char mdfile[])
{
    int   lineNum=0;          /* First line will be 1 */
    FILE *stream;
    char  inbuff[MAXINBUFF];
    int   modelI;             /* Index into the current model. */

    inbuff[MAXREAD+1]  ='\0';  /* Just to protect mem. */
    inbuff[MAXINBUFF-1]='\0';  /* Just to protect mem. */

    if (!(stream=fopen(mdfile,"rt"))) {
        fprintf(stderr,"Failed to open \"%s\" for reading.\n",mdfile);
        return 0;
    }
    rewind(stream);

    modelI=-1; /* First model will be 0 */
    while (fgets(inbuff,MAXREAD,stream)) {
        lineNum++;

        /* Ensure record size is correct. */
        if (strlen(inbuff)!=RECL)
        {
            fprintf(stderr,"Corrupt record in file %s on line %d.\n",mdfile,lineNum);
            fclose(stream);
            return 0;
        }
#if 0
        /* old statement Dec 1999 */
		if (!strncmp(inbuff,"    ",4)){ /* If 1st 4 chars are spaces */
#else
        /* New statement Dec 1999 changed by wmd  required by year 2000 models */
        if (!strncmp(inbuff,"   ",3)) /* If 1st 3 chars are spaces */
#endif
        {
            /* New model. */
            modelI++;

            /* If too many headers */
            if (modelI>MAXMOD)
            {
                fprintf(stderr,"Too many models in file %s on line %d.",mdfile,lineNum);
                fclose(stream);
                return 0;
            }
            model->irec_pos[modelI]=ftell(stream);

            /* Get fields from buffer into individual vars.  */
            sscanf(inbuff,"%s%lg%d%d%d%lg%lg%lg%lg",
                   model->name[modelI],
                   &(model->epoch[modelI]),
                   &(model->max1[modelI]),
                   &(model->max2[modelI]),
                   &(model->max3[modelI]),
                   &(model->yrmin[modelI]),
                   &(model->yrmax[modelI]),
                   &(model->altmin[modelI]),
                   &(model->altmax[modelI]));

            /* Compute date range for all models */
            if (modelI==0) { /* If first model */
                model->minyr=model->yrmin[0];
                model->maxyr=model->yrmax[0];
            } else {
                if (model->yrmin[modelI]<model->minyr) {
                    model->minyr=model->yrmin[modelI];
                }
                if (model->yrmax[modelI]>model->maxyr){
                    model->maxyr=model->yrmax[modelI];
                }
            }
        }
    }
    model->nmodel=modelI+1;

    for (modelI=1;modelI<model->nmodel;modelI++) {
        if (model->max2[modelI]==0) {
            getshc(stream,1,model->irec_pos[modelI  ],model->max1[modelI  ],1,
                   model->gh1[modelI],model->gh2[modelI]);
            getshc(stream,1,model->irec_pos[modelI+1],model->max1[modelI+1],2,
                   model->gh1[modelI],model->gh2[modelI]);
        } else {
            getshc(stream,1,model->irec_pos[modelI],model->max1[modelI],1,
                   model->gh1[modelI],model->gh2[modelI]);
            getshc(stream,0,model->irec_pos[modelI],model->max2[modelI],2,
                   model->gh1[modelI],model->gh2[modelI]);
        }
    }
    fclose(stream);
    return 1;
}
#endif

/**
 * Computes the decimal day of year from month, day, year.
 * @author Daniel Bergstrom
 *
 * References:
 *
 * 1. Nachum Dershowitz and Edward M. Reingold, Calendrical Calculations,
 *    Cambridge University Press, 3rd edition, ISBN 978-0-521-88540-9.
 *
 * 2. Claus Tøndering, Frequently Asked Questions about Calendars,
 *    Version 2.9, http://www.tondering.dk/claus/calendar.html
 *
 * @param  month
 * @param  day
 * @param  year
 */
extern double julday(const int month, const int day, const int year)
{
    int days[12]={0,31,59,90,120,151,181,212,243,273,304,334 };
    int leap_year=(((year%4  )==0)&&
                  (((year%100)!=0)||((year%400)==0)));

    double day_in_year=(days[month-1]+day+(month>2?leap_year:0));
    return ((double)year+(day_in_year/(365.0+leap_year)));
}
#ifndef TARGET_EMBEDDED
/**
 * Reads spherical harmonic coefficients from the specified model into an
 * array.
 *
 * @param  stream     Logical unit number
 * @param  iflag      Flag for SV equal to ) or not equal to 0
 *                    for designated read statements
 * @param  strec      Starting record number to read from model
 * @param  nmax_of_gh Maximum degree and order of model
 *
 * @param  gh1 or 2   Schmidt quasi-normal internal spherical
 *                    harmonic coefficients
 * @return 0 on failure, 1 on success.
 *
 *     FORTRAN
 *           Bill Flanagan
 *           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301
 *
 *     C
 *           C. H. Shaffer
 *           Lockheed Missiles and Space Company, Sunnyvale CA
 *           August 15, 1988
 *
 */
static int getshc(FILE *stream, int iflag, long int strec,
                  int nmax_of_gh, const int gh,
                  double gh1[MAXCOEFF],
                  double gh2[MAXCOEFF])
{
    char inbuff[MAXINBUFF];
    char irat[9];
    int ii,m,n,mm,nn;
    int line_num;
    double g,hh;
    double trash;

    if (!(gh==1||gh==2)) {
        fprintf(stderr,"getshc: Fatal: argument gh may only be 1 or 2\n");
        return 0;
    }
    ii=0;
    fseek(stream,strec,SEEK_SET);

    for (nn=1;nn<=nmax_of_gh;nn++) {
        for (mm=0;mm<=nn;mm++) {
            if (iflag==1) {
                fgets(inbuff,MAXREAD,stream);
                sscanf(inbuff,"%d%d%lg%lg%lg%lg%s%d",
                       &n,&m,&g,&hh,&trash,&trash,irat,&line_num);
            }
            else {
                fgets(inbuff,MAXREAD,stream);
                sscanf(inbuff,"%d%d%lg%lg%lg%lg%s%d",
                       &n,&m,&trash,&trash,&g,&hh,irat,&line_num);
            }
            if ((nn!=n)||(mm != m)) {
                fclose(stream);
                return 0;
            }
            ii=ii+1;

            switch(gh) {
                case 1: gh1[ii]=g; break;
                case 2: gh2[ii]=g; break;
            }
            if (m!=0) {
                ii=ii+1;

                switch (gh) {
                    case 1: gh1[ii]=hh; break;
                    case 2: gh2[ii]=hh; break;
                }
            }
        }
    }
    return 1;
}
#endif

/**
 * Extrapolates linearly a spherical harmonic model with a rate-of-change
 * model.
 *
 * @param  date   date of resulting model (in decimal year)
 * @param  dte1   date of base model
 * @param  nmax1  maximum degree and order of base model
 * @param  gh1    Schmidt quasi-normal internal spherical harmonic coefficients
 *                of base model
 * @param  nmax2  maximum degree and order of rate-of-change model
 * @param  gh2    Schmidt quasi-normal internal spherical harmonic coefficients
 *                of rate-of-change model
 *
 * @param  gha    Schmidt quasi-normal internal spherical harmonic coefficients
 * @param  ghb    Schmidt quasi-normal internal spherical harmonic coefficients
 * @param  nmax   maximum degree and order of resulting model
 *
 */
static int extrapsh(double date, double dte1, int nmax1, int nmax2, int gh,
                    const double gh1[MAXCOEFF], const double gh2[MAXCOEFF],
                    double gha[MAXCOEFF], double ghb[MAXCOEFF])
{
    int    nmax;
    int    k,l;
    int    ii;
    double factor;

    if (!(gh==3||gh==4)) {
        fprintf(stderr,"extrapsh: fatal: gh may only be 3 or 4\n");
    }

    factor=date-dte1;
    if (nmax1==nmax2) {
        k=nmax1*(nmax1+2);
        nmax=nmax1;
    }
    else {
        if (nmax1>nmax2) {
            k=nmax2*(nmax2+2);
            l=nmax1*(nmax1+2);
            if (gh==3) {
                for (ii=k+1;ii<=l;ii++) {
                    gha[ii]=gh1[ii];
                }
            } else if (gh==4) {
                for (ii=k+1;ii<=l;ii++) {
                    ghb[ii]=gh1[ii];
                }
            }
            nmax=nmax1;
        }
        else {
            k=nmax1*(nmax1+2);
            l=nmax2*(nmax2+2);
            if (gh==3) {
                for (ii=k+1;ii<=l;ii++) {
                    gha[ii]=factor*gh2[ii];
                }
            } else if (gh==4) {
                for (ii=k+1;ii<=l;ii++) {
                    ghb[ii]=factor*gh2[ii];
                }
            }
            nmax=nmax2;
        }
    }
    if (gh==3) {
        for (ii=1;ii<=k;ii++) {
            gha[ii]=gh1[ii]+factor*gh2[ii];
        }
    } else if (gh==4) {
        for (ii=1;ii<=k;ii++) {
            ghb[ii]=gh1[ii]+factor*gh2[ii];
        }
    }
    return(nmax);
}

/**
 * Interpolates linearly, in time, between two spherical harmonic models.
 *
 * @param date  date of resulting model (in decimal year)
 * @param dte1  date of earlier model
 * @param nmax1 maximum degree and order of earlier model
 * @param gh1   Schmidt quasi-normal internal spherical harmonic coefficients
 *              of earlier model
 * @param dte2  date of later model
 * @param nmax2 maximum degree and order of later model
 * @param gh2   Schmidt quasi-normal internal spherical harmonic coefficients
 *              of internal model
 *
 * @param gha Coefficients of resulting model
 * @param ghb Coefficients of resulting model
 * @param nmax Maximum degree and order of resulting model
 *
 * @author Original Fortran code by
 *         A. Zunde
 *         USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225
 *
 * @author Conversion to C by
 *         C. H. Shaffer
 *         Lockheed Missiles and Space Company, Sunnyvale CA
 *         August 17, 1988
 */
static int interpsh(const double date, double dte1, int nmax1, double dte2,
                    int nmax2, const int gh,
                    const double gh1[MAXCOEFF], const double gh2[MAXCOEFF],
                    double gha[MAXCOEFF], double ghb[MAXCOEFF])
{
    int   nmax;
    int   k, l;
    int   ii;
    double factor;

    factor=(date-dte1)/(dte2-dte1);
    if (nmax1==nmax2) {
        k=nmax1*(nmax1+2);
        nmax=nmax1;
    }
    else {
        if (nmax1>nmax2) {
            k=nmax2*(nmax2+2);
            l=nmax1*(nmax1+2);
            switch(gh) {
                case 3:
                    for (ii=k+1;ii<=l;ii++) {
                        gha[ii]=gh1[ii]+factor*(-gh1[ii]);
                    }
                    break;
                case 4:
                    for (ii=k+1;ii<=l;ii++) {
                        ghb[ii]=gh1[ii]+factor*(-gh1[ii]);
                    }
                    break;
                default:
                    printf("\nError in subroutine extrapsh");
                    break;
            }
            nmax=nmax1;
        }
        else {
            k=nmax1*(nmax1+2);
            l=nmax2*(nmax2+2);
            switch (gh) {
                case 3:
                    for (ii=k+1;ii<=l;ii++) {
                        gha[ii]=factor*gh2[ii];
                    }
                    break;
                case 4:
                    for (ii=k+1;ii<=l;ii++) {
                        ghb[ii]=factor*gh2[ii];
                    }
                    break;
                default:
                    fprintf(stderr,"\nError in subroutine extrapsh");
                    break;
            }
            nmax=nmax2;
        }
    }
    switch(gh) {
        case 3:
            for (ii=1;ii<=k;ii++) {
                gha[ii]=gh1[ii]+factor*(gh2[ii]-gh1[ii]);
            }
            break;
        case 4:
            for (ii=1;ii<=k;ii++) {
                ghb[ii]=gh1[ii]+factor*(gh2[ii]-gh1[ii]);
            }
            break;
        default:
            printf("\nError in subroutine extrapsh");
            break;
    }
    return(nmax);
}

/**
 * Calculates field components from spherical harmonic (sh) models.
 *
 * @param coordSys  indicates coordinate system used
 * @param latitude  north latitude, in degrees
 * @param longitude east longitude, in degrees
 * @param elev      WGS84 altitude above ellipsoid (coordSys is
 *                  kCoordSysGeodetic), or radial distance from earth's center
 *                  (coordSys is kCoordSysGeocentric)
 * @param a2,b2     squares of semi-major and semi-minor axes of the reference
 *                  spheroid used for transforming between geodetic and
 *                  geocentric coordinates or components
 * @param nmax      maximum degree and order of coefficients
 * @param iext      external coefficients flag (=0 if none)
 * @param ext1,2,3  the three 1st-degree external coefficients
 *                  (not used if iext = 0)
 *
 * @param x Output northward component
 * @param y Output eastward component
 * @param z Output vertically-downward component
 */
static int shval3(const CoordinateSystem coordSys, double flat, double flon,
                  const double elev, const int nmax, const int gh,
                  int iext, double ext1, double ext2, double ext3,
                  const double gha[MAXCOEFF], const double ghb[MAXCOEFF],
                  double *const x, double *const y, double *const z)
{
    double earths_radius=6371.2;
    double dtr=0.01745329;
    double slat;
    double clat;
    double ratio;
    double aa,bb,cc,dd;
    double sd;
    double cd;
    double r;
    double a2;
    double b2;
    double rr;
    double fm,fn;
    double sl[14];
    double cl[14];
    double p[119];
    double q[119];
    int ii,j,k,l,m,n;
    int npq;
    int ios;
    double argument;
    double power;

    if (!(gh==3||gh==4)) {
        fprintf(stderr,"shval3: gh may only be 3 or 4\n");
    }
    a2=40680631.59; /* Square of the semi-major axes of the WGS84 reference
	                   sphereoid used for transforming between geodetic and
	                   geocentric coordinates. */
    b2=40408299.98; /* Square of the semi-minor axes of the WGS84 reference
	                   sphereoid used for transforming between geodetic and
	                   geocentric coordinates. */
    ios=0;
    r  =elev;
    argument=flat*dtr;
    slat=sin(argument);
    if ((90.0-flat)<0.001) {
        aa=89.999;            /*  300 ft. from North pole  */
    }
    else {
        if ((90.0+flat)<0.001) {
            aa=-89.999;        /*  300 ft. from South pole  */
        }
        else {
            aa=flat;
        }
    }
    argument=aa*dtr;
    clat    =cos(argument);
    argument=flon*dtr;
    sl[1]=sin(argument);
    cl[1]=cos(argument);

    *x=*y=*z=0;

    sd=0.0;
    cd=1.0;
    l=1;
    n=0;
    m=1;
    npq=(nmax*(nmax+3))/2;
    if (coordSys == kCoordSysGeodetic) {
        aa=a2*clat*clat;
        bb=b2*slat*slat;
        cc=aa+bb;
        argument=cc;
        dd=sqrt(argument);
        argument=elev*(elev+2.0*dd)+(a2*aa+b2*bb)/cc;
        r =sqrt(argument);
        cd=(elev+dd)/r;
        sd=(a2-b2)/dd*slat*clat/r;
        aa=slat;
        slat=slat*cd-clat*sd;
        clat=clat*cd+aa*sd;
    }
    ratio=earths_radius/r;
    argument= 3.0;
    aa=sqrt(argument);
    p[1]=2.0*slat;
    p[2]=2.0*clat;
    p[3]=4.5*slat*slat-1.5;
    p[4]=3.0*aa*clat*slat;
    q[1]=-clat;
    q[2]=slat;
    q[3]=-3.0*clat*slat;
    q[4]=aa*(slat *slat-clat*clat);
    for (k=1;k<=npq;++k) {
        if (n<m) {
            m=0;
            n=n+1;
            argument=ratio;
            power=n+2;
            rr=pow(argument,power);
            fn=n;
        }
        fm=m;
        if (k>=5) {
            if (m==n) {
                argument=(1.0-0.5/fm);
                aa=sqrt(argument);
                j =k-n-1;
                p[k]=(1.0+1.0/fm)*aa*clat*p[j];
                q[k]=aa*(clat*q[j] +slat/fm*p[j]);
                sl[m]=sl[m-1]*cl[1]+cl[m-1]*sl[1];
                cl[m]=cl[m-1]*cl[1]-sl[m-1]*sl[1];
            }
            else {
                argument=fn*fn-fm*fm;
                aa=sqrt(argument);
                argument=((fn-1.0)*(fn-1.0))-(fm*fm);
                bb=sqrt( argument )/aa;
                cc=(2.0 * fn - 1.0)/aa;
                ii=k-n;
                j =k-2*n+1;
                p[k]=(fn+1.0)*(cc*slat/fn*p[ii]-bb/(fn-1.0)*p[j]);
                q[k]=cc*(slat*q[ii]-clat/fn*p[ii])-bb*q[j];
            }
        }
        switch (gh) {
            case 3: aa=rr*gha[l]; break;
            case 4: aa=rr*ghb[l]; break;
        }
        if (m==0) {
            *x+=aa*q[k];
            *z-=aa*p[k];
            l=l+1;
        }
        else {
            switch(gh) {
                case 3: bb=rr*gha[l+1]; break;
                case 4: bb=rr*ghb[l+1]; break;
            }
            cc=aa*cl[m]+bb*sl[m];
            *x=*x+cc*q[k];
            *z=*z-cc*p[k];
            if (clat>0) {
                *y+=(aa*sl[m]-bb*cl[m])*fm*p[k]/((fn+1.0)*clat);
            } else {
                *y+=(aa*sl[m]-bb*cl[m])*q[k]*slat;
            }
            l=l+2;
        }
        m=m+1;
    }
    if (iext!=0) {
        aa=ext2*cl[1]+ext3*sl[1];
        
        *x=*x-ext1*clat +aa  *slat;
        *y=*y+ext2*sl[1]-ext3*cl[1];
        *z=*z+ext1*slat +aa  *clat;
    }
    aa=*x;
    *x=*x*cd+*z*sd;
    *z=*z*cd-aa*sd;
    return ios;
}

/**
 * Computes the geomagnetic declination (d), inclination (i), horizontal field
 * strength (h), and total field strength (f) from x, y, and z.
 *
 * @param x Input northward component.
 * @param y Input eastward component.
 * @param z Input vertically-downward component.
 * @param d Output declination.
 * @param i Output inclination.
 * @param h Output horizontal intensity.
 * @param f Output total intensity.
 *
 */
static void dihf(const double x, const double y, const double z,
                 double *const d, double *const i, double *const h, double *const f)
{
    double sn;
    double hpx;

    sn=0.0001;

    *h=sqrt(x*x+y*y);       /* calculate horizontal intensity */
    *f=sqrt(x*x+y*y+z*z);   /* calculate total intensity */
    if (*f<sn) {
        /* If d and i cannot be determined, set them equal to NaN. */
        *d=NaN;
        *i=NaN;
    }
    else {
        *i=atan2(z,*h);
        if (*h<sn) {
            *d=NaN;
        }
        else {
            hpx=*h+x;
            if (hpx<sn) {
                *d=PI;
            }
            else {
                *d=2.0*atan2(y,hpx);
            }
        }
    }
}