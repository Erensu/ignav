/*-----------------------------------------------------------------------------
 * pgm.cc
 *
 * Various routines to manipulate PGM files.
 *----------------------------------------------------------------------------*/
#include <carvig.h>

/* global variable -----------------------------------------------------------*/
#define LENGTH 80

/*----------------------------------------------------------------------------*/
static void _getNextString(FILE *fp,char *line)
{
    int i;

    line[0]='\0';
    while (line[0]=='\0') {
        fscanf(fp,"%s",line);
        i=-1;
        do {
            i++;
            if (line[i]=='#') {
                line[i]='\0';
                while (fgetc(fp)!='\n');
            }
        }  while (line[i]!='\0');
    }
}
/*----------------------------------------------------------------------------*/
static void pnmReadHeader(FILE *fp,int *magic,int *ncols, int *nrows,int *maxval)
{
    char line[LENGTH];

    /* Read magic number */
    _getNextString(fp,line);
    if (line[0]!='P')
        trace(2,"Magic number does not begin with 'P', "
                "but with a '%c'",line[0]);
    sscanf(line,"P%d",magic);

    /* Read size, skipping comments */
    _getNextString(fp,line);
    *ncols=atoi(line);
    _getNextString(fp,line);
    *nrows=atoi(line);
    if (*ncols<0||*nrows<0||*ncols>10000||*nrows>10000) {
        trace(2,"The dimensions %d x %d are unacceptable",*ncols,*nrows);
    }
    /* Read maxval, skipping comments */
    _getNextString(fp,line);
    *maxval=atoi(line);
    fread(line,1,1,fp); /* Read newline which follows maxval */

    if (*maxval!=255) {
        trace(2,"Maxval is not 255, but %d",*maxval);
    }
}
/*----------------------------------------------------------------------------*/
static void pgmReadHeader(FILE *fp, int *magic, int *ncols, int *nrows,
                          int *maxval)
{
    pnmReadHeader(fp,magic,ncols,nrows,maxval);
    if (*magic!=5) {
        trace(2,"Magic number is not 'P5', but 'P%d'",*magic);
    }
}
/*----------------------------------------------------------------------------*/
static void ppmReadHeader(FILE *fp,int *magic,int *ncols, int *nrows,int *maxval)
{
    pnmReadHeader(fp,magic,ncols,nrows,maxval);
    if (*magic!=6) {
        trace(2,"Magic number is not 'P6', but 'P%d'",*magic);
    }
}
/*----------------------------------------------------------------------------*/
static void pgmReadHeaderFile(char *fname, int *magic, int *ncols, int *nrows,
                              int *maxval)
{
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"rb"))==NULL) {
        trace(2,"can't open file named '%s' for reading\n",fname);
    }
    /* Read header */
    pgmReadHeader(fp,magic,ncols,nrows,maxval);

    /* Close file */
    fclose(fp);
}
/*----------------------------------------------------------------------------*/
static void ppmReadHeaderFile(char *fname, int *magic, int *ncols, int *nrows,
                              int *maxval)
{
    FILE *fp;

    /* Open file */
    if ((fp = fopen(fname,"rb"))==NULL) {
        trace(2,"Can't open file named '%s' for reading\n",fname);
    }
    /* Read header */
    ppmReadHeader(fp,magic,ncols,nrows,maxval);

    /* Close file */
    fclose(fp);
}
/*----------------------------------------------------------------------------
 * pgmRead
 *
 * NOTE:  If img is NULL, memory is allocated.
 *--------------------------------------------------------------------------*/
extern unsigned char* pgmRead(FILE *fp,unsigned char *img,
                              int *ncols, int *nrows)
{
    unsigned char *ptr;
    int magic, maxval;
    int i;

    /* Read header */
    pgmReadHeader(fp,&magic,ncols,nrows,&maxval);

    /* Allocate memory, if necessary, and set pointer */
    if (img==NULL) {
        ptr=(unsigned char *)malloc(*ncols**nrows*sizeof(char));
        if (ptr==NULL)
            trace(2,"(pgmRead) Memory not allocated");
    }
    else
        ptr=img;

    /* Read binary image data */
    {
        unsigned char *tmpptr=ptr;
        for (i=0;i<*nrows;i++) {
            fread(tmpptr,*ncols,1,fp);
            tmpptr+=*ncols;
        }
    }
    return ptr;
}
/*----------------------------------------------------------------------------
 * pgmReadFile
 *
 * NOTE:  If img is NULL, memory is allocated.
 *----------------------------------------------------------------------------*/
extern unsigned char* pgmReadFile(char *fname,unsigned char *img,
                                  int *ncols, int *nrows)
{
    unsigned char *ptr;
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"rb"))==NULL) {
        trace(2,"Can't open file named '%s' for reading\n",fname);
    }
    /* Read file */
    ptr=pgmRead(fp,img,ncols,nrows);

    /* Close file */
    fclose(fp);
    return ptr;
}
/*----------------------------------------------------------------------------
 * pgmWrite
 *----------------------------------------------------------------------------*/
extern void pgmWrite(FILE *fp,unsigned char *img, int ncols, int nrows)
{
    int i;

    /* Write header */
    fprintf(fp,"P5\n");
    fprintf(fp,"%d %d\n",ncols,nrows);
    fprintf(fp,"255\n");

    /* Write binary data */
    for (i=0;i<nrows;i++) {
        fwrite(img,ncols,1,fp);
        img+=ncols;
    }
}
/*----------------------------------------------------------------------------
 * pgmWriteFile
 *---------------------------------------------------------------------------*/
extern void pgmWriteFile(char *fname, unsigned char *img, int ncols, int nrows)
{
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"wb"))==NULL) {
        trace(2,"(pgmWriteFile) Can't open file named '%s' for writing\n",fname);
    }
    /* Write to file */
    pgmWrite(fp,img,ncols,nrows);

    /* Close file */
    fclose(fp);
}
/*----------------------------------------------------------------------------
 * ppmWrite
 *---------------------------------------------------------------------------*/
extern void ppmWrite(FILE *fp,unsigned char *redimg,unsigned char *greenimg,
                     unsigned char *blueimg,
                     int ncols, int nrows)
{
    int i,j;

    /* Write header */
    fprintf(fp,"P6\n");
    fprintf(fp,"%d %d\n",ncols,nrows);
    fprintf(fp,"255\n");

    /* Write binary data */
    for (j=0;j<nrows;j++) {
        for (i=0;i<ncols;i++) {
            fwrite(redimg,1,1,fp);
            fwrite(greenimg,1,1,fp);
            fwrite(blueimg,1,1,fp);
            redimg++; greenimg++; blueimg++;
        }
    }
}
/*----------------------------------------------------------------------------
 * ppmWriteFileRGB
 *---------------------------------------------------------------------------*/
extern void ppmWriteFileRGB(char *fname, unsigned char *redimg,
                            unsigned char *greenimg,unsigned char *blueimg,
                            int ncols, int nrows)
{
    FILE *fp;

    /* Open file */
    if ((fp=fopen(fname,"wb"))==NULL) {
        trace(2,"(ppmWriteFileRGB) Can't open file named '%s' for writing\n",fname);
    }
    /* Write to file */
    ppmWrite(fp,redimg,greenimg,blueimg,ncols,nrows);

    /* Close file */
    fclose(fp);
}


