/*-----------------------------------------------------------------------------
 * jpeg.cc: read/write jpeg file
*
* version : $Revision:$ $Date:$
* history : 2018/11/27  1.0  new
 *----------------------------------------------------------------------------*/
#include <carvig.h>
#include <jpeglib.h>
#include <setjmp.h>

/* save jpeg file-------------------------------------------------------------*/
extern void savejpg(const char* filename, uchar* body, int h, int w, int ch,
                    int quality)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    FILE * outfile;		        /* target file */
    JSAMPROW row_pointer[1];	/* pointer to JSAMPLE row[s] */
    int row_stride;	        	/* physical row width in image buffer */

    /* Step 1: allocate and initialize JPEG compression object */

    cinfo.err=jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    /* Step 2: specify data destination (eg, a file) */
    /* Note: steps 2 and 3 can be done in either order. */

    if ((outfile=fopen(filename,"wb"))==NULL) {
        fprintf(stderr,"can't open %s\n",filename);
        exit(1);
    }
    jpeg_stdio_dest(&cinfo,outfile);

    /* Step 3: set parameters for compression */

    cinfo.image_width =w;
    cinfo.image_height=h;
    cinfo.input_components=ch; /* # of color components per pixel */

    if (ch==3) {
        cinfo.in_color_space=JCS_RGB;  /* colorspace of input image */
    }
    else
        cinfo.in_color_space=JCS_GRAYSCALE;  /* colorspace of input image */

    jpeg_set_defaults(&cinfo);

    /* Now you can set any non-default parameters you wish to.
     * Here we just illustrate the use of quality (quantization table) scaling:
     */
    jpeg_set_quality(&cinfo,quality,TRUE/* limit to baseline-JPEG values */);

    /* Step 4: Start compressor */
    jpeg_start_compress(&cinfo,TRUE);

    /* Step 5: while (scan lines remain to be written) */
    /*           jpeg_write_scanlines(...); */

    /* Here we use the library's state variable cinfo.next_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     * To keep things simple, we pass one scanline per call; you can pass
     * more if you wish, though.
     */
    row_stride=w*ch; /* JSAMPLEs per row in image_buffer */

    while (cinfo.next_scanline<cinfo.image_height) {
        /* jpeg_write_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could pass
         * more than one scanline at a time if that's more convenient.
         */
        row_pointer[0]=&body[cinfo.next_scanline*row_stride];
        (void)jpeg_write_scanlines(&cinfo,row_pointer,1);
    }
    /* Step 6: Finish compression */
    jpeg_finish_compress(&cinfo);
    fclose(outfile);

    /* Step 7: release JPEG compression object */
    jpeg_destroy_compress(&cinfo);

    /* And we're done! */
}

struct my_error_mgr {
    struct jpeg_error_mgr pub;	/* "public" fields */
    jmp_buf setjmp_buffer;	    /* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
/* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
my_error_ptr myerr = (my_error_ptr)cinfo->err;

/* Always display the message. */
/* We could postpone this until after returning, if we chose. */
(*cinfo->err->output_message) (cinfo);

/* Return control to the setjmp point */
longjmp(myerr->setjmp_buffer, 1);
}

/* load jpeg file-------------------------------------------------------------*/
extern int loadjpg(const char* filename, uchar* &body, int &h, int &w, int &ch)
{
    struct jpeg_decompress_struct cinfo;
    struct my_error_mgr jerr;

    FILE * infile;		/* source file */
    JSAMPARRAY buffer;  /* Output row buffer */
    int row_stride;		/* physical row width in output buffer */

    if ((infile=fopen(filename,"rb"))==NULL) {
        fprintf(stderr, "can't open %s\n",filename);
        return 0;
    }
    /* Step 1: allocate and initialize JPEG decompression object */

    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err=jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit=my_error_exit;
    /* Establish the setjmp return context for my_error_exit to use. */
    if (setjmp(jerr.setjmp_buffer)) {
        /* If we get here, the JPEG code has signaled an error. */
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        return 0;
    }
    /* Now we can initialize the JPEG decompression object. */
    jpeg_create_decompress(&cinfo);

    /* Step 2: specify data source (eg, a file) */
    jpeg_stdio_src(&cinfo,infile);

    /* Step 3: read file parameters with jpeg_read_header() */
    (void)jpeg_read_header(&cinfo,TRUE);

    /* Step 4: set parameters for decompression */

    /* we don't need to change any of the defaults set by jpeg_read_header(), so
     * we do nothing here.  */

    /* Step 5: Start decompressor */

    (void)jpeg_start_decompress(&cinfo);

    /* We may need to do some setup of our own at this point before reading
     * the data.  After jpeg_start_decompress() we have the correct scaled
     * output image dimensions available, as well as the output colormap
     * if we asked for color quantization.
     * In this example, we need to make an output work buffer of the right size.
     */
    /* JSAMPLEs per row in output buffer */
    row_stride=cinfo.output_width*cinfo.output_components;

    w =cinfo.output_width;
    h =cinfo.output_height;
    ch=cinfo.output_components;
    body=new uchar[w*h*ch];

    /* Make a one-row-high sample array that will go away when done with image */
    buffer=(*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo,JPOOL_IMAGE,row_stride,1);

    /* Step 6: while (scan lines remain to be read) */
    /*           jpeg_read_scanlines(...); */

    /* Here we use the library's state variable cinfo.output_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     */
    while (cinfo.output_scanline<cinfo.output_height) {
        /* jpeg_read_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could ask for
         * more than one scanline at a time if that's more convenient.
         */
        (void)jpeg_read_scanlines(&cinfo,buffer,1);
        /* Assume put_scanline_someplace wants a pointer and sample count. */

        int row_number=cinfo.output_scanline-1;
        uchar* row=body+row_number*row_stride;
        for (int k=0;k<row_stride;k++) {
            row[k]=buffer[0][k];
        }
    }
    /* Step 7: Finish decompression */
    (void)jpeg_finish_decompress(&cinfo);

    /* Step 8: Release JPEG decompression object */
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);

    /* At this point you may want to check to see whether any corrupt-data
     * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
     */
    /* And we're done! */
    return 1;
}
/* RGB convert to GRAY level-------------------------------------------------*/
static int rgb2gray(unsigned char r,unsigned char g,unsigned char b)
{
    return (int)(r*0.114+g*0.587+b*0.299);
}
/* read jpeg-format image raw data from given image file path----------------
 * args  : char *imgfile  I   image file path
 *         gtime_t time   I   image observation time
 *         img_t *img     IO  image raw data
 *         int flag       I   0: left,1: right
 * return: 1:ok, 0:fail
 * --------------------------------------------------------------------------*/
extern int readjpeg(const char *imgfile,gtime_t time,img_t *img,int flag)
{
    trace(3,"readjpeg: path=%s\n",imgfile);

    img->time=time;

    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *input_file;
    JSAMPARRAY buffer;
    int row_width,i;
    unsigned char *tmp=NULL;

    cinfo.err=jpeg_std_error(&jerr);

    if ((input_file=fopen(imgfile,"rb"))==NULL) {
        fprintf(stderr,"can't open %s\n",imgfile);
        return 0;
    }
    /* initialization of jpeg compression objects */
    jpeg_create_decompress(&cinfo);

    /* specify data source for decompression */
    jpeg_stdio_src(&cinfo,input_file);

    /* read file header, set default decompression parameters */
    jpeg_read_header(&cinfo,TRUE);

    /* start decompressor */
    jpeg_start_decompress(&cinfo);

    /* save image informations */
    img->w=cinfo.output_width;
    img->h=cinfo.output_height;
    img->id++;

    row_width=cinfo.output_width*cinfo.output_components;

    /* image raw data buffer */
    buffer=(*cinfo.mem->alloc_sarray)((j_common_ptr)&cinfo,JPOOL_IMAGE,row_width,1);

    tmp=img->data;

    /* read image raw data */
    while (cinfo.output_scanline<cinfo.output_height) {
        jpeg_read_scanlines(&cinfo,buffer,1);
        if (cinfo.output_components==3) {
            for (i=0;i<cinfo.output_width;i++) {
                tmp[i]=rgb2gray((*buffer)[3*i],(*buffer)[3*i+1],(*buffer)[3*i+2]);
            }
            tmp+=cinfo.output_width;
        }
        else {
            memcpy(tmp,*buffer,row_width);
            tmp+=row_width;
        }
    }
#if 0
    /* show image for debug */
    dipsplyimg(img);
#endif
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(input_file);
    return 1;
}