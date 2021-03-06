/*
 * template_draw.c - fast line drawing for the template display
 *
 * Andrew Whitwham, January 2010
 * Wellcome Trust Sanger Institute
 *
 */

#include "template_draw.h"

/*
 * On Windows we need ckfree/ckalloc for allocating the XImage as the
 * X emulation layer inside Tk/win internally uses ckfree/ckalloc (we're
 * not meant to be calling these functions directly ourselves) and mixing
 * ckfree with malloc causes crashes.
 *
 * On unix though the XImage calls aren't part of Tk, so we get the same
 * error if we switch to using ckfree/ckalloc instead. Hence we redefine
 * them here so we can choose which functions to call.
 */
#ifndef _WIN32
#    ifdef ckfree
#        undef ckfree
#    endif
#    define ckfree(x) free((x))
#    ifdef ckalloc
#        undef ckalloc
#    endif
#    define ckalloc(x) malloc((x))
#endif

/*
    Byte order is significant for the XImage structure.
*/
static int get_byte_order(void) {
    union {
    	char c[sizeof(short)];
    	short s;
    } order;

    order.s = 1;
    
    if ((1 == order.c[0])) {
    	return LSBFirst;
    } else {
    	return MSBFirst;
    }
}

/*
    Encode colour from the RGB values.
*/
static void code_colour(image_t *image, unsigned int *r, unsigned int *g, unsigned int *b) {
    Visual *vis;
    double r_ratio;
    double g_ratio;
    double b_ratio;

    vis = DefaultVisual(image->dis, image->screen);

#ifdef _WIN32
    // The Tk XImage support uses a BGR layout, but allocating colours
    // here always yields red and blue swaps as if internally it is RGB.
    // I know not why this is the case, only that for now this hack
    // resolves it.
    //
    // http://sourceforge.net/p/tktoolkit/bugs/1049/ touches on this
    // subject, but it was discounted as being a Tk bug.

    unsigned int t = *r;
    *r = *b;
    *b = t;
#endif

    r_ratio = vis->red_mask   / 255.0;
    g_ratio = vis->green_mask / 255.0;
    b_ratio = vis->blue_mask  / 255.0;

    *r *= r_ratio;
    *g *= g_ratio;
    *b *= b_ratio;

    *r &= vis->red_mask;
    *g &= vis->green_mask;
    *b &= vis->blue_mask;
}


static uint32_t get_colour32(image_t *image, unsigned int red, unsigned int green, unsigned int blue) {
    uint32_t colour;

    code_colour(image, &red, &green, &blue);
    colour = red | green | blue;

    return colour;
}


static uint16_t get_colour16(image_t *image, unsigned int red, unsigned int green, unsigned int blue) {
    uint16_t colour;

    code_colour(image, &red, &green, &blue);
    colour = red | green | blue;

    return colour;
}


/*
    Return a colour from the colour array.
*/
static void *return_colour(image_t *image, int col_no) {

    if (image->depth >= 24) {
    	uint32_t *col = (uint32_t *)(image->map_col);
    	return (void *)(col + col_no);
    } else if (image->depth >= 15) {
    	uint16_t *col = (uint16_t *)(image->map_col);
    	return (void *)(col + col_no);
    } else {
    	return NULL; // nothing lower than 16 bit colour
    }
}
		

/*
    Add a colour to the colour array.
    Returns the colour number on success or -1 on failure
*/
int add_colour(image_t *image, unsigned int red, unsigned int green, unsigned int blue) {

    if (image->depth >= 24) {
    	uint32_t *map;

    	if (image->colour_count == image->colour_size) {
    	    image->colour_size *= 2;
    	    image->map_col = realloc(image->map_col, image->colour_size * sizeof(uint32_t));
    	}

    	map = (uint32_t *)(image->map_col);
    	map[image->colour_count] = get_colour32(image, red, green, blue);
    } else if (image->depth >= 15) {
    	uint16_t *map;

    	if (image->colour_count == image->colour_size) {
    	    image->colour_size *= 2;
    	    image->map_col = realloc(image->map_col, image->colour_size * sizeof(uint16_t));
    	}

    	map = (uint16_t *)(image->map_col);
    	map[image->colour_count] = get_colour16(image, red, green, blue);
    } else {
    	return -1;
    }

    image->colour_count++;

    return (image->colour_count - 1);
}


/*
   Initialises an image struct (but not the image itself).
   To free memory, use image_destroy (see below).
*/
image_t *initialise_image(Display *dis) {
    image_t *image = NULL;
    
    if (NULL == (image = malloc(sizeof(image_t)))) return NULL;

    image->dis    = dis;
    image->screen = DefaultScreen(dis);
    image->depth  = DefaultDepth(image->dis, image->screen);
    image->colour_size = 256;
    image->colour_count = 0;
    image->buf = NULL;
	
    if (image->depth >= 24) {
    	uint32_t *map_col = malloc(256 * sizeof(uint32_t));
	image->map_col = map_col;
    } else if (image->depth >= 15) {
    	uint16_t *map_col = malloc(256 * sizeof(uint16_t));
	image->map_col = map_col;
    } else {
    	fprintf(stderr, "Min 16 bit colour needed\n");
	free(image);
	return NULL;
    }
	
    image->img = 0;
    
    return image;
}


/*
   Allocate the memory that drawing is done on to, initialising it with
   a background colour.
   Free using image_remove (see below).
*/	
int create_image_buffer(image_t *image, int width, int height, int bg_colour) {
    size_t i;
    size_t size;

    size = width * height;
    
    image->width  = width;
    image->height = height;
	
    if (image->depth >= 24) {
    	uint32_t *buf = (uint32_t *)ckalloc(size * sizeof(uint32_t));
	uint32_t *col = return_colour(image, bg_colour);

    	if (NULL == buf) return 0;

    	for (i = 0; i < size; i++) {
	    buf[i] = *(uint32_t *)col;
    	}

	image->buf = buf;
    } else if (image->depth >= 15) {
    	uint16_t *buf = (uint16_t *)ckalloc(size * sizeof(uint16_t));
	uint16_t *col = return_colour(image, bg_colour);

    	if (NULL == buf) return 0;

    	for (i = 0; i < size; i++) {
	    buf[i] = *(uint16_t *)col;
    	}

	image->buf = buf;
    } else { /* less than 16-bit color */
	return 0;
    }

    return 1;
}


int draw_line(image_t *image, int x1, int x2, int y, int colour) {
    int i;
    int length;

    // check boundary conditions	
    if (y >= image->height || y < 0) return 0;
    if ((x1 < 0 && x2 < 0) || (x1 >= image->width && x2 >= image->width)) return 0;

    if (x1 > x2) {
    	int tmp = x1;
	x1 = x2;
	x2 = tmp;
    }
    
    if (x1 < 0) x1 = 0;
    if (x2 >= image->width) x2 = image->width - 1;

    length = x2 - x1;

    y *= image->width; // y offset
    x1 += y;

    // draw our line
    i = x1;
	
    if (image->depth >= 24) {
	uint32_t *buf = (uint32_t *)image->buf;
    	uint32_t *col = ((uint32_t *)(image->map_col) + colour);
    
    	do {
    	    buf[i] = *col;
    	} while (i++ < (x1 + length));
		
    } else if (image->depth >= 15) {
	uint16_t *buf = (uint16_t *)image->buf;
	uint16_t *col = ((uint16_t *)(image->map_col) + colour);
    
    	do {
    	    buf[i] = *col;
    	} while (i++ < (x1 + length));
    } else { /* less than 16-bit color */ 
	return 0;
    }

    return 1;
}


void create_image_from_buffer(image_t *image) {

    if (image->depth >= 24) {
    	image->img = XCreateImage(image->dis, CopyFromParent, image->depth, ZPixmap, 0, 
		    (char *)image->buf, image->width, image->height, 32, 0);
    } else if (image->depth >= 15) {
    	image->img = XCreateImage(image->dis, CopyFromParent, image->depth, ZPixmap, 0, 
		    (char *)image->buf, image->width, image->height, 16, 0);
    }

    //XInitImage(image->img);

    // set the byte order for use with XPutImage
    if ((LSBFirst == get_byte_order())) {
    	image->img->byte_order = LSBFirst;
    } else {
	image->img->byte_order = MSBFirst;
    }

    // bitmap_bit_order doesn't matter with ZPixmap
    image->img->bitmap_bit_order = MSBFirst;

    /*
     * Most X servers use 32-bit arrays for 24-bit depth, which is easy
     * and fast to fill out.
     *
     * We assume this too, so force 32-bit data if depth is 24.
     */
    if (image->img->depth >= 24)
        image->img->bits_per_pixel = 32;

    image->img->bytes_per_line =
	image->img->width * image->img->bits_per_pixel / 8;
	
    return;

}


/*
    Removes the image without destroying the image structure.
*/		
void image_remove(image_t *image) {
    if (image) {
	if (image->img) {
    	    XDestroyImage(image->img); // should de-allocate buffer too
	    image->img = 0;
	} else if (image->buf) { // to prevent startup memory leaks
	    ckfree(image->buf);
	    image->buf = 0;
	}
    }
}


/*
   Destroy the image structure and frees all its memory.
*/
void image_destroy(image_t *image) {

    if (image) {
	image_remove(image);

	if (image->map_col) {
	    free(image->map_col);
	}
		
	free(image);
    }
}
	

