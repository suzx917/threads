#include "bitmap.h"

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>

#define MAX_THREAD_COUNT 128

double xcenter = 0;
double ycenter = 0;
double scale = 4;
int    image_width = 500;
int    image_height = 500;
int    max = 1000;
int    thread_count = 1;

struct bitmap *bm;

int iteration_to_color( int i, int max );
int iterations_at_point( double x, double y, int max );
void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int tid );

void show_help()
{
	printf("Use: mandel [options]\n");
	printf("Where options are:\n");
	printf("-m <max>    The maximum number of iterations per point. (default=1000)\n");
	printf("-x <coord>  X coordinate of image center point. (default=0)\n");
	printf("-y <coord>  Y coordinate of image center point. (default=0)\n");
	printf("-s <scale>  Scale of the image in Mandlebrot coordinates. (default=4)\n");
	printf("-W <pixels> Width of the image in pixels. (default=500)\n");
	printf("-H <pixels> Height of the image in pixels. (default=500)\n");
	printf("-o <file>   Set output file. (default=mandel.bmp)\n");
	printf("-h          Show this help text.\n");
	printf("\nSome examples are:\n");
	printf("mandel -x -0.5 -y -0.5 -s 0.2\n");
	printf("mandel -x -.38 -y -.665 -s .05 -m 100\n");
	printf("mandel -x 0.286932 -y 0.014287 -s .0005 -m 1000\n\n");
}

void *start_routine(void *arg) 
{
    int t = *((int*)arg);
    //printf("\\\\executing start routine #%d\n", t);
    // Compute the Mandelbrot image
	compute_image(bm,xcenter-scale,xcenter+scale,ycenter-scale,ycenter+scale,max,t);
    pthread_exit(NULL);
}

int main( int argc, char *argv[] )
{
	char c;


	// These are the default configuration values used
	// if no command line arguments are given.

	const char *outfile = "mandel.bmp";

	// For each command line argument given,
	// override the appropriate configuration value.

	while((c = getopt(argc,argv,"x:y:s:W:H:n:m:o:h"))!=-1) {
		switch(c) {
			case 'x':
				xcenter = atof(optarg);
				break;
			case 'y':
				ycenter = atof(optarg);
				break;
			case 's':
				scale = atof(optarg);
				break;
			case 'W':
				image_width = atoi(optarg);
				break;
			case 'H':
				image_height = atoi(optarg);
				break;
            case 'n':
                thread_count = atoi(optarg);
                //printf("thread count set as %d\n", thread_count);
				break;
			case 'm':
				max = atoi(optarg);
				break;
			case 'o':
				outfile = optarg;
				break;
			case 'h':
				show_help();
				exit(1);
                break;
		}
	}

    bm = bitmap_create(image_width,image_height);

	// Display the configuration of the image.
	printf("mandel: x=%lf y=%lf scale=%lf max=%d outfile=%s\n",xcenter,ycenter,scale,max,outfile);

	// Fill it with a dark blue, for debugging
	bitmap_reset(bm,MAKE_RGBA(0,0,255,0));
    
///////////////////////////////////////////////////////////////////////////////////////////////////
    struct timeval start, stop;

    // *** Multithreads ***
    pthread_t ptt[MAX_THREAD_COUNT];
    int* temp = (int*)malloc(thread_count*sizeof(int));

    // start timer
    gettimeofday(&start,NULL);

    for (int i = 0; i < thread_count; ++i) {
        temp[i] = i;
        if (pthread_create(&ptt[i], NULL, start_routine, (void*)&(temp[i]) )) 
        {
            perror("Error creating thread: ");
            exit( EXIT_FAILURE ); 
        }
        // else
        // {
        //     printf("\\\\created tid = %d\n", ptt[i]);
        // }
        
    }

    for (int i = 0; i < thread_count; ++i) {
        if (pthread_join(ptt[i], NULL))
        {
            perror("Error creating thread: ");
            exit( EXIT_FAILURE ); 
        }        
    }
    // stop timer
    gettimeofday(&stop,NULL);
    double sec = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
    
    // save result to log file
    FILE* log = fopen("log.csv", "a");
    fprintf(log, "%d, %.4f\n", thread_count, sec);

	// Save the image in the stated file.
	if(!bitmap_save(bm,outfile)) {
		fprintf(stderr,"mandel: couldn't write to %s: %s\n",outfile,strerror(errno));
		return 1;
	}
	return 0;
}

/*
Compute an entire Mandelbrot image, writing each point to the given bitmap.
Scale the image to the range (xmin-xmax,ymin-ymax), limiting iterations to "max"
*/

void compute_image( struct bitmap *bm, double xmin, double xmax, double ymin, double ymax, int max, int t )
{
    //printf("\\\\ computing_image t = %d, tcount = %d\n", t, thread_count);
	int i,j;

	int width = bitmap_width(bm);
	int height = bitmap_height(bm);

	// For every pixel in the image...

	for(j= t*height/thread_count ;j < (t+1)*height/thread_count ; j++) {

		for(i=0;i<width;i++) {

			// Determine the point in x,y space for that pixel.
			double x = xmin + i*(xmax-xmin)/width;
			double y = ymin + j*(ymax-ymin)/height;

			// Compute the iterations at that point.
			int iters = iterations_at_point(x,y,max);

			// Set the pixel in the bitmap.
			bitmap_set(bm,i,j,iters);
		}
	}
}

/*
Return the number of iterations at point x, y
in the Mandelbrot space, up to a maximum of max.
*/

int iterations_at_point( double x, double y, int max )
{
	double x0 = x;
	double y0 = y;

	int iter = 0;

	while( (x*x + y*y <= 4) && iter < max ) {

		double xt = x*x - y*y + x0;
		double yt = 2*x*y + y0;

		x = xt;
		y = yt;

		iter++;
	}

	return iteration_to_color(iter,max);
}

/*
Convert a iteration number to an RGBA color.
Here, we just scale to gray with a maximum of imax.
Modify this function to make more interesting colors.
*/

int iteration_to_color( int i, int max )
{
	int gray = 255*i/max;
	return MAKE_RGBA(gray,gray,gray,0);
}
