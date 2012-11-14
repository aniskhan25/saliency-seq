
//#include <cmath>
#include <iostream>

#include "cutil_inline.h"

#include "UTY_HelperProcs.h"

using namespace std;

UTY_HelperProcs::UTY_HelperProcs(){};
UTY_HelperProcs::~UTY_HelperProcs(){};

siz_t UTY_HelperProcs::im_size;

float UTY_HelperProcs::UTY_rint( float x)
{
	//middle value point test
	if( ceil( x+0.5)== floor( x+0.5))
	{
		int a =( int)ceil( x);
		if( a%2 == 0)
		{return ceil( x);}
		else
		{return floor( x);}
	}

	else return floor( x+0.5f);
}

int UTY_HelperProcs::UTY_rotate_90( float *in, siz_t map_size, float *out)
{
	int i, j;

	if( !in || !out) return 0;

	for( i = 0 ; i < map_size.y ; i++)
		for( j = 0 ; j < map_size.x ; j++)
			out[i*map_size.x + j] = in[j*map_size.y + i];

	return 1;
}

int UTY_HelperProcs::UTY_r2c( float *in, siz_t map_size, fftw_complex *out)
{
	int i, j;

	if( !in)return 0;

		for( i = 0; i < map_size.y*map_size.x ; i++){
		out[i][0] = in[i];
		out[i][1] = 0.0f;
	}

	return 1;
}

int UTY_HelperProcs::UTY_c2r( fftw_complex *in, siz_t map_size, float *out)
{
	int i, j;

	if( !in)return 0;

		for( i = 0; i < map_size.y*map_size.x ; i++)
		out[i] =( float)in[i][0] /( float)( map_size.y * map_size.x);

	return 1;
}

int UTY_HelperProcs::UTY_shift(fftw_complex *in ,siz_t map_size)
{	
	fftw_complex tmp;
	unsigned int centerW, centerH;

	centerW = map_size.x/2;
	centerH	= map_size.y/2;

	for( unsigned int i = 0 ; i < centerH ; i++)
	{
		for( unsigned int j = 0 ; j < map_size.x ; j++)
		{
			if (j < centerW)
			{
				tmp[0] = in[i*map_size.x + j][0];
				tmp[1] = in[i*map_size.x + j][1];

				in[i*map_size.x + j][0] = in[( i+centerH)*map_size.x + j + centerW][0];
				in[i*map_size.x + j][1] = in[( i+centerH)*map_size.x + j + centerW][1];

				in[( i+centerH)*map_size.x + j + centerW][0] = tmp[0];
				in[( i+centerH)*map_size.x + j + centerW][1] = tmp[1];
			}
			else
			{
				tmp[0] = in[i*map_size.x + j][0];
				tmp[1] = in[i*map_size.x + j][1];

				in[i*map_size.x + j][0] = in[( i+centerH)*map_size.x + j - centerW][0];
				in[i*map_size.x + j][1] = in[( i+centerH)*map_size.x + j - centerW][1];

				in[( i+centerH)*map_size.x + j - centerW][0] = tmp[0];
				in[( i+centerH)*map_size.x + j - centerW][1] = tmp[1];
			}

		}
	}

	return 1;
}

int UTY_HelperProcs::UTY_ishift(fftw_complex *in ,siz_t map_size)
{
	int i, j;
	int widthEO, heightEO ,centerW ,centerH ,xx ,yy;

	fftw_complex *temp;

	if( !in)return 0;

	temp = ( fftw_complex *)fftw_malloc ( sizeof ( fftw_complex ) * map_size.y * map_size.x );

	centerW		=	map_size.x/2;
	centerH		=	map_size.y/2;
	widthEO		= ( map_size.x%2 == 0) ? 0 : 1;
	heightEO	= ( map_size.y%2 == 0) ? 0 : 1;

	for(i = 0 ; i < map_size.y ; i++)
	{
		for(j = 0 ; j < map_size.x ; j++)
		{
			yy = (i < centerH) ? ( i+centerH + heightEO ) : ( i-centerH );
			xx = (j < centerW) ? ( j+centerW + widthEO  ) : ( j-centerW );

			temp[i*map_size.x + j][0] = in[yy*map_size.x + xx][0];
			temp[i*map_size.x + j][1] = in[yy*map_size.x + xx][1];       
		}
	}

	memcpy(in,temp, sizeof ( fftw_complex ) * map_size.y * map_size.x);

	temp = NULL;
	free(temp);

	return 1;
}

int UTY_HelperProcs::UTY_fft_2d( fftw_complex *in, fftw_complex *out, int rows, int cols)
{     
	fftw_plan plan_forward;

	if( !in)return 0;

	plan_forward = fftw_plan_dft_2d( rows, cols, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute( plan_forward);

	fftw_destroy_plan( plan_forward);

	return 1;
}

int UTY_HelperProcs::UTY_ifft_2d( fftw_complex *in, fftw_complex *out, int rows, int cols)
{  	
	fftw_plan plan_backward;

	if( !in)return 0;

	plan_backward = fftw_plan_dft_2d( rows, cols, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute( plan_backward);

	fftw_destroy_plan( plan_backward);

	fftw_cleanup_threads();

	return 1;
}

int UTY_HelperProcs::UTY_convolve_2d( float* in, float* out, int dataSizeX, int dataSizeY, 
					float* kernel, int kernelSizeX, int kernelSizeY)
{
	int i, j, m, n;
	float *inPtr, *inPtr2, *outPtr, *kPtr;
	int kCenterX, kCenterY;
	int rowMin, rowMax;                             // to check boundary of input array
	int colMin, colMax;                             //

	// check validity of params
	if( !in || !out || !kernel)return 0;
	if( dataSizeX <= 0 || kernelSizeX <= 0)return 0;

	// find center position of kernel( half of kernel size)
	kCenterX = kernelSizeX >> 1;
	kCenterY = kernelSizeY >> 1;

	// init working  pointers
	inPtr = inPtr2 = &in[dataSizeX * kCenterY + kCenterX];  // note that  it is shifted( kCenterX, kCenterY), 
	outPtr = out;
	kPtr = kernel;

	// start convolution
	for( i= 0; i < dataSizeY; ++i)                  // number of rows
	{
		// compute the range of convolution, the current row of kernel should be between these
		rowMax = i + kCenterY;
		rowMin = i - dataSizeY + kCenterY;

		for( j = 0; j < dataSizeX; ++j)             // number of columns
		{
			// compute the range of convolution, the current column of kernel should be between these
			colMax = j + kCenterX;
			colMin = j - dataSizeX + kCenterX;

			*outPtr = 0;                            // set to 0 before accumulate

			// flip the kernel and traverse all the kernel values
			// multiply each kernel value with underlying input data
			for( m = 0; m < kernelSizeY; ++m)       // kernel rows
			{
				// check if the index is out of bound of input array
				if( m <= rowMax && m > rowMin)
				{
					for( n = 0; n < kernelSizeX; ++n)
					{
						// check the boundary of array
						if( n <= colMax && n > colMin)
							*outPtr += *( inPtr - n)* *kPtr;
						++kPtr;                     // next kernel
					}
				}
				else
					kPtr += kernelSizeX;            // out of bound, move to next row of kernel

				inPtr -= dataSizeX;                 // move input data 1 raw up
			}

			kPtr = kernel;                          // reset kernel to( 0, 0)
			inPtr = ++inPtr2;                       // next input
			++outPtr;                               // next output
		}
	}

	return 1;
}

int UTY_HelperProcs::UTY_gaussian_kernel( float *h, siz_t kernel_size, float sig_p)
{
	//case 'gaussian' % Gaussian filter
	//  siz   =( p2-1)/2;
	int siz =( kernel_size.x - 1)/ 2;
	float std   = sig_p;

	//  [x, y] = meshgrid( -siz( 2):siz( 2), -siz( 1):siz( 1));
	float *X, *Y;

	float hmax;
	float temp;
	int i, j, index;	

	float sumh;

	if( !h)return 0;

	X =( float *)calloc( sizeof( float),( size_t)kernel_size.x*kernel_size.y);
	Y =( float *)calloc( sizeof( float),( size_t)kernel_size.x*kernel_size.y);

	for( i = 0; i < kernel_size.y ; i++){
		for( j = 0 ; j < kernel_size.x ; j++){           
			index = i*kernel_size.x + j;

			X[index] =( float)( j - siz);
			Y[index] =( float)( i - siz);					
		}
	}

	//  arg   = -( x.*x + y.*y)/( 2*std*std);
	//  h     = exp( arg);	
	hmax = exp( -( ( ( X[0] * X[0])+( Y[0] * Y[0]))/( 2 * std * std)));

	for( i = 0 ; i < kernel_size.y ; i++){
		for( j = 0 ; j < kernel_size.x ; j++){      
			index = i*kernel_size.x + j;

			h[index] = exp( -( ( ( X[index] * X[index])+( Y[index] * Y[index]))/( 2 * std * std)));					

			if( hmax < h[index])hmax = h[index];
		}
	}

	//  h( h<eps*max( h( :)))= 0;
	temp =( float)EPS * hmax;

	for( i = 0 ; i < kernel_size.y*kernel_size.x ; i++)
		if( h[i] < temp) h[i] = 0;

	//  sumh = sum( h( :));
	sumh = 0.0;

	
	
	for( i = 0 ; i < kernel_size.y*kernel_size.x ; i++)
			sumh += h[i];

	//  if sumh ~= 0, 
	//    h  = h/sumh;
	//  end;
	if( sumh != 0){
	for( i = 0 ; i < kernel_size.y*kernel_size.x ; i++)
				h[i] /= sumh;
	}

	return 1;
}

float UTY_HelperProcs::UTY_mean( float *a, unsigned int SIZE_W, unsigned int SIZE_H)
{
	float sum = 0.0f;

	unsigned int i;
	for (i=1 ; i<SIZE_W*SIZE_H ; i++)
		sum += a[i];

	return( sum / (float)(SIZE_W*SIZE_H))	;
}

float UTY_HelperProcs::UTY_standard_deviation( float *a, unsigned int SIZE_W, unsigned int SIZE_H)
{
	unsigned int i;
	float avg, sum01 = 0.0f, sum02 = 0.0f;

	for (i=1 ; i<SIZE_W*SIZE_H ; i++)
	{
		sum01 += a[i];
		sum02 += (a[i]*a[i]);
	}

	return ( sqrt(sum02/(SIZE_W*SIZE_H) - powf(sum01/(SIZE_W*SIZE_H), 2)));
}

float UTY_HelperProcs::UTY_skewness( float *a, unsigned int SIZE_W, unsigned int SIZE_H, float ma, float sa)
{
	unsigned int i;
	float sum = 0.0f;

	for (i=1 ; i<SIZE_W*SIZE_H ; i++)
		sum += powf((a[i] - ma),3);

	return( (float)( sum/(SIZE_W*SIZE_H)) / powf(sa,3));
}

float UTY_HelperProcs::UTY_max( float *a, unsigned int SIZE_W, unsigned int SIZE_H)
{
	unsigned int i;
	float max = a[0];
	
	for (i=1 ; i<SIZE_W*SIZE_H ; i++)
		if( a[i]> max )
			max = a[i];

	return max;
}

int UTY_HelperProcs::UTY_fusion(float *mapS, float *mapD, unsigned int SIZE_W, unsigned int SIZE_H, float *out)
{	
	unsigned int i;	
	float maxS, skwD;

	maxS = UTY_max(mapS, SIZE_W, SIZE_H);
	skwD = UTY_skewness(mapD, SIZE_W, SIZE_H, UTY_mean( mapD, SIZE_W, SIZE_H), UTY_standard_deviation( mapD, SIZE_W, SIZE_H));

	for (i=0 ; i<SIZE_W*SIZE_H ; i++)
		out[i] = maxS*mapS[i] + skwD*mapD[i] + maxS*skwD*mapS[i]*mapD[i];

	return 1;
}

/**************************************************************************
Compensation
**************************************************************************/

#define DATA_LINE_LEN 16
#define LINE_MAX 256

int UTY_HelperProcs::iof_get_compensation (const char *filename, float list[], unsigned int frameNumber)
{
	int n;

	FILE *file = fopen ( filename, "r" );
	if ( file != NULL ){
		char line[LINE_MAX] ;

		while( fgets( line, sizeof line, file ) != NULL )
		{
			if (atoi(&line[0]) == frameNumber)
			{	
				char *p;

				p = strtok (line," ");

				if (p != NULL)
				{
					sscanf(p, "%d", &n);
					//printf("Frame No. = %d\n", n);

					for (unsigned int i=0 ; i<DATA_LINE_LEN ; i++)
					{
						p = strtok (NULL, " ");
						sscanf(p, "%f ", &list[i]);

						//printf("list[%d] = %f\n", i, list[i]);
					}
				}
			}
		}

		fclose ( file );
	}
	else {
		perror ( filename );
	}

	return 0;
}

int UTY_HelperProcs::BilinearInterpolation(float *in, float *vx, float *vy, int w, int h, float *out)
{		
	float fraction_x, fraction_y, one_minus_x, one_minus_y;
	int ceil_x, ceil_y, floor_x, floor_y;

	float pix[4];

	for (int y = 0; y < h ; ++y)
		for (int x = 0; x < w ; ++x)
		{
			if (vx[y*w+x]<0.0f || vx[y*w+x]>w-1 || vy[y*w+x]<0.0f || vy[y*w+x]>h-1) { out[y*w+x] = 0.0f; continue; }

			floor_x = (int)floor(vx[y*w+x]);			
			floor_y = (int)floor(vy[y*w+x]);

			if (floor_x < 0) floor_x = 0;
			if (floor_y < 0) floor_y = 0;

			ceil_x = floor_x + 1;
			if (ceil_x >= w) ceil_x = floor_x;

			ceil_y = floor_y + 1;
			if (ceil_y >= h) ceil_y = floor_y;
			
			fraction_x = vx[y*w+x] - (float)floor_x;
			fraction_y = vy[y*w+x] - (float)floor_y;
			
			one_minus_x = 1.0 - fraction_x;
			one_minus_y = 1.0 - fraction_y;

			pix[0] = in[floor_y*w + floor_x];
			pix[1] = in[floor_y*w + ceil_x];
			pix[2] = in[ceil_y*w + floor_x];
			pix[3] = in[ceil_y*w + ceil_x];

			out[y*w + x] = one_minus_y * 
				(one_minus_x * pix[0] + fraction_x * pix[1]) + fraction_y * (one_minus_x * pix[2] + fraction_x * pix[3]);
		}

		return 0;
}