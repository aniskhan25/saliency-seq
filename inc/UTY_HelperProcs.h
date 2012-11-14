
#ifndef UTILITY_H
#define UTILITY_H

#include "types.h"
#include "fftw3.h"

class UTY_HelperProcs
{
private:
	UTY_HelperProcs();
	~UTY_HelperProcs();

public:		
	static float UTY_rint( float x);
	
	static int UTY_rotate_90( float *in, siz_t map_size, float *out);
	
	static int UTY_r2c( float *in, siz_t map_size, fftw_complex *out);
	
	static int UTY_c2r( fftw_complex *in, siz_t map_size, float *out);
	
	static int UTY_shift(fftw_complex *in ,siz_t map_size);
	
	static int UTY_ishift(fftw_complex *in ,siz_t map_size);
	
	static int UTY_fft_2d( fftw_complex *in, fftw_complex *out, int rows, int cols);
	
	static int UTY_ifft_2d( fftw_complex *in, fftw_complex *out, int rows, int cols);
	
	static int UTY_convolve_2d( float* in, float* out, int dataSizeX, int dataSizeY
		, float* kernel, int kernelSizeX, int kernelSizeY);
	
	static int UTY_gaussian_kernel( float *h, siz_t kernel_size, float sig_p);
	
	static siz_t UTY_imread( const char *name, float *h_idata);

	static int UTY_imwrite( float *in, siz_t im_size, char *image_path);

	static siz_t UTY_HelperProcs::UTY_imread( float *h_idata);

	static void UTY_iof_save_pgm( const char *name, float *mat, siz_t im_size);
	
	static float* UTY_iof_load_pgm( const char *name, siz_t im_size);

	static siz_t UTY_iof_read_size( const char *name);

	static siz_t im_size;

	static float UTY_mean( float *a, unsigned int SIZE_W, unsigned int SIZE_H);

	static float UTY_standard_deviation( float *a, unsigned int SIZE_W, unsigned int SIZE_H);

	static float UTY_skewness( float *a, unsigned int SIZE_W, unsigned int SIZE_H, float ma, float sa);

	static float UTY_max( float *a, unsigned int SIZE_W, unsigned int SIZE_H);

	static int UTY_fusion(float *mapS, float *mapD, unsigned int SIZE_W, unsigned int SIZE_H, float *out);

	static int iof_get_compensation (const char *filename, float list[], unsigned int frameNumber);

	static int BilinearInterpolation(float *in, float *vx, float *vy, int w, int h, float *out);
};

#endif