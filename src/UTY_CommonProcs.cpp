
#include "UTY_HelperProcs.h"
#include "UTY_CommonProcs.h"

//#include <cstdlib>
//#include <cmath>

UTY_CommonProcs::UTY_CommonProcs(){};
UTY_CommonProcs::~UTY_CommonProcs(){};

int UTY_CommonProcs::UTY_retina_gaussion_filter( float *kernel, float *im
							   , siz_t kernel_size, siz_t im_size
							   , float sigma, float *out)
{
	float *tmp;

	if( !kernel || !im || !out) return 0;

	tmp = ( float *)malloc( sizeof( float)* kernel_size.x*kernel_size.y);

	UTY_HelperProcs::UTY_gaussian_kernel( kernel, kernel_size, sigma);

	UTY_HelperProcs::UTY_rotate_90( kernel, kernel_size, tmp);
	UTY_HelperProcs::UTY_rotate_90( kernel, kernel_size, tmp);

	//if( stencilSize.y*stencilSize.x > im_size.y*im_size.x)   
	UTY_HelperProcs::UTY_convolve_2d( im, out, im_size.x, im_size.y, kernel, kernel_size.x, kernel_size.y);

	free( tmp);

	return 1;
}

int UTY_CommonProcs::UTY_normalization_nl( float *map, siz_t map_size)
{   
	int i;

	float map_min, map_max;

	if( !map)return 0;

	map_min = map_max = map[0];

	for( i = 0 ; i < map_size.y*map_size.x ; i++){

		if( map_min > map[i])map_min = map[i];
		if( map_max < map[i])map_max = map[i];
	}

	for( i = 0 ; i < map_size.y*map_size.x ; i++){

		if( map_min != map_max){
			map[i] -= map_min;
			map[i] /=( map_max - map_min);
		}

		if( LOWER_BOUND > map[i])map[i] = LOWER_BOUND;
		else if( UPPER_BOUND < map[i])map[i] = UPPER_BOUND;
		else{				
			// TODO :: Initially no effect				
			//image[i*map_size.x + j] -= LOWER_BOUND;
			//image[i*map_size.x + j] /=( UPPER_BOUND - LOWER_BOUND);
		}	
	}

	return 1;
}

int UTY_CommonProcs::UTY_normalization_itti( float *map, siz_t map_size)
{
	int i;

	float sum, _max, temp, mean;

	if( !map)return 0;

	sum = 0.0; _max = map[0];

	for( i = 0 ; i < map_size.y*map_size.x ; i++){
		sum += map[i];

		if( _max < map[i])
			_max = map[i];			
	}

	mean = sum /( map_size.y*map_size.x);

	temp =( _max - mean)*( _max - mean);

	for( i = 0 ; i < map_size.y*map_size.x ; i++)
		map[i] *= temp;

	return 1;
}

int UTY_CommonProcs::UTY_normalization_pc( float *map, siz_t map_size)
{
	int i;

	float level, _max;

	if( !map)return 0;

	_max = map[0];
	for( i = 0 ; i < map_size.y*map_size.x ; i++)		
		if( _max < map[i])_max = map[i];

	level = PC_LEVEL * _max;

	for( i = 0 ; i < map_size.y*map_size.x ; i++)
		if( map[i] > level)	map[i] =( map[i] - level)/( _max - level)* _max;			
		else				map[i] = 0.0f;

		return 1;
}

int UTY_CommonProcs::UTY_frequency_to_spatial( fftw_complex *io, siz_t map_size)
{
	if( !io)return 0;

	UTY_HelperProcs::UTY_shift( io, map_size);

	UTY_HelperProcs::UTY_ifft_2d( io, io, map_size.y, map_size.x); 

	return 1;
}

int UTY_CommonProcs::UTY_spatial_to_frequency( float *in, siz_t map_size, fftw_complex *out)
{     
	if( !in)return 0;

	UTY_HelperProcs::UTY_r2c( in, map_size, out);

	UTY_HelperProcs::UTY_fft_2d( out, out, map_size.y, map_size.x);

	UTY_HelperProcs::UTY_shift( out, map_size);

	return 1;
}

int	UTY_CommonProcs::UTY_retina_bipolar( float *in_photo, float *in_hor, siz_t map_size, int beta, float *out)
{
	int i;

	float tempON, tempOFF;

	if( !in_photo || !in_hor)return 0;

	for( i = 0 ; i < map_size.y*map_size.x ; i++){	
		tempON = in_photo[i] - beta * in_hor[i];

		tempOFF = beta * in_hor[i] - in_photo[i];                   

		out[i] =( ( tempON < 0.0f)? 0.0f : tempON)-( ( tempOFF < 0.0f)? 0.0f : tempOFF);    
	}    

	return 1;
}


int UTY_CommonProcs::UTY_mask_create( float *mask, siz_t map_size, int mu)
{    
	float var1 =( float)floor( map_size.x / 2.0f);
	float var2 = pow( var1, mu);

	float var3 =( float)floor( map_size.y / 2.0f);
	float var4 = pow( var3, mu);

	int i, j;

	if( !mask)return 0;

	for( i = 0 ; i < map_size.y ; i++){
		for( j = 0 ; j < map_size.x ; j++){
			mask[i*map_size.x + j] =( 
				( 1.0f -( pow( ( j - var1), mu)/ var2))*
				( 1.0f -( pow( ( i - var3), mu)/ var4))
				);
		}
	}

	return 1;
}

int UTY_CommonProcs::UTY_mesh_create( float *in, float d1, float d2, int n)
{   
	float diff = d2 - d1;

	int i;

	if( !in)return 0;

	for( i = 0; i < n - 1 ; i++)
	{            
		in[i] = i*diff /( floor( ( float)n)- 1)+ d1;
	}

	in[i] = d2;

	return 1;
}

int UTY_CommonProcs::UTY_gabor_masks( float *u1, float *v1, double *teta, siz_t map_size)
{
	double theta;

	int j, y, x;

	float *u, *v;

	if( !u1 || !v1)return 0;

	u = ( float *)malloc( sizeof( float)*( map_size.x + 1));
	v = ( float *)malloc( sizeof( float)*( map_size.y + 1));

	UTY_mesh_create( u, -0.5f,  0.5f, map_size.x + 1);
	UTY_mesh_create( v,  0.5f, -0.5f, map_size.y + 1);

	for( j = 0; j < NO_OF_ORIENTS ; j++){

		theta = teta[j] * PI / 180.0;

		for( y = 0 ; y < map_size.y ; y++){
			for( x = 0 ; x < map_size.x ; x++){

				u1[j *( map_size.y*map_size.x)+( y*map_size.x + x)] = 
					u[x]*( float)cos( theta)+ v[y]*( float)sin( theta);
				v1[j *( map_size.y*map_size.x)+( y*map_size.x + x)] = 
					v[y]*( float)cos( theta)- u[x]*( float)sin( theta);				
			}
		}	
	}

	free( u);
	free( v);

	return 1;
}

int UTY_CommonProcs::MotionCompensation(float *in, int w, int h, float *lsMC, float *out)
{
	int nh = 2*lsMC[0], nw = 2*lsMC[1];

	float *vx,*vy;
	vx = (float*)malloc(nw*nh * sizeof(float));
	vy = (float*)malloc(nw*nh * sizeof(float));

	float t1, t2;
	for (int i=0 ; i<nh ; ++i)
	{
		for (int j=0 ; j<nw ; ++j)
		{
			t1 = (float)(j - lsMC[1]);
			t2 = (float)(i - lsMC[0]);

			vx[i*nw+j] = lsMC[2] + lsMC[4]*t1 + lsMC[5]*t2 + lsMC[8 ]*(t1*t1) + lsMC[9 ]*(t1*t2) + lsMC[10]*(t2*t2) + j;
			vy[i*nw+j] = lsMC[3] + lsMC[6]*t1 + lsMC[7]*t2 + lsMC[11]*(t1*t1) + lsMC[12]*(t1*t2) + lsMC[13]*(t2*t2) + i;
		}
	}
	
	UTY_HelperProcs::BilinearInterpolation(in, vx, vy, nw, nh, out);

	return 0;
}
