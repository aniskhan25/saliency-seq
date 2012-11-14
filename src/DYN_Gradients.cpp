
#include "DYN_Gradients.h"

int DYN_Gradients::calculate_gradients( gradient_t grad, float* level_data, gabor_mask_t mask, gabor_mask_t mask_prev, mod_t mod, siz_t level_size)
{
	#pragma omp for nowait
	for( int n=0 ; n<N ; ++n)
	{
		DYN_Gradients::calculate_gradients_x_y ( grad.gx[n], grad.gy[n], mask.masks[n], level_size);
		DYN_Gradients::calculate_gradients_t( grad.gt[n], mask.masks[n], mask_prev.masks[n], level_size);
	}

	return 0;
}

void DYN_Gradients::calculate_gradients_x_y ( cmplx_t gx, cmplx_t gy, cmplx_t in, siz_t level_size) 
{
	int i, j, idx;	

	#pragma omp for nowait
	for( j=0 ; j<level_size.y ; j++){
		for( i=0 ; i<level_size.x ; i++){

			idx = i + level_size.x*j;

			// default constructor sets to (0,0)
			gx.re[idx] = gx.im[idx] = gy.re[idx] = gy.im[idx] = 0.0f;			

			if ( !( i<2 || i>=level_size.x-2 || j<2 || j>=level_size.y-2))
			{
				gx.re[idx] = ( in.re[idx-2] - 8.0f*in.re[idx-1] + 8.0f*in.re[idx+1] - in.re[idx+2]) / 12.0f;
				gx.im[idx] = ( in.im[idx-2] - 8.0f*in.im[idx-1] + 8.0f*in.im[idx+1] - in.im[idx+2]) / 12.0f;

				gy.re[idx] = ( in.re[idx-2*level_size.x] - 8.0f*in.re[idx-level_size.x] + 8.0f*in.re[idx+level_size.x] - in.re[idx+2*level_size.x]) / 12.0f;
				gy.im[idx] = ( in.im[idx-2*level_size.x] - 8.0f*in.im[idx-level_size.x] + 8.0f*in.im[idx+level_size.x] - in.im[idx+2*level_size.x]) / 12.0f;
			}
		}
	}
}

void DYN_Gradients::calculate_gradients_t ( cmplx_t gt, cmplx_t curr, cmplx_t prev, siz_t level_size) 
{
	#pragma omp for nowait
	for(int m=0 ; m<level_size.x*level_size.y ; ++m){
		gt.re[m] = curr.re[m] - prev.re[m];
		gt.im[m] = curr.im[m] - prev.im[m];
	}
}