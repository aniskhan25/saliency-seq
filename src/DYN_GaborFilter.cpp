
#include "DYN_GaussianRecursive.h"

#include "DYN_GaborFilter.h"

DYN_GaborFilter::DYN_GaborFilter(){}
DYN_GaborFilter::~DYN_GaborFilter(){}

int DYN_GaborFilter::apply_gabor_filtering( gabor_mask_t gabor_mask, float* level_data, mod_t mod, siz_t level_size, siz_t mod_size)
{
	DYN_GaborFilter::high_pass_prefilter( gabor_mask.masks[N-1], level_data, level_size, S1);		

	for( int n=0 ; n<N ; n++)
	{
		DYN_GaborFilter::modulation( gabor_mask.masks[n], gabor_mask.masks[N-1], mod.mt_mod[n], level_size, mod_size, 2.0, S1);

		DYN_GaussianRecursive::gaussian_recursive (gabor_mask.masks[n].re, gabor_mask.masks[n].re, level_size.x, level_size.y, S1);
		DYN_GaussianRecursive::gaussian_recursive (gabor_mask.masks[n].im, gabor_mask.masks[n].im, level_size.x, level_size.y, S1);

		DYN_GaborFilter::demodulation( gabor_mask.masks[n], gabor_mask.masks[n], mod.mt_mod[n], level_size, mod_size, 2.0, S1);
	}

	return 0;
}

void DYN_GaborFilter::projection( motion_vect_t v_curr, motion_vect_t v_prev, siz_t v_curr_size, siz_t v_prev_size) 
{
	int i,j,m,l;

#pragma omp for nowait
	for( j=0 ; j<v_curr_size.y ; j++){
		for( i=0 ; i<v_curr_size.x ; i++){

			m =  i + j*v_curr_size.x;
			l = (i/2) + (j/2)*v_prev_size.x;

			v_curr.vx[m] = v_prev.vx[l] * 2;
			v_curr.vy[m] = v_prev.vy[l] * 2;			

			//pyr->dec[k].lum[m]=0;
			//temp[m]=0;
		}
	}
}

void DYN_GaborFilter::interpolation( float *dec, motion_vect_t v_proj, float *niv, siz_t level_size)
{
	int i,j,m,l;
	int a,b,sa,sb;
	float u,v,val;
	float k00,k10,k01,k11;
	float i00,i10,i01,i11;
	float j00,j10,j01,j11;

	float *temp = (float *) malloc( level_size.x*level_size.y * sizeof(float));

	#pragma omp for nowait
	for( j=0 ; j<level_size.y ; j++){
		for( i=0 ; i<level_size.x ; i++){

			m = i + j*level_size.x;

			val = niv[m];

			u = - v_proj.vx[m];
			v = - v_proj.vy[m];

			a = (int)u;
			b = (int)v;

			if( u>=0) sa=1; else sa=-1;
			if( v>=0) sb=1; else sb=-1;

			u = fabs( u - (float)a);
			v = fabs( v - (float)b);

			i00 = i + (float)a; i01 = i00; i10 = i00 + (float)sa; i11 = i10;
			j00 = j + (float)b; j10 = j00; j01 = j00 + (float)sb; j11 = j01;

			k00 = (1-u) * (1-v);
			k10 =	 u  * (1-v);
			k01 = (1-u) *	 v;
			k11 =	 u  *	 v;

			if( i00>=0 && i00<level_size.x && j00>=0 && j00<level_size.y){
				l = (int)(i00 + j00*level_size.x);

				dec [l] += k00*val;
				temp[l] += k00;
			}
			if( i10>=0 && i10<level_size.x && j10>=0 && j10<level_size.y){
				l = (int)(i10 + j10*level_size.x);

				dec [l] += k10*val;
				temp[l] += k10;
			}
			if( i01>=0 && i01<level_size.x && j01>=0 && j01<level_size.y){
				l = (int)(i01 + j01*level_size.x);

				dec [l] += k01*val;
				temp[l] += k01;
			}
			if( i11>=0 && i11<level_size.x && j11>=0 && j11<level_size.y){
				l = (int)(i11 + j11*level_size.x);

				dec [l] += k11*val;
				temp[l] += k11;
			}
		}
	}

	#pragma omp for nowait
	for( m=0 ; m<level_size.x*level_size.y ; m++)
		( temp[m]!=0) ? dec[m] /= temp[m] : dec[m] = 127.0f;		
}

void DYN_GaborFilter::modulation( cmplx_t out, cmplx_t prev, cmplx_t mod, siz_t level_size, siz_t mod_size, float dpf0, float teta)
{
	int i, j, l, m;

	/* Modulation */
#pragma omp parallel shared(out,prev,mod,level_size,mod_size) private(j,i,m,l)
	{
#pragma omp for nowait
		for( j=0 ; j<level_size.y ; j++){
			for( i=0 ; i<level_size.x ; i++){

				m = i + j*level_size.x;
				l = i + j*mod_size.x;

				out.re[m] = prev.re[m] * mod.re[l];
				out.im[m] = prev.re[m] * mod.im[l];
			}
		}
	}
}

void DYN_GaborFilter::demodulation( cmplx_t out, cmplx_t in, cmplx_t mod, siz_t level_size, siz_t mod_size, float dpf0, float teta)
{
	int i, j, l, m;
	float norm;

	/* Demodulation */
#pragma omp parallel shared(out,in,mod,level_size,mod_size) private(j,i,m,l,norm)
	{
#pragma omp for nowait
		for( j=0 ; j<level_size.y ; j++){
			for( i=0 ; i<level_size.x ; i++){

				m = i + j*level_size.x;
				l = i + j*mod_size.x;

				norm = sqrtf( powf( in.re[m], 2) + powf( in.im[m], 2)); if( norm==0) ++norm;

				// (a, ib)*(c, id) = [(ac+bd), i(-ad+bc)] / sqrt(a^2 + b^2)
				out.re[m] = (	 in.re[m]*mod.re[l] + in.im[m]*mod.im[l]) / norm;
				out.im[m] = ( -( in.re[m]*mod.im[l] - in.im[m]*mod.re[l])) / norm;
			}
		}
	}
}

int DYN_GaborFilter::high_pass_prefilter( cmplx_t out, float *in, siz_t level_size, float sigma)
{	
	int m;

	DYN_GaussianRecursive::gaussian_recursive( out.re, in, level_size.x, level_size.y, sigma);

#pragma omp parallel shared(out,in,level_size) private(m)
	{
#pragma omp for nowait
		for( int m=0 ; m<level_size.x*level_size.y ; ++m)	
			out.re[m] = in[m] - out.re[m];
	}
	return 1;
}