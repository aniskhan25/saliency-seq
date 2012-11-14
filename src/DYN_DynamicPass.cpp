
#include "DYN_GaborFilter.h"
#include "DYN_Gradients.h"
#include "DYN_MotionEstimation.h"
#include "DYN_GaussianRecursive.h"

#include "DYN_DynamicPass.h"

#include "STA_Pathway.h"

DYN_DynamicPass::DYN_DynamicPass( float *im_data_01, float *im_data_02, siz_t im_size)
{
	oPyramid.init( &pyramid, im_data_01, im_data_02, im_size);
}

DYN_DynamicPass::DYN_DynamicPass(){}
DYN_DynamicPass::~DYN_DynamicPass(){}

float* DYN_DynamicPass::calculate_saliency_map()
{			
	STA_Pathway oStaticPathway( pyramid.im_size);

	for (unsigned int id = 0;id<pyramid.im_size.x*pyramid.im_size.y;++id)
		pyramid.seuil_tmp[id] = 100 * ( ( ( fabs( pyramid.im_data01[id]-pyramid.im_data02[id] ) ) < 0.2f ) ? 1.0f : 0.0f );

	for( int k=K-1 ; k>=0 ; k--){

		oStaticPathway.retina_filter( pyramid.pre_levels[k].im_data, pyramid.levels[k].level_size, 0, 0, 0, pyramid.pre_levels[k].im_data);
		oStaticPathway.retina_filter( pyramid.levels[k].im_data, pyramid.levels[k].level_size, 0, 0, 0, pyramid.levels[k].im_data);

		if( k<K-1){
			DYN_GaborFilter::projection( pyramid.levels[k].motion_vect, pyramid.levels[k+1].motion_vect, pyramid.levels[k].level_size, pyramid.levels[k+1].level_size);			

			DYN_GaborFilter::interpolation( pyramid.levels[k].dec, pyramid.levels[k].motion_vect, pyramid.levels[k].im_data, pyramid.levels[k].level_size);

			DYN_GaborFilter::apply_gabor_filtering( pyramid.pre_levels[k].gabor_mask, pyramid.pre_levels[k].im_data, pyramid.mod, pyramid.pre_levels[k].level_size, pyramid.im_size);
			DYN_GaborFilter::apply_gabor_filtering( pyramid.levels[k].gabor_mask,	  pyramid.levels[k].im_data,	 pyramid.mod, pyramid.levels[k].level_size, pyramid.im_size);
		}
		else 
			DYN_GaborFilter::apply_gabor_filtering( pyramid.levels[k].gabor_mask, pyramid.levels[k].im_data, pyramid.mod, pyramid.levels[k].level_size, pyramid.im_size);

		//iof_save_txt_simple ("mask.txt", pyramid.levels[k].gabor_mask.masks[0].re, pyramid.levels[k].level_size.x, pyramid.levels[k].level_size.y);

		DYN_Gradients::calculate_gradients( pyramid.grad, pyramid.levels[k].im_data, pyramid.levels[k].gabor_mask, pyramid.pre_levels[k].gabor_mask, pyramid.mod, pyramid.levels[k].level_size);			

		//iof_save_txt_simple ("Gx.txt", pyramid.grad.gx[0].re, pyramid.levels[k].level_size.x, pyramid.levels[k].level_size.y);
		//iof_save_txt_simple ("Gy.txt", pyramid.grad.gy[0].re, pyramid.levels[k].level_size.x, pyramid.levels[k].level_size.y);
		//iof_save_txt_simple ("Gt.tx", pyramid.grad.gt[0].re, pyramid.levels[k].level_size.x, pyramid.levels[k].level_size.y);

		if (k==K-1) 
			DYN_MotionEstimation::estimate_motion( pyramid.levels[k].motion_vect, pyramid.grad, pyramid.levels[k].level_size, true);
		else
			DYN_MotionEstimation::estimate_motion( pyramid.levels[k].motion_vect, pyramid.grad, pyramid.levels[k].level_size, false);

		DYN_GaussianRecursive::gaussian_recursive( pyramid.levels[k].motion_vect.vx, pyramid.levels[k].motion_vect.vx, pyramid.levels[k].level_size.x, pyramid.levels[k].level_size.y, S2);
		DYN_GaussianRecursive::gaussian_recursive( pyramid.levels[k].motion_vect.vy, pyramid.levels[k].motion_vect.vy, pyramid.levels[k].level_size.x, pyramid.levels[k].level_size.y, S2);
	}

	//iof_save_txt_simple ("Vx.txt", pyramid.levels[0].motion_vect.vx, pyramid.im_size.x, pyramid.im_size.y);
	//iof_save_txt_simple ("Vy.txt", pyramid.levels[0].motion_vect.vy, pyramid.im_size.x, pyramid.im_size.y);

	float *out = ( float *)malloc( pyramid.im_size.x*pyramid.im_size.y * sizeof(float));

	for (unsigned int id = 0;id<pyramid.im_size.x*pyramid.im_size.y;++id){

		if (pyramid.seuil_tmp[id]<=VP)
			out[id] = sqrt( 
			powf( pyramid.levels[0].motion_vect.vx[id], 2) + 
			powf( pyramid.levels[0].motion_vect.vy[id], 2)
			);
		else
			out[id] = 0.0f;		
	}
	
	return out;
}
