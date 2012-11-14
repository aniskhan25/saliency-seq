#ifndef DYNAMIC_PASS_H
#define DYNAMIC_PASS_H

#include "types.h"

#include "DYN_Pyramid.h"

class DYN_DynamicPass
{

public:	
	pyramid_t pyramid;

	DYN_Pyramid oPyramid;

	DYN_DynamicPass( float *im_data_01, float *im_data_02, siz_t im_size);

	DYN_DynamicPass();
	~DYN_DynamicPass();

	float* calculate_saliency_map();

	void calculate_modulation_matrix (mod_t mod, float dpf0, float teta, siz_t im_size);

};

#endif