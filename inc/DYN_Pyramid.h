#ifndef PYRAMID_H
#define PYRAMID_H

#include "types.h"

class DYN_Pyramid
{

public:
	pyramid_t *pyramid;
	
	DYN_Pyramid();

	~DYN_Pyramid();

	void init( pyramid_t *_pyramid, float *_im_data01, float *_im_data02, siz_t _im_size);

	void init_pyramid();

	void create_pyramid();

	void create_pyramid_precedent();

	void reinit_pyramid_precedent();

	void calculate_modulation_matrix (mod_t out, float f0, siz_t im_size);

	void gabor_filter (float* im_data, std::complex<float> *b, mod_t mod, float* temp);

	void free_pyramid();

};

#endif