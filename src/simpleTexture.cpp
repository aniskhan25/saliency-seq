////
////#include <cstdlib>
////#include <iostream>
////
////#include "cutil_inline.h"
////
////#include "types.h"
////#include "DYN_DynamicPass.h"
////
////char *image_filename001 = "lena_bw.pgm";
////char *image_filename002 = "lena_bw.pgm";
////
////////////////////////////////////////////////////////////////////////////////////
////// declaration, forward
////void runTest( int argc, char** argv);
////
////extern "C"
////void computeGold( float* reference, float* idata, const unsigned int len);
////
////////////////////////////////////////////////////////////////////////////////////
////// Program main
////////////////////////////////////////////////////////////////////////////////////
////int
////main( int argc, char** argv) 
////{
////    runTest( argc, argv);
////
////    cutilExit(argc, argv);
////}
////
////
////////////////////////////////////////////////////////////////////////////////////
//////! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////////
////void
////runTest( int argc, char** argv) 
////{	
////    // load image from disk
////    float* h_data001 = NULL;
////	float* h_data002 = NULL;
////
////    unsigned int width, height;
////
////	char* image_path;
////
////	image_path = cutFindFilePath(image_filename001, argv[0]);
////    if (image_path == 0) {
////        printf("Unable to source file file %s\n", image_filename001);
////        exit(EXIT_FAILURE);
////    }
////    cutilCheckError( cutLoadPGMf(image_path, &h_data001, &width, &height));
////
////	image_path = cutFindFilePath(image_filename002, argv[0]);
////    if (image_path == 0) {
////        printf("Unable to source file file %s\n", image_filename002);
////        exit(EXIT_FAILURE);
////    }
////    cutilCheckError( cutLoadPGMf(image_path, &h_data002, &width, &height));
////
////    unsigned int size = width * height * sizeof(float);
////    printf("Loaded '%s' & '%s', %d x %d pixels\n", image_filename001, image_filename002, width, height);
//// 
////	siz_t im_size;
////	im_size.x = width;
////	im_size.y = height;
////
////	DYN_DynamicPass oDynamicPass( h_data001, h_data002, im_size);
////
////	float* h_odata = (float*) malloc( size);
////
////    cutFree(image_path);
////}
