#include <omp.h>
#include <time.h>
#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    int i, nthreads;
    clock_t clock_timer;
    double wall_timer;
    for (nthreads = 1; nthreads <=8; nthreads++) {
        clock_timer = clock();
        wall_timer = omp_get_wtime();
        #pragma omp parallel for private(i) num_threads(nthreads)
        for (i = 0; i < 100000000; i++) cos(i*1.0f);
        printf("%d threads: time on clock() = %.3f, on wall = %.3f\n", \
            nthreads, \
            (double) (clock() - clock_timer) / CLOCKS_PER_SEC, \
            omp_get_wtime() - wall_timer);
    }
}

/*

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector_types.h>
#include <omp.h>

#include "cutil_inline.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "types.h"

#include "UTY_CommonProcs.h"
#include "UTY_HelperProcs.h"
#include "STA_StaticPass.h"

#include "DYN_DynamicPass.h"

#define NUM_THREADS 4

using namespace std;

void temporal_filtering( float *out, float** in, int w, int h );

#define QUEUE_LEN 13

inline void exit(){	cout << endl << "Press enter to continue..." << endl; cin.get();}

int main( int argc, char** argv) 
{	
	omp_set_num_threads( omp_get_max_threads());

	float mvt_lst[16];

	cv::Mat	im001, im002;

	float *h_idata, *h_idata_01, *h_idata_02, *h_idata_com
		, *h_odata_sta, *h_odata_dyn, *h_odata;

	float **maps = NULL;

	std::vector<int> params;

	params.push_back( CV_IMWRITE_PXM_BINARY );
	params.push_back( 100 );

	int width, height;

	im001 = cv::imread( "C:\\Sophie_these\\Program_Marat_21-01-10\\saillance\\data\\original\\pgm\\ClipRMX_3_0000i.pgm", 0 );

	siz_t im_size;
	im_size.x = im001.cols;
	im_size.y = im001.rows;

	h_idata	   = (float *)calloc( im_size.x*im_size.y, sizeof(float) );	
	h_idata_01 = (float *)calloc( im_size.x*im_size.y, sizeof(float) );	
	h_idata_02 = (float *)calloc( im_size.x*im_size.y, sizeof(float) );	
	h_idata_com = (float *)calloc( im_size.x*im_size.y, sizeof(float) );
	h_odata	   = (float *)calloc( im_size.x*im_size.y, sizeof(float) );	

	maps = (float**) malloc ( QUEUE_LEN * sizeof(float*) );
	for( unsigned int i=0;i<QUEUE_LEN ;++i )
		maps[i] = (float*) malloc ( im_size.x*im_size.y * sizeof(float) );

	int NO_OF_IMAGES[] = { 
		0, //0
		0, //1
		0, //2
		2, //745,
		743,
		752,
		775,
		761,
		735,
		0, //9
		765,
		747,
		739,
		775,
		738,
		725,
		758,
		739,
		784,
		773,
		791,
		727,
		778,
		770,
		734,
		735
	};

	char filename[100];
	char command[100];

	unsigned int queue_idx, queue_full;

	int id;

	std::cout << "Starting calculating saliency maps ..." << std::endl;

	clock_t start = clock();

	clock_t clock_timer;
    double wall_timer;

	for( unsigned int k=0;k<26;++k )
	{
		switch( k )
		{
		case 3:
			/*		case 4:
			case 5:
			case 6:
			case 7:
			case 8:
			case 10:
			case 11:
			case 12:
			case 13:
			case 14:
			case 15:
			case 16:
			case 17:
			case 18:
			case 19:
			case 20:
			case 21:
			case 22:
			case 23:
			case 24:
			case 25:
			*/
			{
				std::cout << "Snippet " << k << std::endl;

				sprintf( filename, "%s%d%s", "C:\\Sophie_these\\Program_Marat_21-01-10\\saillance\\data\\original\\pgm\\ClipRMX_", k, "_0000i.pgm" );
				im001 = cv::imread( filename, 0 );

				queue_full = 0;
				queue_idx = 0;

				for( unsigned int i=1;i<NO_OF_IMAGES[k];++i ){

					std::cout << "Image " << i << "\r" ;

					sprintf( filename, "%s%d%s%s%d%s", "C:\\Sophie_these\\Program_Marat_21-01-10\\saillance\\data\\original\\pgm\\ClipRMX_", k, "_", (i<10)?"000":(i<100)?"00":(i<1000)?"0":"", i, "i.pgm" );
					im002 = cv::imread( filename, 0 );		

#pragma omp parallel shared(im001,im002,h_idata,h_idata_01,h_idata_02,im_size) private(id)
					{
#pragma omp for nowait
						for ( id=0 ; id<im_size.x*im_size.y ; ++id)
						{
							h_idata[id]	   = im002.data[id];

							h_idata_01[id] = im001.data[id];
							h_idata_02[id] = im002.data[id];
						}
					}  /* end of parallel section */

					sprintf( filename, "%s%d%s", "C:\\Sophie_these\\Program_Marat_21-01-10\\saillance\\data\\ClipRMX_", k, ".txt" );
					UTY_HelperProcs::iof_get_compensation( filename, mvt_lst, i);
					UTY_CommonProcs::MotionCompensation( h_idata_01, im_size.x, im_size.y, mvt_lst, h_idata_com);

					clock_timer = clock();
					wall_timer = omp_get_wtime();

					STA_StaticPass oStaticPass( im_size);
					h_odata_sta = oStaticPass.calculate_saliency_map( h_idata);

					cout << "(STA) on clock = " << ( (double)clock() - clock_timer) / CLOCKS_PER_SEC << ", on wall = " << omp_get_wtime() - wall_timer << endl;

					clock_timer = clock();
					wall_timer = omp_get_wtime();

					DYN_DynamicPass oDynamicPass( h_idata_01, h_idata_02, im_size);
					h_odata_dyn = oDynamicPass.calculate_saliency_map();

					cout << "(DYN) on clock = " << ( (double)clock() - clock_timer) / CLOCKS_PER_SEC << ", on wall = " << omp_get_wtime() - wall_timer << endl;

					float max_s = std::numeric_limits<float>::min();
					float max_d = std::numeric_limits<float>::min();

#pragma omp parallel
					{
						float t_max_s = std::numeric_limits<float>::min();
						float t_max_d = std::numeric_limits<float>::min();
#pragma omp for nowait
						for (int id = 0;id<im_size.x*im_size.y;++id)
						{	
							t_max_s = std::max(h_odata_sta[id], t_max_s);
							t_max_d = std::max(h_odata_dyn[id], t_max_d);
						}
#pragma omp critical 
						{
							max_s = std::max(max_s, t_max_s);
							max_d = std::max(max_d, t_max_d);
						}
					}  /* end of parallel section */

#pragma omp parallel shared(h_odata_sta,h_odata_dyn,max_s,max_d,im_size) private(id)
					{
#pragma omp for nowait
						for (id = 0;id<im_size.x*im_size.y;++id){

							h_odata_sta[id] = h_odata_sta[id]/max_s * 255.0f;
							h_odata_dyn[id] = h_odata_dyn[id]/max_d * 255.0f;
						}
					}  /* end of parallel section */

					UTY_HelperProcs::UTY_fusion( h_odata_sta, h_odata_dyn, im_size.x, im_size.y, h_odata );

					float max_sd = std::numeric_limits<float>::min();

#pragma omp parallel 
					{
						float t_max_sd = std::numeric_limits<float>::min();						
#pragma omp for nowait
						for (int id = 0;id<im_size.x*im_size.y;++id)
						{	
							t_max_sd = std::max(h_odata[id], t_max_sd);							
						}
#pragma omp critical 
						{
							max_sd = std::max(max_sd, t_max_sd);							
						}
					}  /* end of parallel section */

					if( !queue_full ){

#pragma omp parallel shared(maps,queue_idx,max_sd,im001,im002,im_size) private(id)
						{
#pragma omp for nowait
							for (id = 0;id<im_size.x*im_size.y;++id)
								maps[queue_idx][id] = h_odata[id]/max_sd * 255.0f;						

#pragma omp for nowait
							for (id = 0;id<im_size.x*im_size.y;++id){

								im001.data[id] = im002.data[id];
								im002.data[id] = (char)(unsigned char)( maps[queue_idx][id] );
							}
						}  /* end of parallel section */

						++queue_idx;

						if( queue_idx == QUEUE_LEN ){
							queue_full = 1;
							queue_idx = 0;
						}
					}
					else{

						if( queue_idx == QUEUE_LEN ){
							queue_idx = 0;
						}

#pragma omp parallel shared(maps,queue_idx,max_sd,im_size) private(id)
						{
#pragma omp for nowait
							for (id = 0;id<im_size.x*im_size.y;++id)			
								maps[queue_idx][id] = h_odata[id]/max_sd * 255.0f;						
						}  /* end of parallel section */

						temporal_filtering( h_odata, maps, im_size.x, im_size.y );


#pragma omp parallel shared(h_odata,im001,im002,im_size) private(id)
						{
#pragma omp for nowait
							for (id = 0;id<im_size.x*im_size.y;++id){

								im001.data[id] = im002.data[id];
								im002.data[id] = (char)(unsigned char)( h_odata[id] );
							}
						}  /* end of parallel section */

						++queue_idx;
					}

					sprintf( filename, "%s%d%s%s%d%s", "D:\\tmp\\dynamic\\ClipRMX_", k, "_", (i<10)?"000":(i<100)?"00":(i<1000)?"0":"", i, "_fus.pgm" );
					cv::imwrite( filename, im002, params );
				}

				std::cout << "Calculation of saliency maps completed." << std::endl;
			}
		}
	}

	cout << "Time elapsed: " << ( (float)clock() - start) / CLOCKS_PER_SEC << endl;

	for( unsigned int i=0;i<QUEUE_LEN;++i )
		free(maps[i]);
	free(maps);

	h_idata = h_odata = h_odata_sta = h_odata_dyn = h_idata_01 = h_idata_02 = NULL;

	free( h_idata);
	free( h_odata);
	free( h_odata_sta);		
	free( h_odata_dyn);		
	free( h_idata_01);
	free( h_idata_02);

	exit();
}

void temporal_filtering( float *out, float** in, int w, int h )
{
	float arr[QUEUE_LEN];

	float tmp;

	unsigned int min_idx;

	int middle = QUEUE_LEN/2;

	float average;

	for( unsigned int idx=0;idx<w*h;++idx )
	{
		for( unsigned int i=0;i<QUEUE_LEN;++i )
			arr[i] = in[i][idx];

		for( unsigned int i=0;i<QUEUE_LEN-1;++i ){

			min_idx = i;
			for( unsigned int j=i+1;j<QUEUE_LEN;++j ){

				if( arr[min_idx]>arr[j] )
					min_idx = j;
			}

			tmp = arr[i];
			arr[i] = arr[min_idx];
			arr[min_idx] = tmp;
		}

		if( QUEUE_LEN%2==0 ) average = (float)(arr[middle-1]+arr[middle])/2;
		else average = arr[middle];

		out[idx] = average;
	}
}
*/