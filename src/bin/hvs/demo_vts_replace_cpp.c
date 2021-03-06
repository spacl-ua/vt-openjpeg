// demo for noise masked FOV block based VT calculations
// Author: Feng Liu, Electrical and Computer Engineering, the University of Arizona, liuf@email.arizona.edu
// 3/31/2017
//Subbands appear in the sequence, LL_D, HL_D, LH_D, ...,HL_1, LH_1, HH_1,
#include<stdio.h>
#include<math.h>
#include "covertfunctions.h"
void main(int argc, char *argv[])
{
	double jnd=-1.0;
	int level=-1;
	int orient=-1;
	int comp=-1;
	double variance=-1.0;
	int extrap_method=-1;
	
	int ic;
	for (ic = 1; ic<argc ; ic++){
		if (sscanf(argv[ic],"jnd=%lf",&jnd)==1){
//			fprintf(stderr,"JND=%lf\n",jnd);
		}
		if (sscanf(argv[ic],"variance=%lf",&variance)==1){
//			fprintf(stderr,"variance=%lf\n",variance);
		}
		/*
		if (sscanf(argv[ic],"comp=%d",&comp)==1){
//			fprintf(stderr,"comp=%d\n",comp);
		}
		if (sscanf(argv[ic],"level=%d",&level)==1){
//			fprintf(stderr,"level=%d\n",level);
		}
		if (sscanf(argv[ic],"orient=%d",&orient)==1){
//			fprintf(stderr,"orient=%d\n",orient);
		}
		*/
		if (sscanf(argv[ic],"extrap_method=%d",&extrap_method)==1){
//			fprintf(stderr,"extrap_method=%d\n",extrap_method);
		}
	}
	if (level<0 || orient< 0 || comp<0 || jnd <0){
		//fprintf(stderr,"missing cmd options\n");
		//return;	
	}
	int i=0;
	/*
	double** tb= (double**)malloc(15* sizeof(double*));
	for(i=0;i<15;i++)
		tb[i]= (double*)malloc(3 * sizeof(double));
	*/
	double tb[1];
	double pitch = 0.1845;	// Monitor PA328Q
	double dist = 600.0; // 60cm obervation distance
	/*find_block_tb(pitch, dist, jnd, tb);*/

	for (comp=0; comp<3; comp++){
		switch(comp){
			case 0:
				printf("double b_b1[5][2]={");
				break;

			case 1:
				printf("double b1_c1[5][2]={");
				break;

			case 2:
				printf(" double b1_c2[5][2]={");
				break;
		}
		for (level=0; level<=4; level++){
			for (orient=3; orient>=1; orient--){
			/* LL */
			if (level==4 && orient==0)
				continue;
			if (orient==2)
				continue;
			find_block_tb_given_sig2(pitch, dist, jnd, tb,variance,comp,level,orient,extrap_method,500);
//			printf("The resulted block based VTs are:\n");
			if (orient == 3)
				printf("{%f,",tb[0]);
			else 
				printf("%f}",tb[0]);
			fflush(stdout);
			} /*orient*/

			if (level!=4)
				printf(",");
		}/*level*/
		printf("};\n");
	} /*comp*/
/* replcae block_encoder.cpp*/
printf("//REPLACE the block_encoder.cpp\n");
printf("double ll5_jnd_thresholds[3]={");
	for (comp=0; comp<3; comp++){
		level=4;orient=0;
		find_block_tb_given_sig2(pitch, dist, jnd, tb,variance,comp,level,orient,extrap_method,500);
		if (comp!=2)
			printf("%f,",tb[0]);
		else
			printf("%f};\n",tb[0]);
	}
	/*
	printf("Y\tCb\tCr\n");
	for (i = 0; i < 15; i++)
		printf("%f\t%f\t\%f\n",tb[i][0],tb[i][1],tb[i][2]);
	for (i = 0; i<15; i++)
		free(tb[i]);
	free(tb);
	*/
//	system("pause");
//
}
