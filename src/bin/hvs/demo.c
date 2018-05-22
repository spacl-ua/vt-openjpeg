// demo for noise masked FOV block based VT calculations
// Author: Feng Liu, Electrical and Computer Engineering, the University of Arizona, liuf@email.arizona.edu
// 3/31/2017

#include<stdio.h>
#include<math.h>
#include "covertfunctions.h"
void main(int argc, char *argv[])
{
	double jnd;
	if (sscanf(argv[1],"%lf",&jnd)!=1){
		fprintf(stderr,"[ERROR] Provide jnd !\n");
		return; 
	}
	int i=0;
	double** tb= (double**)malloc(15* sizeof(double*));
	for(i=0;i<15;i++)
		tb[i]= (double*)malloc(3 * sizeof(double));
	double pitch = 0.1845;	// Monitor PA328Q
	double dist = 600.0; // 60cm obervation distance
	find_block_tb(pitch, dist, jnd, tb);
	printf("The resulted block based VTs are:\n");
	printf("Y\tCb\tCr\n");
	for (i = 0; i < 15; i++)
		printf("%f\t%f\t\%f\n",tb[i][0],tb[i][1],tb[i][2]);
	for (i = 0; i<15; i++)
		free(tb[i]);
	free(tb);
//	system("pause");
}
