// head file for the codes to calculate noise masked FOV block based VTs for JPEG2000 HVS based encoder, given a JND level in a 3AFC experiement
// Author: Feng Liu, Electrical and Computer Engineering, the University of Arizona, liuf@email.arizona.edu
// 3/31/2017

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double distj2k(int band, double sig2, double tb, double e, double pin);
double findpcan(double beta, int band, double sig2, double tb, double t0);
double findyb(int band, double sig2, double tb);
double t02tb(double t0, double jnd, double beta, int band, double sig2, int level, double pitch, double dist,int extrap_method);
double randn();
double randl(double sigma);
double masking(double sig2, double tb, double omega, double phro, double beta, int band, int level, double pitch, double dist,int num_samples);
double round(double x);
double trunc(double x);
void find_block_tb(double pitch, double dist, double jnd, double** tb);
void find_block_tb_given_sig2(double pitch,double dist,double jnd,double* tb,double variance,int comp, int level, int orient,int extrap_method,int num_samples,int flag_mask);
