// 9/7 2D Wavelet Transform and Inverse Wavelet Transform
// Author: Han Oh, the University of Arizona 11/1/2008

#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>
#include "mex.h"
#define HORIZONTAL 0
#define VERTICAL 1

typedef unsigned char uint8;
typedef unsigned int uint32;



void forward_wavelet(double* x, double* y, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, bool bFirst);
void inverse_wavelet(double* x, double* y, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, bool bFirst);

void lazy_transform(double* x, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1);
void inverse_lazy_transform(double* x, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1);
void write_wavelet_coeff(double* x, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1);
void read_wavelet_coeff(double* x, double* y, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1);
void lifting_down(double* y0, double* y1, uint32 n_y0, uint32 n_y1, double coeff);
void lifting_up(double* y0, double* y1, uint32 n_y0, uint32 n_y1, double coeff);
void lifting_normalize(double* y0, double* y1, uint32 n_y0, uint32 n_y1, double K0, double K1);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[])
{


	double *x;
	
	int nLevel;

	if (nrhs < 1 || nrhs > 2)
	{
		mexErrMsgTxt("Invalid number of input arguments");
	}

	//number of DWT Level

	nLevel = (nrhs == 1)? 5: (uint32)(mxGetScalar(prhs[1]));

	//Getting the pointer of x

	x = mxGetPr(prhs[0]);
	uint32 nCols = (uint32)mxGetN(prhs[0]);
	uint32 nRows = (uint32)mxGetM(prhs[0]);
	uint32 nLevelCols;
	uint32 nLevelRows;

	//Allocate memory and assign output pointer
	plhs[0] = mxCreateDoubleMatrix(nCols, nRows, mxREAL);

	double* y = mxGetPr(plhs[0]);
	bool bFirst;

	//forward lifting transform
	if (nLevel > 0)
	{
		for (int i=0; i < nLevel; i++)
		{
			bFirst = (i==0)?1:0;
			nLevelRows = nRows/pow(2.0, (int)i);
			nLevelCols = nCols/pow(2.0, (int)i);
			forward_wavelet(x, y, nRows, nCols, nLevelRows, nLevelRows, bFirst);
		}
	}

	//inverse lifting transform
	if (nLevel < 0)
	{
		for (int i=-nLevel-1; i >= 0 ; i--)
		{
			bFirst = (i==-nLevel-1)?1:0;
			nLevelRows = nRows/pow(2.0, (int)i);
			nLevelCols = nCols/pow(2.0, (int)i);
			inverse_wavelet(x, y, nRows, nCols, nLevelRows, nLevelRows, bFirst);
		}
		
	}
	

}

void forward_wavelet(double* x, double* y, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, bool bFirst)
{

	double K = 1.230174105;
	double K0 = 1.0/K;
	double K1 = K/2.0;
	

	//horizontal
	uint32 n_y0 = ceil(nLevelCols/2.0);
	uint32 n_y1 = floor(nLevelCols/2.0);
	double* y0 = new double[n_y0];
	double* y1 = new double[n_y1];


	
	for (int line=0; line < nLevelRows; line++)
	{
		//lazy transform
		if (bFirst == 1)   //for the first time, read x. Otherwise, start with y
			lazy_transform(x, nRows, nCols, nLevelRows, nLevelCols, line, HORIZONTAL, y0, y1);
		else
			lazy_transform(y, nRows, nCols, nLevelRows, nLevelCols, line, HORIZONTAL, y0, y1);


		//lifting stage 1
		lifting_down(y0, y1, n_y0, n_y1, -1.586134342);

		//lifting stage 2
		lifting_up(y0, y1, n_y0, n_y1, -0.052980118);

		//lifting stage 3
		lifting_down(y0, y1, n_y0, n_y1, 0.882911075);

		//lifting stage 4
		lifting_up(y0, y1, n_y0, n_y1, 0.443506852);

		//normalization
		lifting_normalize(y0,y1, n_y0, n_y1, K0, K1);

		//write wavelet coefficients to y
		write_wavelet_coeff(y, nRows, nCols, nLevelRows, nLevelCols, line, HORIZONTAL, y0, y1);

	}

	delete[] y0;
	delete[] y1;


	//vertical
	n_y0 = ceil(nLevelRows/2.0);
	n_y1 = floor(nLevelRows/2.0);
	y0 = new double[n_y0];
	y1 = new double[n_y1];

	for (int line=0; line < nLevelCols; line++)
	{
		//lazy transform
		lazy_transform(y, nRows, nCols, nLevelRows, nLevelCols, line, VERTICAL, y0, y1);

		//lifting stage 1
		lifting_down(y0, y1, n_y0, n_y1, -1.586134342);

		//lifting stage 2
		lifting_up(y0, y1, n_y0, n_y1, -0.052980118);

		//lifting stage 3
		lifting_down(y0, y1, n_y0, n_y1, 0.882911075);

		//lifting stage 4
		lifting_up(y0, y1, n_y0, n_y1, 0.443506852);

		//normalization
		lifting_normalize(y0,y1, n_y0, n_y1, K0, K1);

		//write wavelet coefficients to y
		write_wavelet_coeff(y, nRows, nCols, nLevelRows, nLevelCols, line, VERTICAL, y0, y1);

	}

	delete[] y0;
	delete[] y1;


}

void inverse_wavelet(double* x, double* y, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, bool bFirst)
{

	double K = 1.230174105;
	double K0 = K;
	double K1 = 2.0/K;

	//horizontal
	uint32 n_y0 = ceil(nLevelCols/2.0);
	uint32 n_y1 = floor(nLevelCols/2.0);
	double* y0 = new double[n_y0];
	double* y1 = new double[n_y1];

	uint32 nLLRows, nLLCols;


	if (bFirst == 1)
	{
		nLLRows = ceil(nLevelRows/2.0);
		nLLCols = ceil(nLevelCols/2.0);
	}
	

	for (int line=0; line < nLevelRows; line++)
	{

		//read wavelet coefficients
		if (bFirst == 1 && line < nLLRows)
			memcpy(&y[line*nCols], &x[line*nCols], sizeof(double)*nLLCols);
		
		read_wavelet_coeff(x, y, nRows, nCols, nLevelRows, nLevelCols, line, HORIZONTAL, y0, y1);

		//normalization
		lifting_normalize(y0,y1, n_y0, n_y1, K0, K1);

		//lifting stage 1
		lifting_up(y0, y1, n_y0, n_y1, -0.443506852);

	
		//lifting stage 2
		lifting_down(y0, y1, n_y0, n_y1, -0.882911075);
		

		//lifting stage 3
		lifting_up(y0, y1, n_y0, n_y1, 0.052980118);
		

		//lifting stage 4
		lifting_down(y0, y1, n_y0, n_y1, 1.586134342);

		//inverse lazy transform
		inverse_lazy_transform(y, nRows, nCols, nLevelRows, nLevelCols, line, HORIZONTAL, y0, y1);
		
	}
	delete[] y0;
	delete[] y1;

	//vertical
	n_y0 = ceil(nLevelRows/2.0);
	n_y1 = floor(nLevelRows/2.0);
	y0 = new double[n_y0];
	y1 = new double[n_y1];

	for (int line=0; line < nLevelCols; line++)
	{
		//back to y
		read_wavelet_coeff(x, y, nRows, nCols, nLevelRows, nLevelCols, line, VERTICAL, y0, y1);

		//normalization
		lifting_normalize(y0,y1, n_y0, n_y1, K0, K1);

		//lifting stage 1
		lifting_up(y0, y1, n_y0, n_y1, -0.443506852);


		//lifting stage 2
		lifting_down(y0, y1, n_y0, n_y1, -0.882911075);


		//lifting stage 3
		lifting_up(y0, y1, n_y0, n_y1, 0.052980118);


		//lifting stage 4
		lifting_down(y0, y1, n_y0, n_y1, 1.586134342);

		//inverse lazy transform
		inverse_lazy_transform(y, nRows, nCols, nLevelRows, nLevelCols, line, VERTICAL, y0, y1);

	}
	delete[] y0;
	delete[] y1;
	

}


void lazy_transform(double* x, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1)
{
	
	if (dir == HORIZONTAL)
	{
		uint32 rowline = line*nCols;
		int j = 0;
		for (int i=0; i < nLevelCols; i+=2)
		{			
			y0[j] = x[rowline+i];
			y1[j] = x[rowline+(i+1)];
			j++;
		}
	}
	else //VERTICAL
	{
		int j = 0;
		for (int i=0; i < nLevelRows; i+=2)
		{
			y0[j] = x[i*nCols+line];
			y1[j] = x[(i+1)*nCols+line];
			j++;
		}		
	}
}

void inverse_lazy_transform(double* y, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1)
{
	if (dir == HORIZONTAL)
	{
		uint32 rowline = line*nCols;
		int j = 0;
		for (int i=0; i < nLevelCols; i+=2)
		{			
			y[rowline+i] = y0[j];
			y[rowline+(i+1)] = y1[j];
			j++;
		}
	}
	else //VERTICAL
	{
		int j = 0;
		for (int i=0; i < nLevelRows; i+=2)
		{
			y[i*nCols+line] = y0[j];
			y[(i+1)*nCols+line] = y1[j];
			j++;
		}		
	}

}

void write_wavelet_coeff(double* x, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1)
{
	
	int j = 0; //index for x
	if (dir == HORIZONTAL)
	{
		uint32 rowline = line*nCols;
		uint32 n_y0 = ceil(nLevelCols/2.0);
		uint32 n_y1 = floor(nLevelCols/2.0);

		//low band
		for (int i=0; i < n_y0; i++)
		{
			x[rowline+j] = y0[i];
			j++;
		}
		//high band
		for (int i=0; i < n_y1; i++)
		{
			x[rowline+j] = y1[i];
			j++;
		}		
	}
	else
	{
		uint32 n_y0 = ceil(nLevelRows/2.0);
		uint32 n_y1 = floor(nLevelRows/2.0);

		//low band
		for (int i=0; i < n_y0; i++)
		{
			x[j*nRows+line] = y0[i];
			j++;
		}
		//high band
		for (int i=0; i < n_y1; i++)
		{
			x[j*nRows+line] = y1[i];
			j++;
		}		

	}

}

void read_wavelet_coeff(double* x, double* y, uint32 nRows, uint32 nCols, uint32 nLevelRows, uint32 nLevelCols, uint32 line, bool dir, double* y0, double* y1)
{
	int j = 0; //index for x
	if (dir == HORIZONTAL)
	{
		uint32 rowline = line*nCols;
		uint32 n_y0 = ceil(nLevelCols/2.0);
		uint32 n_y1 = floor(nLevelCols/2.0);
		uint32 nLLRows = ceil(nLevelRows/2.0);

		//low band
		for (int i=0; i < n_y0; i++)
		{
			if (line < nLLRows)
				y0[i] = y[rowline+j];
			else
				y0[i] = x[rowline+j];
			j++;
		}
		//high band
		for (int i=0; i < n_y1; i++)
		{
			y1[i] = x[rowline+j];
			j++;
		}		
	}
	else
	{
		uint32 n_y0 = ceil(nLevelRows/2.0);
		uint32 n_y1 = floor(nLevelRows/2.0);

		//low band
		for (int i=0; i < n_y0; i++)
		{
			y0[i] = y[j*nCols+line];
			j++;
		}
		//high band
		for (int i=0; i < n_y1; i++)
		{
			y1[i] = y[j*nCols+line];
			j++;
		}		

	}

}

void lifting_down(double* y0, double* y1, uint32 n_y0, uint32 n_y1, double coeff)
{
	int i;
	for (i=0; i < (n_y1-1); i++)
	{
		y1[i] = y1[i] + coeff*(y0[i] + y0[i+1]);
	}
	//boundary(end) 
	y1[i] = y1[i] + coeff*(y0[i] + y0[i-1]);
}

void lifting_up(double* y0, double* y1, uint32 n_y0, uint32 n_y1, double coeff)
{
	int i = 0;
	//boundary(beginning)
	y0[i] = y0[i] + coeff*(y1[i] + y1[i+1]);
	for (i = 1; i < n_y0; i++)
	{
		y0[i] = y0[i] + coeff*(y1[i] + y1[i-1]);
	}

}

void lifting_normalize(double* y0, double* y1, uint32 n_y0, uint32 n_y1, double K0, double K1)
{
	//low band
	if (K0 != 1)
	{
		for (int i=0; i < n_y0; i++)
			y0[i] *= K0;
	}
	//high band
	if (K1 != 1)
	{
		for (int i=0; i < n_y1; i++)
			y1[i] *= K1;
	}

}


