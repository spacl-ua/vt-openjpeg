// codes to calculate noise masked FOV block based VTs for JPEG2000 HVS based encoder, given a JND level in a 3AFC experiement
// Author: Feng Liu, Electrical and Computer Engineering, the University of Arizona, liuf@email.arizona.edu
// codes to calculate noise masked FOV block based VTs for JPEG2000 HVS based encoder, given a JND level in a 3AFC experiement
// // Author: Feng Liu, Electrical and Computer Engineering, the University of Arizona, liuf@email.arizona.edu
// 3/31/2017

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double distj2k(int band, double sig2, double tb, double e, double pin);	// function to calculate the distortion PDF
// band=0: LL, band=1: HL/LH/HH; sig2: wavelet coefficent variance; tb: block based threhsold; e: quantization distortion; pin: a precalculated variable to accelerate the calculation
double findpcan(double beta, int band, double sig2, double tb, double t0); // function to find the probability of distortion detection in a 3AFC FOV block based experiment
// beta=3.2 for 3AFC; band=0: LL, band=1: HL/LH/HH; sig2: wavelet coefficient variance; tb: block based VT; t0: single coefficient threshold leading to 75% probability of detection
double findyb(int band, double sig2, double tb); // function to calculate the mean of the absolute value of the quantizaion distortion
// band=0: LL, band=1: HL/LH/HH; sig2: wavelet coefficent variance; tb: block based VT
double t02tb(double t0, double jnd, double beta, int band, double sig2, int level, double pitch, double dist,int extrap_method); // function to calculate the theoretical block based VT (no noise masking)
// t0: single coefficient threshold leading to 75% probability of detection; jnd: JND level; beta=3.2 for 3AFC; band=0: LL, band=1: HL/LH/HH; sig2: wavelet coefficient variance; level: wavelet decomposition level; pitch: monitor pitch size; dist: distance from the subject to the monitor
double randn(); // normal random
double randl(double sigma); // Laplacian random with zero mean
// sigma: standard deviation
double masking(double sig2, double tb, double omega, double phro, double beta, int band, int level, double pitch, double dist,int num_samples); // masking factor calculation for Luminance Component
// sig2: wavelet coefficent variance; tb: block based VT; omega, phro, beta: masking parameters; band=0: LL, band=1: HL/LH/HH; level: wavelet decomposition level; pitch: monitor pitch size; dist: distance from the subject to the monitor
double round(double x); // rounding
double trunc(double x); // truncation to the nearest integer toward 0
void find_block_tb(double pitch, double dist, double jnd, double** tb); // main function to using the aove function to calculate masked noise masked block based VTs given any JND level input
// pitch: monitor pitch size; dist: distance from the subject to the monitor; jnd: JND level; tb: address to put the resulted VTs, tb[15][3] is organized in the form of:
/*
YHH1		CbHH1		CrHH1
YHL/LH1		CbHL/LH1	CrHL/LH1
YLL1		CbLL1		CrLL1
YHH2		CbHH2		CrHH2
YHL/LH2		CbHL/LH2	CrHL/LH2
YLL2		CbLL2		CrLL2
YHH3		CbHH3		CrHH3
YHL/LH3		CbHL/LH3	CrHL/LH3
YLL3		CbLL3		CrLL3
YHH4		CbHH4		CrHH4
YHL/LH4		CbHL/LH4	CrHL/LH4
YLL4		CbLL4		CrLL4
YHH5		CbHH5		CrHH5
YHL/LH5		CbHL/LH5	CrHL/LH5
YLL5		CbLL5		CrLL5
*/
void find_block_tb_given_sig2(double pitch,double dist,double jnd,double* tb,double variance,int comp, int level, int orient,int extrap_method,int num_samples,int flag_mask);
/*calculating VT using the provided variance instead of a fixed variance
 * Variance: the variance of the code block. If variance=-1.0, then use the default variance. 
 * level: range [0,4]
 * orient:{0,1,2,3}	
 * comp:{0,1,2}
 * flag_mask: 1 - turn masking on ; 0 - turn masking off
 */


double distj2k(int band, double sig2, double tb, double e, double pin)
{
	double z=0.0;
	if (band==0)
	{
		if(fabs(e)<=tb/2)
			z=1/sqrt(12*sig2)+(1-pin)/tb;
		else if((fabs(e)>tb/2)&&(fabs(e)<=tb))
			z=1/sqrt(12*sig2);
		else
			z=0.0;
	}
	else
	{
		if(fabs(e)<=tb/2)
			z=exp(-sqrt(2/sig2)*fabs(e))/sqrt(2*sig2)+(1-pin)/tb;
		else if((fabs(e)>tb/2)&&(fabs(e)<=tb))
			z=exp(-sqrt(2/sig2)*fabs(e))/sqrt(2*sig2);
		else
			z=0.0;
	}
	return z;
}

double findpcan(double beta, int band, double sig2, double tb, double t0)
{
	int n=128;
	int i=0;
	double pin=0.0;
	if(band==0)
		pin=tb/sqrt(3*sig2);
	else
		pin=1-exp(-sqrt(2/sig2)*tb);
	double pcan1=0.0;
	double pcan2=0.0;
	double update_err=100.0;
	double stepsize=2.0*tb/((double)(n));
	for(i=0;i<n;i++)
	{
		pcan1+=distj2k(band,sig2,tb,(i+0.5)*stepsize-tb,pin)*exp(-pow(fabs(((i+0.5)*stepsize-tb)/t0),beta));
	}
	pcan1=pcan1*stepsize;
	while(update_err>1.0e-7)
	{
		pcan2=pcan1;
		n*=2;
		stepsize/=2.0;
		pcan1=0.0;
		for(i=0;i<n;i++)
		{
			pcan1+=distj2k(band,sig2,tb,(i+0.5)*stepsize-tb,pin)*exp(-pow(fabs(((i+0.5)*stepsize-tb)/t0),beta));
		}
		pcan1=pcan1*stepsize;
		update_err=fabs(pcan1-pcan2);
	}
	return pcan1;
}

double findyb(int band, double sig2, double tb)
{
	int n=128;
	int i=0;
	double pin=0.0;
	if(band==0)
		pin=tb/sqrt(3*sig2);
	else
		pin=1-exp(-sqrt(2/sig2)*tb);
	double yb1=0.0;
	double yb2=0.0;
	double update_err=100.0;
	double stepsize=2.0*tb/((double)(n));
	for(i=0;i<n;i++)
	{
		yb1+=distj2k(band,sig2,tb,(i+0.5)*stepsize-tb,pin)*fabs((i+0.5)*stepsize-tb);
	}
	yb1=yb1*stepsize;
	while(update_err>1.0e-7)
	{
		yb2=yb1;
		n*=2;
		stepsize/=2.0;
		yb1=0.0;
		for(i=0;i<n;i++)
		{
			yb1+=distj2k(band,sig2,tb,(i+0.5)*stepsize-tb,pin)*fabs((i+0.5)*stepsize-tb);
		}
		yb1=yb1*stepsize;
		update_err=fabs(yb1-yb2);
	}
	return yb1;
}

double t02tb(double t0, double jnd, double beta, int band, double sig2, int level, double pitch, double dist,int extrap_method)
{
/*Extrap method
 * 	0: JND is defined as subband JND
 * 	3: JND is defined as image-level JND(need to address non 5-level decomposition(fractional JNDs allowed) 
 * 	2: same as 1, but if subband JND is lower than 1, 1 is used to avoid fractional subband JND*/

	//total number of coefficients 
	int n_total=(int)round((2*dist*tan(3.1415926/180.0))/pitch);
	//number of coefficients in the subband
	int n=(int)round((2*dist*tan(3.1415926/180.0))/(pitch*pow(2.0,(double)level)));
	//probability of no-detection for each coeffcients(uniform for all coeffs.)
	double ptrans_unit = exp((log(1.5) + jnd*log(0.25)) / (n_total*n_total));
	double log_ptrans_unit = (log(1.5) + jnd*log(0.25)) / (n_total*n_total);
	double ptrans;
	//prob. of no-detection at subband_JND=1
	double ptrans_sub_jnd1 = exp((log(1.5) + log(0.25)) / (n*n)); 

        double ptrans_def1 = exp((log(1.5) + jnd*log(0.25)) / (n*n)); //using subband jnd defininition

	double ptrans_def2 = ptrans_unit;
//exp(log_p_trans_unit*((double)(n_total)/pow(2.0,(double)(level)))); //using image-jnd def.

	switch (extrap_method){
		case 0:
			ptrans = ptrans_def1;
			break;
		case 2:
			ptrans= (ptrans_def2>ptrans_sub_jnd1)?ptrans_sub_jnd1:ptrans_def2;
			break;
		case 3:
			ptrans = ptrans_def2;
			break;
		default:
			fprintf(stderr,"extrap_method should only be from {0,2,3}\n");
			exit(0);
		
		
	}


	double tbcan1=t0*0.01;
	double pcan1=findpcan(beta,band,sig2,tbcan1,t0);
	double tbcan2=t0*100;
	double pcan2=findpcan(beta,band,sig2,tbcan2,t0);
	while((pcan2>ptrans)||(pcan1<ptrans))
	{
		if(pcan2>ptrans)
		{
			tbcan2=tbcan2*2;
			pcan2=findpcan(beta,band,sig2,tbcan2,t0);
		}
		if(pcan1<ptrans)
		{
			tbcan1=tbcan1/2;
			pcan1=findpcan(beta,band,sig2,tbcan1,t0);
		}
	}
	double tbcanm=0;
	double pcanm=0;
	while(tbcan2-tbcan1>((double)(1.0e-7)))
	{
		tbcanm=(tbcan1+tbcan2)/2.0;
		pcanm=findpcan(beta,band,sig2,tbcanm,t0);
		if(pcanm>ptrans)
			tbcan1=tbcanm;
		else
			tbcan2=tbcanm;
	}
	return (tbcan1+tbcan2)/2.0;
}

double randn()
{
	static double V1,V2,S;
	static int phase=0;

	if(phase==1)
	{
		phase=1-phase;
		return V2*sqrt(-2*log(S)/S);
	}
	else
	{
		do
		{
			double U1=(double)rand()/RAND_MAX;
			double U2=(double)rand()/RAND_MAX;
			V1=2.0*U1-1.0;
			V2=2.0*U2-1.0;
			S=V1*V1+V2*V2;
		}
		while((S>=1)||(S==0));
		phase=1-phase;
		return V1*sqrt(-2*log(S)/S);
	}
}

double randl(double sigma)
{
	double x;
	do
	{
		x=(double)rand()/RAND_MAX;
	}
	while((x==1)||(x==0));
	return -sigma*(2.0*round(x)-1.0)*log(1.0-2.0*fabs(x-0.5))/sqrt(2.0);
}

double masking(double sig2, double tb, double omega, double phro, double beta, int band, int level, double pitch, double dist,int num_samples)
{
	if (omega!=0.0)
	{
		int n = (int)round((2 * dist*tan(3.1415926 / 180.0)) / (pitch*pow(2.0, (double)level)));
		int i = 0;
		int j = 0;
		double mb = 0.0;
		double mb_sum = 0.0;
		double tmpval = 0.0;
		time_t t;
		srand((unsigned)time(&t));
		double* randarray = (double*)malloc(n*n * sizeof(double));
		double yb = findyb(band, sig2, tb);
		for (i = 0; i < num_samples; i++)
		{
			if (band == 0)
			{
				mb = 0.0;
				for (j = 0; j < n*n; j++)
				{
					tmpval = sqrt(sig2)*randn();
					if (fabs(trunc(tmpval / tb)) < 0.1)
						randarray[j] = tmpval;
					else if (trunc(tmpval / tb) > 0.0)
						randarray[j] = tmpval - tb*(trunc(tmpval / tb) + 0.5);
					else
						randarray[j] = tmpval - tb*(trunc(tmpval / tb) - 0.5);
					randarray[j] = omega*pow((fabs(randarray[j]) / yb), phro);
					if (randarray[j] < 1.0)
						randarray[j] = 1.0;
					mb += pow(randarray[j], beta);
				}
				mb = pow(mb / ((double)(n*n)), 1 / beta);
			}
			else
			{
				mb = 0.0;
				for (j = 0; j < n*n; j++)
				{
					tmpval = randl(sqrt(sig2));
					if (fabs(trunc(tmpval / tb)) < 0.1)
						randarray[j] = tmpval;
					else if (trunc(tmpval / tb) > 0.0)
						randarray[j] = tmpval - tb*(trunc(tmpval / tb) + 0.5);
					else
						randarray[j] = tmpval - tb*(trunc(tmpval / tb) - 0.5);
					randarray[j] = omega*pow((fabs(randarray[j]) / yb), phro);
					if (randarray[j] < 1.0)
						randarray[j] = 1.0;
					mb += pow(randarray[j], beta);
				}
				mb = pow(mb / ((double)(n*n)), 1 / beta);
			}
			mb_sum += mb;
		}
		free(randarray);
		return tb*mb_sum / ((double)num_samples);
	}
	else
		return 128.0;
}

double round(double x)
{
	if(x>=0)
		return (double)((int)(x+0.5));
	else
		return (double)((int)(x-0.5));
}

double trunc(double x)
{
	return (double)((int)x);
}


void find_block_tb_given_sig2(double pitch,double dist,double jnd,double* tb,double variance,int comp, int _level, int orient,int extrap_method,int num_samples,int flag_mask){
	/*fprintf(stderr,"find_block_threshold(pitch=%lf,dist=%lf,jnd=%lf,tb[0]=%lf,variance=%lf,comp=%d,_level=%d,orient=%d)......\n",pitch,dist,jnd,tb[0],variance,comp,_level,orient);*/
	/*level: range [0,4]*/
	/*orient:{0,1,2,3}*/	
	/*comp:{0,1,2}*/
	int i = 0;
	i = _level*3 + 
		((orient>=2) ? (3-orient) : (2-orient));
	
	double t0[45] = { 15.1589500000000,5.78200300000000,3.86422600000000,2.27105000000000,1.48838200000000,1.44031800000000,0.881776000000000,0.629355000000000,0.547459000000000,0.536842000000000,0.413154000000000,0.491147000000000,0.342918000000000,0.304238000000000,0.361628000000000,78.5283700000000,37.0890900000000,33.6068000000000,12.2486400000000,5.63615200000000,12.9292700000000,5.13597100000000,4.22272500000000,6.33303100000000,2.85380200000000,2.27571000000000,2.50087000000000,1.59646400000000,1.01583900000000,1.34284100000000,19.0255900000000,8.38992600000000,7.82262700000000,4.84437300000000,2.69327000000000,3.25015900000000,2.19162700000000,1.45638500000000,1.60518700000000,1.02453000000000,0.777205000000000,0.763466000000000,0.564523000000000,0.452921000000000,0.418489000000000 };
	double sig2[45] = { 2.16099832608020,8.17260263615171,1530.09098781364,9.41378984356092,42.3703693274298,1375.18351640656,39.2712161904762,122.172799357143,980.864081044705,82.3596455736434,162.263606674419,890.691695096900,70.1153483372094,114.123851445736,589.763500317829,0.646588069940939,0.907185713630710,22.4146445174075,0.547895114952463,1.09062621447084,19.0811568916703,0.705490948493683,1.72557702186589,13.4564817074830,1.09502434883721,2.12036664922481,12.6044835775194,0.923182294573643,1.51636455038760,8.57970437596899,1.42834087674853,2.59700577758783,134.706178802145,2.14184157346587,5.41104890539960,125.386408244713,4.05179737609329,9.88127590281827,97.7506646705540,5.56551618604651,11.2003352093023,96.8042756046511,5.20223722480620,9.13866604457364,73.2788081124031 };
	double omega[15] = { 0.0,2.51188643150958,1.58489319246111,0.0398107170553497,0.398107170553497,0.00158489319246111,0.000630957344480193,2.51188643150958,0.630957344480194,1.0,1.58489319246111,0.00158489319246111,1.58489319246111,0.158489319246111,0.251188643150958 };
	double phro[15] = { 0.0,0.630957344480194,0.0100000000000000,3.98107170553497,1.58489319246111,15.8489319246111,15.8489319246111,0.251188643150958,3.98107170553497,2.51188643150958,0.251188643150958,15.8489319246111,0.00630957344480193,3.98107170553497,6.30957344480193 };
	double beta[15] = { 0.0,0.0158489319246111,25.1188643150958,6.30957344480193,15.8489319246111,0.000398107170553497,0.00100000000000000,0.158489319246111,1.0,0.0630957344480194,630.957344480193,0.0158489319246111,0.630957344480194,1.58489319246111,0.00251188643150958 };
	double scl_cb[15] = { 0.0,0.0,0.987443835660405,0.0,0.0,0.802171644410405,0.0,2.07309221817751,0.762900884574164,1.32734416516483,0.925111062832940,0.797787522347863,0.610967428488595,0.928270215536784,0.905828222074527 };
	double scl_cr[15] = { 0.0,3.45498410289289,1.87825091785739,3.47599894937224,2.56390718416733,1.57674397937559,1.83407623664399,1.91030648245327,1.39603770369435,1.57449905506255,1.56770562115591,1.37504912388072,1.23015550287436,1.29569497078972,1.54105745038586 };	// parameters from measurements
	int level[15] = { 1,1,1,2,2,2,3,3,3,4,4,4,5,5,5 };
	int band[15] = { 1,1,0,1,1,0,1,1,0,1,1,0,1,1,0 };

	variance = (variance<0) ? (sig2[i+comp*15]) : (variance) ;
	switch (comp){
	case 0: /*Luma*/ 
	{
		if (omega[i] != 0)
		{
			tb[0] = t02tb(t0[i], jnd, 3.2, band[i], variance, level[i], pitch, dist,extrap_method);
			if (flag_mask)
			tb[0] = masking(variance, tb[0], omega[i], phro[i], beta[i], band[i], level[i], pitch, dist,num_samples);
		}
		else
			tb[0] = 128.0;	// infinity, the same case in the following 
		//fprintf(stderr,"tb[0]=%lf\n",tb[0]);
	}
	break;

	case 1: /*Cb*/
	{
		if (scl_cb[i] != 0.0)
		{
			tb[0] = t02tb(t0[15 + i], jnd, 3.2, band[i], variance, level[i], pitch, dist,extrap_method);
			if (flag_mask)
			tb[0] = (scl_cb[i])*(tb[0]);
		}
		else
			tb[0] = 128.0;
		//fprintf(stderr,"tb[0]=%lf\n",tb[0]);
	}
	break;

	case 2: /*Cr*/
	{
		if (scl_cr[i] != 0.0)
		{
			tb[0] = t02tb(t0[30 + i], jnd, 3.2, band[i], variance, level[i], pitch, dist,extrap_method);
			if (flag_mask)
			tb[0] = (scl_cr[i])*(tb[0]);
		}
		else
			tb[0] = 128.0;
		//fprintf(stderr,"tb[0]=%lf\n",tb[0]);
	}
	break;
	}/*switch(comp)*/
}
	

void find_block_tb(double pitch,double dist,double jnd,double** tb)
{
	int i = 0;
	double t0[45] = { 15.1589500000000,5.78200300000000,3.86422600000000,2.27105000000000,1.48838200000000,1.44031800000000,0.881776000000000,0.629355000000000,0.547459000000000,0.536842000000000,0.413154000000000,0.491147000000000,0.342918000000000,0.304238000000000,0.361628000000000,78.5283700000000,37.0890900000000,33.6068000000000,12.2486400000000,5.63615200000000,12.9292700000000,5.13597100000000,4.22272500000000,6.33303100000000,2.85380200000000,2.27571000000000,2.50087000000000,1.59646400000000,1.01583900000000,1.34284100000000,19.0255900000000,8.38992600000000,7.82262700000000,4.84437300000000,2.69327000000000,3.25015900000000,2.19162700000000,1.45638500000000,1.60518700000000,1.02453000000000,0.777205000000000,0.763466000000000,0.564523000000000,0.452921000000000,0.418489000000000 };
	double sig2[45] = { 2.16099832608020,8.17260263615171,1530.09098781364,9.41378984356092,42.3703693274298,1375.18351640656,39.2712161904762,122.172799357143,980.864081044705,82.3596455736434,162.263606674419,890.691695096900,70.1153483372094,114.123851445736,589.763500317829,0.646588069940939,0.907185713630710,22.4146445174075,0.547895114952463,1.09062621447084,19.0811568916703,0.705490948493683,1.72557702186589,13.4564817074830,1.09502434883721,2.12036664922481,12.6044835775194,0.923182294573643,1.51636455038760,8.57970437596899,1.42834087674853,2.59700577758783,134.706178802145,2.14184157346587,5.41104890539960,125.386408244713,4.05179737609329,9.88127590281827,97.7506646705540,5.56551618604651,11.2003352093023,96.8042756046511,5.20223722480620,9.13866604457364,73.2788081124031 };
	double omega[15] = { 0.0,2.51188643150958,1.58489319246111,0.0398107170553497,0.398107170553497,0.00158489319246111,0.000630957344480193,2.51188643150958,0.630957344480194,1.0,1.58489319246111,0.00158489319246111,1.58489319246111,0.158489319246111,0.251188643150958 };
	double phro[15] = { 0.0,0.630957344480194,0.0100000000000000,3.98107170553497,1.58489319246111,15.8489319246111,15.8489319246111,0.251188643150958,3.98107170553497,2.51188643150958,0.251188643150958,15.8489319246111,0.00630957344480193,3.98107170553497,6.30957344480193 };
	double beta[15] = { 0.0,0.0158489319246111,25.1188643150958,6.30957344480193,15.8489319246111,0.000398107170553497,0.00100000000000000,0.158489319246111,1.0,0.0630957344480194,630.957344480193,0.0158489319246111,0.630957344480194,1.58489319246111,0.00251188643150958 };
	double scl_cb[15] = { 0.0,0.0,0.987443835660405,0.0,0.0,0.802171644410405,0.0,2.07309221817751,0.762900884574164,1.32734416516483,0.925111062832940,0.797787522347863,0.610967428488595,0.928270215536784,0.905828222074527 };
	double scl_cr[15] = { 0.0,3.45498410289289,1.87825091785739,3.47599894937224,2.56390718416733,1.57674397937559,1.83407623664399,1.91030648245327,1.39603770369435,1.57449905506255,1.56770562115591,1.37504912388072,1.23015550287436,1.29569497078972,1.54105745038586 };	// parameters from measurements
	int level[15] = { 1,1,1,2,2,2,3,3,3,4,4,4,5,5,5 };
	int band[15] = { 1,1,0,1,1,0,1,1,0,1,1,0,1,1,0 };
	for (i = 0; i < 15; i++)
	{
		if (omega[i] != 0)
		{
			tb[i][0] = t02tb(t0[i], jnd, 3.2, band[i], sig2[i], level[i], pitch, dist,0);
			tb[i][0] = masking(sig2[i], tb[i][0], omega[i], phro[i], beta[i], band[i], level[i], pitch, dist,500);
		}
		else
			tb[i][0] = 128.0;	// infinity, the same case in the following 
	}
	for (i = 0; i < 15; i++)
	{
		if (scl_cb[i] != 0.0)
		{
			tb[i][1] = t02tb(t0[15 + i], jnd, 3.2, band[i], sig2[15 + i], level[i], pitch, dist,0);
			tb[i][1] = (scl_cb[i])*(tb[i][1]);
		}
		else
			tb[i][1] = 128.0;
	}
	for (i = 0; i < 15; i++)
	{
		if (scl_cr[i] != 0.0)
		{
			tb[i][2] = t02tb(t0[30 + i], jnd, 3.2, band[i], sig2[30 + i], level[i], pitch, dist,0);
			tb[i][2] = (scl_cr[i])*(tb[i][2]);
		}
		else
			tb[i][2] = 128.0;
	}
}
