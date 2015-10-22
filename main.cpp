#include "ct.h"
#include <time.h>

int main()
{
	int idGeometry = Scanner::GEOMETRY_EQUIANGULAR;

	unsigned int nDetU = 749;
	unsigned int nDetV = 1;
	double dDetUSize = 0.1*PI/180;
	double dDetVSize = 1;

	double dFocalLength = 600;
	double dSourceDetectorLength = 600;

	unsigned int nView = 360;
	double startView = -PI/2;
	double rotInc = -2*PI/nView;

	Scanner *pscanner = new Scanner(idGeometry, nDetU, nDetV, dDetUSize, dDetVSize,nView,startView, rotInc,dFocalLength, dSourceDetectorLength);

	Volume vol;
	vol.nx = 512;
	vol.ny = 512;
	vol.nz = 1;

	ReconstructionROI roi;
	roi.m_dXMin = -256;
	roi.m_dXMax = 256;
	roi.m_dYMin = -256;
	roi.m_dYMax = 256;
	roi.m_dZMin = -0;
	roi.m_dZMax = 0;

	double *pSino = new double[pscanner->nDetU()*pscanner->nDetV()*pscanner->nView()];
	memset(pSino, 0, sizeof(double)*pscanner->nDetU()*pscanner->nDetV()*pscanner->nView());

	FILE *fp = fopen("data720.txt","rb");
	if(fp == NULL)
	{
		printf("Cannot Open sinogram file.\n");
		return -1;
	}
	fread(pSino, sizeof(double), pscanner->nDetU()*pscanner->nDetV()*pscanner->nView(), fp);
	fclose(fp);


	Projections* rProj = new Projections(pscanner, pSino); 
	if(rProj == NULL || rProj->fail() == true)
	{
		printf("Error ProjectionsData.\n");
	}

	double *pReconData = new double[vol.nx*vol.ny*vol.nz];
	memset(pReconData, 0, sizeof(double)*vol.nx*vol.ny*vol.nz);

	int filterID = SignalFilter::FILTER_ABS_BANDLIMIT;
	double filt_param = 0;
	int filterMethodID = ProcessSignal :: FILTER_METHOD_FFTW;
	int zeropad = 2;
	int filterGenerationID = ProcessSignal::FILTER_GENERATION_INVERSE_FOURIER;

	int interpID = Backprojector::INTERP_LINEAR;
	int interpFactor = 1;
	bool bRebinToParallel = false;

	clock_t start, finish;   
	double duration;   

	start = clock();     


	Reconstructor *pRec = new Reconstructor(rProj, pReconData, filterID, filt_param, filterMethodID, zeropad, filterGenerationID, interpID, interpFactor, vol, &roi, bRebinToParallel);
	if(pRec == NULL || pRec->fail())
	{
		printf("Error new Reconstructor.\n");
		return -1;
	}
	pRec->reconstructAllViews();

	finish = clock();   
	duration = (double)(finish - start) / CLOCKS_PER_SEC;   
	printf( "%f seconds\n", duration );  

	FILE *ffp = fopen("Circular_2880.rcn","wb");
	if(ffp == NULL)
	{
		printf("Cannot Open file TO save reconstruction results.\n");
		return -1;
	}
	fwrite(pReconData, sizeof(double), vol.nx*vol.ny*vol.nz, ffp);
	fclose(ffp);

	delete[] pReconData;
	delete[] pSino;
	delete rProj;
	delete pRec;
	delete pscanner;
}

