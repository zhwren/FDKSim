#ifndef PROCSIGNAL_H
#define PROCSIGNAL_H


#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include <complex>

class SignalFilter;

typedef std::complex<double> CTSimComplex;


class ProcessSignal 
{
 public:
    static const int FILTER_METHOD_INVALID;
    static const int FILTER_METHOD_CONVOLUTION;
    static const int FILTER_METHOD_FOURIER;
    static const int FILTER_METHOD_FOURIER_TABLE;
    static const int FILTER_METHOD_FFT;

#ifdef HAVE_FFTW
	static const int FILTER_METHOD_FFTW;
	static const int FILTER_METHOD_RFFTW;
#endif

    static const int FILTER_GENERATION_INVALID;
    static const int FILTER_GENERATION_DIRECT;
    static const int FILTER_GENERATION_INVERSE_FOURIER;

    enum 
	{
      FORWARD = -1,
      BACKWARD = 1,
    };

    ProcessSignal (int FilterID, int FilterMethodID, double bw, double signalIncrement,
      int nSignal, double param, int FilterDomainID, int FilterGenerationID,
      const int zeropad, const int preinterpolationFactor, 
	  int iGeometry, double dFocalLength, double dSourceDetectorLength);

    ~ProcessSignal();

	bool fail(void) const       {return m_fail;}
	const std::string& failMessage(void) const {return m_failMessage;}

    void filterSignal (double input[], double output[]) const;


    int getNFilterPoints (void) const  { return m_nFilterPoints; }
    const double getFilterMin(void) const {return m_dFilterMin;}
    const double getFilterMax(void) const {return m_dFilterMax;}
    const double getFilterIncrement(void) const {return m_dFilterInc;}

    double* getFilter(void) {return m_adFilter;}
    const double* getFilter(void) const {return m_adFilter;}

    const int idFilterGeneration() const { return m_idFilterGeneration;}

    // transforms using direct trigometric calculation
    static void finiteFourierTransform (const double input[], double output[], const int n, const int direction);
    static void finiteFourierTransform (const double input[], std::complex<double> output[], const int n, const int direction);
    static void finiteFourierTransform (const std::complex<double> input[], std::complex<double> output[], const int n, const int direction);
    static void finiteFourierTransform (const std::complex<double> input[], double output[], const int n, const int direction);

    static int addZeropadFactor (int n, int iZeropad);

 private:
    int m_idFilterMethod;
    int m_idFilterGeneration;

    int m_nSignalPoints;

    double* m_adFourierCosTable;
    double* m_adFourierSinTable;

    int m_nFilterPoints;
    double m_dSignalInc;
    double m_dFilterInc;
    double m_dFilterMin;
    double m_dFilterMax;
    double* m_adFilter;
    bool m_bFrequencyFiltering;

    // Variables also kept in SignalFilter class
    int m_idFilter;
    int m_idDomain;

    double m_dBandwidth;
    double m_dFilterParam;
    int m_iZeropad;
    int m_nOutputPoints;
    int m_iPreinterpolationFactor;
    int m_idGeometry;
    double m_dFocalLength;
    double m_dSourceDetectorLength;

    static const char* const s_aszFilterMethodName[];
    static const char* const s_aszFilterMethodTitle[];
    static const int s_iFilterMethodCount;
    static const char* const s_aszFilterGenerationName[];
    static const char* const s_aszFilterGenerationTitle[];
    static const int s_iFilterGenerationCount;

	bool m_fail;
	std::string m_failMessage;

#ifdef HAVE_FFTW
    double *m_adRealFftInput, *m_adRealFftOutput, *m_adRealFftSignal, *m_adRealFftBackwardOutput;
    fftw_plan m_realPlanForward, m_realPlanBackward;
    fftw_complex *m_adComplexFftInput, *m_adComplexFftOutput, *m_adComplexFftSignal, *m_adComplexFftBackwardOutput;
    fftw_plan m_complexPlanForward, m_complexPlanBackward;
#endif

    void init (const int idFilter, int idFilterMethod, double dBandwidth, double dSignalIncrement,
      int nSignalPoints, double dFilterParam, const int idDomain, int idFilterGeneration, const int iZeropad,
      const int iPreinterpolationFactor, const int iGeometry, double dFocalLength,
      double dSourceDetectorLength);

    // transforms that use precalculated trig tables, therefore don't
    // require number of data points (n) as an argument
    void finiteFourierTransform (const double input[], std::complex<double> output[], const int direction) const;
    void finiteFourierTransform (const std::complex<double> input[], std::complex<double> output[], const int direction) const;
    void finiteFourierTransform (const std::complex<double> input[], double output[], const int direction) const;

    double convolve (const double func[], const double filter[], const double dx, const int n, const int np) const;
    double convolve (const double f[], const double dx, const int n, const int np) const;
    double convolve (const float f[], const double dx, const int n, const int np) const;

};


#endif
