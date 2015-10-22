#include "ct.h"
#include <fftw3.h>

// FilterMethod ID/Names
const int ProcessSignal::FILTER_METHOD_INVALID = -1;
const int ProcessSignal::FILTER_METHOD_CONVOLUTION = 0;
const int ProcessSignal::FILTER_METHOD_FOURIER = 1;
const int ProcessSignal::FILTER_METHOD_FOURIER_TABLE = 2;
const int ProcessSignal::FILTER_METHOD_FFT = 3;
#ifdef HAVE_FFTW
const int ProcessSignal::FILTER_METHOD_FFTW = 4;
const int ProcessSignal::FILTER_METHOD_RFFTW =5 ;
#endif


// FilterGeneration ID/Names
const int ProcessSignal::FILTER_GENERATION_INVALID = -1;
const int ProcessSignal::FILTER_GENERATION_DIRECT = 0;
const int ProcessSignal::FILTER_GENERATION_INVERSE_FOURIER = 1;



// CLASS IDENTIFICATION
//   ProcessSignal
//
ProcessSignal:: ProcessSignal (int FilterID, int FilterMethodID, double dBandwidth, double dSignalIncrement,
								int nSignalPoints, double dFilterParam, const int FilterDomainID, int FilterGenerationID,
								const int iZeropad, const int iPreinterpolationFactor, 
								int iGeometry, double dFocalLength, double dSourceDetectorLength)
    : m_adFourierCosTable(NULL), m_adFourierSinTable(NULL), m_adFilter(NULL)
{
	m_fail = false;

#ifdef HAVE_FFTW
	if( FilterMethodID < FILTER_METHOD_CONVOLUTION || FilterMethodID > FILTER_METHOD_RFFTW)
	{
		m_fail = true;
		m_failMessage = "Error FilterMethod ";
		return;
	}
#else
	if( FilterMethodID < FILTER_METHOD_CONVOLUTION || FilterMethodID > FILTER_METHOD_FFT)
	{
		m_fail = true;
		m_failMessage = "Error FilterMethod ";
		return;
	}
#endif
	if((FilterGenerationID != FILTER_GENERATION_DIRECT) && (FilterGenerationID != FILTER_GENERATION_INVERSE_FOURIER)) 
	{
		m_fail = true;
		m_failMessage = "Error FilterGeneration ";
		return;
	}

	if(FilterID < SignalFilter::FILTER_ABS_BANDLIMIT || FilterID > SignalFilter::FILTER_ABS_COSINE)
	{
		m_fail = true;
		m_failMessage = "Error FilterID ";
		return;
	}

    if((FilterDomainID != SignalFilter::DOMAIN_FREQUENCY) && (FilterDomainID != SignalFilter::DOMAIN_SPATIAL))
	{
		m_fail = true;
		m_failMessage = "Error FilterDomain ";
		return;
	}


    init (FilterID, FilterMethodID, dBandwidth, dSignalIncrement, 
			nSignalPoints, dFilterParam, FilterDomainID, FilterGenerationID, 
			iZeropad, iPreinterpolationFactor, 
			iGeometry, dFocalLength, dSourceDetectorLength);
}


void ProcessSignal::init (const int idFilter, const int idFilterMethod, double dBandwidth, double dSignalIncrement,
                     int nSignalPoints, double dFilterParam, const int idDomain, const int idFilterGeneration,
                     const int iZeropad, const int iPreinterpolationFactor, 
					 int iGeometry, double dFocalLength, double dSourceDetectorLength)
{
    int i;

    m_idFilter = idFilter;
    m_idDomain = idDomain;
    m_idFilterMethod = idFilterMethod;
    m_idFilterGeneration = idFilterGeneration;

    m_idGeometry = iGeometry;
    m_dFocalLength = dFocalLength;
    m_dSourceDetectorLength = dSourceDetectorLength;

    m_dBandwidth = dBandwidth;	 
    m_nSignalPoints = nSignalPoints;
    m_dSignalInc = dSignalIncrement;

    m_dFilterParam = dFilterParam;
    m_iZeropad = iZeropad;
    m_iPreinterpolationFactor = iPreinterpolationFactor;

    // scale signalInc/BW to adjust for imaginary detector through origin of phantom
    // see Kak-Slaney Fig 3.22, for Collinear diagram
    if (m_idGeometry == Scanner::GEOMETRY_EQUILINEAR)
    {
        double dEquilinearScale = m_dSourceDetectorLength / m_dFocalLength;
        m_dSignalInc /= dEquilinearScale;
        m_dBandwidth *= dEquilinearScale;
    }

    if (m_idFilterMethod == FILTER_METHOD_FFT)
    {
#ifdef HAVE_FFTW
        m_idFilterMethod = FILTER_METHOD_RFFTW;
#else
        m_fail = true;
        m_failMessage = "FFT not yet implemented";
        return;
#endif
    }

    bool m_bFrequencyFiltering = true;
    if (m_idFilterMethod == FILTER_METHOD_CONVOLUTION)
        m_bFrequencyFiltering = false;

    // Spatial-based filtering
    if (! m_bFrequencyFiltering)
    {
        if (m_idFilterGeneration == FILTER_GENERATION_DIRECT)
        {
            m_nFilterPoints = 2 * (m_nSignalPoints - 1) + 1;
            m_dFilterMin = -m_dSignalInc * (m_nSignalPoints - 1);
            m_dFilterMax = m_dSignalInc * (m_nSignalPoints - 1);
            m_dFilterInc = (m_dFilterMax - m_dFilterMin) / (m_nFilterPoints - 1);

            SignalFilter filter (m_idFilter, m_dFilterMin, m_dFilterMax, m_nFilterPoints, m_dBandwidth, m_dFilterParam, SignalFilter::DOMAIN_SPATIAL);
            m_adFilter = new double[ m_nFilterPoints ];
            filter.copyFilterData (m_adFilter, 0, m_nFilterPoints);
        }
        else if (m_idFilterGeneration == FILTER_GENERATION_INVERSE_FOURIER)
        {
            m_nFilterPoints = 2 * (m_nSignalPoints - 1) + 1;
            m_dFilterMin = -1. / (2 * m_dSignalInc);
            m_dFilterMax = 1. / (2 * m_dSignalInc);
            m_dFilterInc = (m_dFilterMax - m_dFilterMin) / (m_nFilterPoints - 1);

            SignalFilter filter (m_idFilter, m_dFilterMin, m_dFilterMax, m_nFilterPoints, m_dBandwidth, m_dFilterParam, SignalFilter::DOMAIN_FREQUENCY);
            m_adFilter = new double[ m_nFilterPoints ];
            double* adFrequencyFilter = new double [m_nFilterPoints];
            filter.copyFilterData (adFrequencyFilter, 0, m_nFilterPoints);

            Fourier::shuffleNaturalToFourierOrder (adFrequencyFilter, m_nFilterPoints);

            ProcessSignal::finiteFourierTransform (adFrequencyFilter, m_adFilter, m_nFilterPoints, FORWARD);
            delete adFrequencyFilter;

            Fourier::shuffleFourierToNaturalOrder (m_adFilter, m_nFilterPoints);

            for (i = 0; i < m_nFilterPoints; i++)
            {
                m_adFilter[i] /= m_dSignalInc;
            }
        }

        if (m_idGeometry == Scanner::GEOMETRY_EQUILINEAR)
        {
            for (i = 0; i < m_nFilterPoints; i++)
                m_adFilter[i] *= 0.5;
        }
        else if (m_idGeometry == Scanner::GEOMETRY_EQUIANGULAR)
        {
            for (i = 0; i < m_nFilterPoints; i++)
            {
                int iDetFromZero = i - ((m_nFilterPoints - 1) / 2);

                double sinScale = 1 / SignalFilter::sinc (iDetFromZero * m_dSignalInc);
                double dScale = 0.5 * sinScale * sinScale;

                m_adFilter[i] *= dScale;
            }
        } // if (geometry)
    } // if (spatial filtering)

    else if (m_bFrequencyFiltering)    // Frequency-based filtering
    {

        if (m_idFilterGeneration == FILTER_GENERATION_DIRECT)
        {
            // calculate number of filter points with zeropadding
            m_nFilterPoints = addZeropadFactor (m_nSignalPoints, m_iZeropad);
            m_nOutputPoints = m_nFilterPoints * m_iPreinterpolationFactor;

            if (isOdd (m_nFilterPoints))   // Odd
            {
                m_dFilterMin = -1. / (2 * m_dSignalInc);
                m_dFilterMax = 1. / (2 * m_dSignalInc);
                m_dFilterInc = (m_dFilterMax - m_dFilterMin) / (m_nFilterPoints - 1);
            }
            else     // Even
            {
                m_dFilterMin = -1. / (2 * m_dSignalInc);
                m_dFilterMax = 1. / (2 * m_dSignalInc);
                m_dFilterInc = (m_dFilterMax - m_dFilterMin) / m_nFilterPoints;
                m_dFilterMax -= m_dFilterInc;
            }

            SignalFilter filter (m_idFilter, m_dFilterMin, m_dFilterMax, m_nFilterPoints, m_dBandwidth,
                                 m_dFilterParam, SignalFilter::DOMAIN_FREQUENCY);
            m_adFilter = new double [m_nFilterPoints];
            filter.copyFilterData (m_adFilter, 0, m_nFilterPoints);

            // This works fairly well. I'm not sure why since scaling for geometries is done on
            // frequency filter rather than spatial filter as it should be.
            // It gives values slightly off than freq/inverse filtering
            if (m_idGeometry == Scanner::GEOMETRY_EQUILINEAR)
            {
                for (i = 0; i < m_nFilterPoints; i++)
                    m_adFilter[i] *= 0.5;
            }
            else if (m_idGeometry == Scanner::GEOMETRY_EQUIANGULAR)
            {
                for (i = 0; i < m_nFilterPoints; i++)
                {
                    int iDetFromZero = i - ((m_nFilterPoints - 1) / 2);

                    double sinScale = 1 / SignalFilter::sinc (iDetFromZero * m_dSignalInc);
                    double dScale = 0.5 * sinScale * sinScale;

                    m_adFilter[i] *= dScale;
                }
            }

            Fourier::shuffleNaturalToFourierOrder (m_adFilter, m_nFilterPoints);

            // FILTERING:  FREQUENCY - INVERSE FOURIER
        }
        else if (m_idFilterGeneration == FILTER_GENERATION_INVERSE_FOURIER)
        {
            // calculate number of filter points with zeropadding
            int nSpatialPoints = 2 * (m_nSignalPoints - 1) + 1;
            m_dFilterMin = -m_dSignalInc * (m_nSignalPoints - 1);
            m_dFilterMax = m_dSignalInc * (m_nSignalPoints - 1);
            m_dFilterInc = (m_dFilterMax - m_dFilterMin) / (nSpatialPoints - 1);
            m_nFilterPoints = nSpatialPoints;

            if (m_iZeropad > 0)
            {
                double logBase2 = log((double)nSpatialPoints) / log((double)2);
                int nextPowerOf2 = static_cast<int>(floor(logBase2));

                if (logBase2 != floor(logBase2))
                    nextPowerOf2++;

                nextPowerOf2 += (m_iZeropad - 1);
                m_nFilterPoints = 1 << nextPowerOf2;
            }

            m_nOutputPoints = m_nFilterPoints * m_iPreinterpolationFactor;

            double* adSpatialFilter = new double [m_nFilterPoints];
            SignalFilter filter (m_idFilter, m_dFilterMin, m_dFilterMax, nSpatialPoints, m_dBandwidth,
                                 m_dFilterParam, SignalFilter::DOMAIN_SPATIAL);
            filter.copyFilterData (adSpatialFilter, 0, nSpatialPoints);

            if (m_idGeometry == Scanner::GEOMETRY_EQUILINEAR)
            {
                for (i = 0; i < nSpatialPoints; i++)
                    adSpatialFilter[i] *= 0.5;
            }
            else if (m_idGeometry == Scanner::GEOMETRY_EQUIANGULAR)
            {
                for (i = 0; i < nSpatialPoints; i++)
                {
                    int iDetFromZero = i - ((nSpatialPoints - 1) / 2);
                    double sinScale = sin (iDetFromZero * m_dSignalInc);

                    if (fabs(sinScale) < 1E-7)
                        sinScale = 1;
                    else
                        sinScale = (iDetFromZero * m_dSignalInc) / sinScale;

                    double dScale = 0.5 * sinScale * sinScale;

                    adSpatialFilter[i] *= dScale;
                }
            }

            for (i = nSpatialPoints; i < m_nFilterPoints; i++)
                adSpatialFilter[i] = 0;

#if 0
		FILE *fp = fopen("filter.txt", "wb");
		fwrite(adSpatialFilter, sizeof(double),m_nFilterPoints,fp);
		fclose(fp);
#endif

            m_adFilter = new double [m_nFilterPoints];

            std::complex<double>* acInverseFilter = new std::complex<double> [m_nFilterPoints];
            finiteFourierTransform (adSpatialFilter, acInverseFilter, m_nFilterPoints, BACKWARD);
            
			delete adSpatialFilter;
            
			for (i = 0; i < m_nFilterPoints; i++)
                m_adFilter[i] = std::abs (acInverseFilter[i]) * m_dSignalInc;

            delete acInverseFilter;
        }
    }

    // precalculate sin and cosine tables for fourier transform
    if (m_idFilterMethod == FILTER_METHOD_FOURIER_TABLE)
    {
        int nFourier = imax (m_nFilterPoints,m_nOutputPoints) * imax (m_nFilterPoints, m_nOutputPoints) + 1;
        double angleIncrement = (2. * PI) / m_nFilterPoints;

        m_adFourierCosTable = new double[ nFourier ];
        m_adFourierSinTable = new double[ nFourier ];

        double angle = 0;

        for (i = 0; i < nFourier; i++)
        {
            m_adFourierCosTable[i] = cos (angle);
            m_adFourierSinTable[i] = sin (angle);
            angle += angleIncrement;
        }
    }

#ifdef HAVE_FFTW
    if (m_idFilterMethod == FILTER_METHOD_FFTW || m_idFilterMethod == FILTER_METHOD_RFFTW)
    {
        for (i = 0; i < m_nFilterPoints; i++)  //fftw uses unnormalized fft
            m_adFilter[i] /= m_nFilterPoints;
    }

    if (m_idFilterMethod == FILTER_METHOD_RFFTW)
    {
        m_adRealFftInput = static_cast<double*>(fftw_malloc (sizeof(double) * m_nFilterPoints));
        m_adRealFftOutput = static_cast<double*>(fftw_malloc (sizeof(double) * m_nFilterPoints));
        m_realPlanForward = fftw_plan_r2r_1d (m_nFilterPoints, m_adRealFftInput, m_adRealFftOutput, FFTW_R2HC, FFTW_ESTIMATE);
        m_adRealFftSignal = static_cast<double*>(fftw_malloc (sizeof(double) *  m_nOutputPoints));
        m_adRealFftBackwardOutput = static_cast<double*>(fftw_malloc (sizeof(double) * m_nOutputPoints));
        m_realPlanBackward = fftw_plan_r2r_1d (m_nOutputPoints, m_adRealFftSignal, m_adRealFftBackwardOutput, FFTW_HC2R, FFTW_ESTIMATE);
        for (i = 0; i < m_nFilterPoints; i++)
            m_adRealFftInput[i] = 0;
    }
    else if (m_idFilterMethod == FILTER_METHOD_FFTW)
    {
        m_adComplexFftInput = static_cast<fftw_complex*>(fftw_malloc (sizeof(fftw_complex) * m_nFilterPoints));
        m_adComplexFftOutput = static_cast<fftw_complex*>(fftw_malloc (sizeof(fftw_complex) * m_nFilterPoints));
        m_complexPlanForward = fftw_plan_dft_1d (m_nFilterPoints, m_adComplexFftInput, m_adComplexFftOutput, FFTW_FORWARD,  FFTW_ESTIMATE);
        m_adComplexFftSignal = static_cast<fftw_complex*>(fftw_malloc (sizeof(fftw_complex) * m_nOutputPoints));
        m_adComplexFftBackwardOutput = static_cast<fftw_complex*>(fftw_malloc (sizeof(fftw_complex) * m_nOutputPoints));
        m_complexPlanBackward = fftw_plan_dft_1d (m_nOutputPoints, m_adComplexFftSignal, m_adComplexFftBackwardOutput, FFTW_BACKWARD,  FFTW_ESTIMATE);

        for (i = 0; i < m_nFilterPoints; i++)
            m_adComplexFftInput[i][0] = m_adComplexFftInput[i][1] = 0;
        for (i = 0; i < m_nOutputPoints; i++)
            m_adComplexFftSignal[i][0] = m_adComplexFftSignal[i][1] = 0;
    }
#endif

}

ProcessSignal::~ProcessSignal (void)
{
    delete [] m_adFourierSinTable;
    delete [] m_adFourierCosTable;
    delete [] m_adFilter;

#ifdef HAVE_FFTW
    if (m_idFilterMethod == FILTER_METHOD_FFTW)
    {
        fftw_destroy_plan(m_complexPlanForward);
        fftw_destroy_plan(m_complexPlanBackward);
        fftw_free (m_adComplexFftInput);
        fftw_free (m_adComplexFftOutput);
        fftw_free (m_adComplexFftSignal);
        fftw_free (m_adComplexFftBackwardOutput);
    }
    if (m_idFilterMethod == FILTER_METHOD_RFFTW)
    {
        fftw_destroy_plan(m_realPlanForward);
        fftw_destroy_plan(m_realPlanBackward);
        fftw_free (m_adRealFftInput);
        fftw_free (m_adRealFftOutput);
        fftw_free (m_adRealFftSignal);
        fftw_free (m_adRealFftBackwardOutput);
    }
#endif
}


void ProcessSignal::filterSignal (double constInput[], double output[]) const
{
    double* input = new double [m_nSignalPoints];
    int i;

    for (i = 0; i < m_nSignalPoints; i++)
        input[i] = constInput[i];

    if (m_idFilterMethod == FILTER_METHOD_CONVOLUTION)
    {
        for (i = 0; i < m_nSignalPoints; i++)
            output[i] = convolve (input, m_dSignalInc, i, m_nSignalPoints);
    }
    else if (m_idFilterMethod == FILTER_METHOD_FOURIER)
    {
        double* inputSignal = new double [m_nFilterPoints];

        for (i = 0; i < m_nSignalPoints; i++)
            inputSignal[i] = input[i];
        for (i = m_nSignalPoints; i < m_nFilterPoints; i++)
            inputSignal[i] = 0;  // zeropad

        std::complex<double>* fftSignal = new std::complex<double> [m_nFilterPoints];

        finiteFourierTransform (inputSignal, fftSignal, m_nFilterPoints, FORWARD);

        delete inputSignal;

        for (i = 0; i < m_nFilterPoints; i++)
            fftSignal[i] *= m_adFilter[i];

        double* inverseFourier = new double [m_nFilterPoints];

        finiteFourierTransform (fftSignal, inverseFourier, m_nFilterPoints, BACKWARD);

        delete fftSignal;

        for (i = 0; i < m_nSignalPoints; i++)
            output[i] = inverseFourier[i];

        delete inverseFourier;
    }
    else if (m_idFilterMethod == FILTER_METHOD_FOURIER_TABLE)
    {
        double* inputSignal = new double [m_nFilterPoints];

        for (i = 0; i < m_nSignalPoints; i++)
            inputSignal[i] = input[i];

        for (i = m_nSignalPoints; i < m_nFilterPoints; i++)
            inputSignal[i] = 0;  // zeropad

        std::complex<double>* fftSignal = new std::complex<double> [m_nFilterPoints];

        finiteFourierTransform (inputSignal, fftSignal, FORWARD);

        delete inputSignal;

        for (i = 0; i < m_nFilterPoints; i++)
            fftSignal[i] *= m_adFilter[i];

        double* inverseFourier = new double [m_nFilterPoints];

        finiteFourierTransform (fftSignal, inverseFourier, BACKWARD);

        delete fftSignal;

        for (i = 0; i < m_nSignalPoints; i++)
            output[i] = inverseFourier[i];
        delete inverseFourier;
    }
#ifdef HAVE_FFTW
    else if (m_idFilterMethod == FILTER_METHOD_RFFTW)
    {
        for (i = 0; i < m_nSignalPoints; i++)
            m_adRealFftInput[i] = input[i];

        fftw_execute (m_realPlanForward);
        for (i = 0; i < m_nFilterPoints; i++)
            m_adRealFftSignal[i] = m_adFilter[i] * m_adRealFftOutput[i];
        for (i = m_nFilterPoints; i < m_nOutputPoints; i++)
            m_adRealFftSignal[i] = 0;

        fftw_execute (m_realPlanBackward);
        for (i = 0; i < m_nSignalPoints * m_iPreinterpolationFactor; i++)
            output[i] = m_adRealFftBackwardOutput[i];
    }
    else if (m_idFilterMethod == FILTER_METHOD_FFTW)
    {
        for (i = 0; i < m_nSignalPoints; i++)
            m_adComplexFftInput[i][0] = input[i];

        fftw_execute (m_complexPlanForward);
        for (i = 0; i < m_nFilterPoints; i++)
        {
            m_adComplexFftSignal[i][0] = m_adFilter[i] * m_adComplexFftOutput[i][0];
            m_adComplexFftSignal[i][1] = m_adFilter[i] * m_adComplexFftOutput[i][1];
        }
        fftw_execute (m_complexPlanBackward);
        for (i = 0; i < m_nSignalPoints * m_iPreinterpolationFactor; i++)
            output[i] = m_adComplexFftBackwardOutput[i][0];
    }
#endif
    delete input;
}


/* NAME
*    convolve                   Discrete convolution of two functions
*
* SYNOPSIS
*    r = convolve (f1, f2, dx, n, np, func_type)
*    double r                   Convolved result
*    double f1[], f2[]          Functions to be convolved
*    double dx                  Difference between successive x values
*    int n                      Array index to center convolution about
*    int np                     Number of points in f1 array
*    int func_type              EVEN or ODD or EVEN_AND_ODD function f2
*
* NOTES
*    f1 is the projection data, its indices range from 0 to np - 1.
*    The index for f2, the filter, ranges from -(np-1) to (np-1).
*    There are 3 ways to handle the negative vertices of f2:
*       1. If we know f2 is an EVEN function, then f2[-n] = f2[n].
*          All filters used in reconstruction are even.
*      2. If we know f2 is an ODD function, then f2[-n] = -f2[n]
*      3. If f2 is both ODD AND EVEN, then we must store the value of f2
*          for negative indices.  Since f2 must range from -(np-1) to (np-1),
*          if we add (np - 1) to f2's array index, then f2's index will
*          range from 0 to 2 * (np - 1), and the origin, x = 0, will be
*          stored at f2[np-1].
*/

double ProcessSignal::convolve (const double func[], const double dx, const int n, const int np) const
{
    double sum = 0.0;

#if UNOPTIMIZED_CONVOLUTION
    for (int i = 0; i < np; i++)
        sum += func[i] * m_adFilter[n - i + (np - 1)];
#else
    double* f2 = m_adFilter + n + (np - 1);
    for (int i = 0; i < np; i++)
        sum += *func++ * *f2--;
#endif

    return (sum * dx);
}


double ProcessSignal::convolve (const float func[], const double dx, const int n, const int np) const
{
    double sum = 0.0;

#if UNOPTIMIZED_CONVOLUTION
    for (int i = 0; i < np; i++)
        sum += func[i] * m_adFilter[n - i + (np - 1)];
#else
    double* f2 = m_adFilter + n + (np - 1);
    for (int i = 0; i < np; i++)
        sum += *func++ * *f2--;
#endif

    return (sum * dx);
}


void ProcessSignal::finiteFourierTransform (const double input[], double output[], const int n, int direction)
{
    std::complex<double>* complexOutput = new std::complex<double> [n];

    finiteFourierTransform (input, complexOutput, n, direction);

    for (int i = 0; i < n; i++)
        output[i] = complexOutput[i].real();
    delete [] complexOutput;
}


void ProcessSignal::finiteFourierTransform (const double input[], std::complex<double> output[], const int n, int direction)
{
    if (direction < 0)
        direction = -1;
    else
        direction = 1;

    double angleIncrement = direction * 2 * PI / n;

    for (int i = 0; i < n; i++)
    {
        double sumReal = 0;
        double sumImag = 0;

        for (int j = 0; j < n; j++)
        {
            double angle = i * j * angleIncrement;
            sumReal += input[j] * cos(angle);
            sumImag += input[j] * sin(angle);
        }

        if (direction < 0)
        {
            sumReal /= n;
            sumImag /= n;
        }

        output[i] = std::complex<double> (sumReal, sumImag);
    }
}


void ProcessSignal::finiteFourierTransform (const std::complex<double> input[], std::complex<double> output[], const int n, int direction)
{
    if (direction < 0)
        direction = -1;
    else
        direction = 1;

    double angleIncrement = direction * 2 * PI / n;

    for (int i = 0; i < n; i++)
    {
        std::complex<double> sum (0,0);

        for (int j = 0; j < n; j++)
        {
            double angle = i * j * angleIncrement;
            std::complex<double> exponentTerm (cos(angle), sin(angle));
            sum += input[j] * exponentTerm;
        }

        if (direction < 0)
        {
            sum /= n;
        }
        output[i] = sum;
    }
}

void ProcessSignal::finiteFourierTransform (const std::complex<double> input[], double output[], const int n, int direction)
{
    if (direction < 0)
        direction = -1;
    else
        direction = 1;

    double angleIncrement = direction * 2 * PI / n;

    for (int i = 0; i < n; i++)
    {
        double sumReal = 0;

        for (int j = 0; j < n; j++)
        {
            double angle = i * j * angleIncrement;
            sumReal += input[j].real() * cos(angle) - input[j].imag() * sin(angle);
        }

        if (direction < 0)
        {
            sumReal /= n;
        }

        output[i] = sumReal;
    }
}

// Table-based routines
void ProcessSignal::finiteFourierTransform (const double input[], std::complex<double> output[], int direction) const
{
    if (direction < 0)
        direction = -1;
    else
        direction = 1;

    for (int i = 0; i < m_nFilterPoints; i++)
    {
        double sumReal = 0, sumImag = 0;

        for (int j = 0; j < m_nFilterPoints; j++)
        {
            int tableIndex = i * j;

            if (direction > 0)
            {
                sumReal += input[j] * m_adFourierCosTable[tableIndex];
                sumImag += input[j] * m_adFourierSinTable[tableIndex];
            }
            else
            {
                sumReal += input[j] * m_adFourierCosTable[tableIndex];
                sumImag -= input[j] * m_adFourierSinTable[tableIndex];
            }
        }

        if (direction < 0)
        {
            sumReal /= m_nFilterPoints;
            sumImag /= m_nFilterPoints;
        }

        output[i] = std::complex<double> (sumReal, sumImag);
    }
}

// (a+bi) * (c + di) = (ac - bd) + (ad + bc)i
void ProcessSignal::finiteFourierTransform (const std::complex<double> input[], std::complex<double> output[], int direction) const
{
    if (direction < 0)
        direction = -1;
    else
        direction = 1;

    for (int i = 0; i < m_nFilterPoints; i++)
    {
        double sumReal = 0, sumImag = 0;

        for (int j = 0; j < m_nFilterPoints; j++)
        {
            int tableIndex = i * j;

            if (direction > 0)
            {
                sumReal += input[j].real() * m_adFourierCosTable[tableIndex]
                           - input[j].imag() * m_adFourierSinTable[tableIndex];
                sumImag += input[j].real() * m_adFourierSinTable[tableIndex]
                           + input[j].imag() * m_adFourierCosTable[tableIndex];
            }
            else
            {
                sumReal += input[j].real() * m_adFourierCosTable[tableIndex]
                           - input[j].imag() * -m_adFourierSinTable[tableIndex];
                sumImag += input[j].real() * -m_adFourierSinTable[tableIndex]
                           + input[j].imag() * m_adFourierCosTable[tableIndex];
            }
        }

        if (direction < 0)
        {
            sumReal /= m_nFilterPoints;
            sumImag /= m_nFilterPoints;
        }
        output[i] = std::complex<double> (sumReal, sumImag);
    }
}

void ProcessSignal::finiteFourierTransform (const std::complex<double> input[], double output[], int direction) const
{
    if (direction < 0)
        direction = -1;
    else
        direction = 1;

    for (int i = 0; i < m_nFilterPoints; i++)
    {
        double sumReal = 0;

        for (int j = 0; j < m_nFilterPoints; j++)
        {
            int tableIndex = i * j;

            if (direction > 0)
            {
                sumReal += input[j].real() * m_adFourierCosTable[tableIndex]
                           - input[j].imag() * m_adFourierSinTable[tableIndex];
            }
            else
            {
                sumReal += input[j].real() * m_adFourierCosTable[tableIndex]
                           - input[j].imag() * -m_adFourierSinTable[tableIndex];
            }
        }

        if (direction < 0)
        {
            sumReal /= m_nFilterPoints;
        }

        output[i] = sumReal;
    }
}


int ProcessSignal::addZeropadFactor (int n, int iZeropad)
{
    if (iZeropad > 0)
    {
        double dLogBase2 = log((double)n) / log((double)2);
        int iLogBase2 = static_cast<int>(floor (dLogBase2));
        int iPaddedN = 1 << (iLogBase2 + iZeropad);

        return iPaddedN;
    }

    return n;
}
