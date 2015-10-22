#include "ct.h"

Correction::Correction(int iGeometry, int nDetU, int nDetV, double dUSize, double dVSize, double dFocalLength, double dSourceDetectorLength)
{
	m_dCorrectionMatix = NULL;
	m_fail = false;

	init(iGeometry, nDetU, nDetV, dUSize, dVSize, dFocalLength, dSourceDetectorLength);
}

Correction::~Correction()
{
	if(m_dCorrectionMatix != NULL)
	{
		delete[] m_dCorrectionMatix;
		m_dCorrectionMatix = NULL;
	}
}

void Correction::init(int iGeometry, int nDetU, int nDetV, double dUSize, double dVSize, double dFocalLength, double dSourceDetectorLength)
{
	if( iGeometry != Scanner::GEOMETRY_EQUIANGULAR && iGeometry != Scanner::GEOMETRY_EQUILINEAR)
	{
		m_fail = true;
		m_failMessage = "Invalid Geometry";
		return;
	}

	m_idGeometry = iGeometry;

	m_nDetU = nDetU;
	m_nDetV = nDetV;
	m_dDetUSize = dUSize;
	m_dDetVSize = dVSize;

	m_dFocalLength = dFocalLength;
	m_dSourceDetectorLength = dSourceDetectorLength;

	m_dCorrectionMatix = new double[m_nDetU*m_nDetV];

	if( m_dCorrectionMatix == NULL)
	{
		m_fail = true;
		m_failMessage = "Cannot allocate Memory for CorrectionMatrix";
		return;
	}

	if(m_idGeometry == Scanner::GEOMETRY_EQUIANGULAR)
	{
		for(int iU= 0; iU<m_nDetU; iU++)
		{
			double cosa = cos((iU - m_nDetU/2.0 + 0.5)*m_dDetUSize); // iU center to column center angle,because m_dDetUSize is angle.

			for (int iV= 0; iV<m_nDetV; iV++)
			{
				double dV = ((iV - m_nDetV/2.0 + 0.5)*m_dDetVSize);//m_dDetVSize is length!!!!!!!!!!
				double cosc = m_dSourceDetectorLength/sqrt(dV * dV + m_dSourceDetectorLength * m_dSourceDetectorLength);

				*(m_dCorrectionMatix + m_nDetU * iV + iU) = cosa * cosc;
			}
		}
	}
	else if(m_idGeometry == Scanner::GEOMETRY_EQUILINEAR)
	{
		for(int iU= 0; iU<m_nDetU; iU++)
		{
			for (int iV= 0; iV<m_nDetV; iV++)
			{
				double dU = ((iU - m_nDetU/2.0 + 0.5)*m_dDetUSize);
				double dV = ((iV - m_nDetV/2.0 + 0.5)*m_dDetVSize);

				double cosv = m_dSourceDetectorLength/sqrt(dU * dU + dV * dV + m_dSourceDetectorLength * m_dSourceDetectorLength);

				*(m_dCorrectionMatix + m_nDetU * iV + iU) = cosv;
			}
		}
	}

}

void Correction::ProcessSignal(double *pSino)
{
	if( m_dCorrectionMatix == NULL)
	{
		m_fail = true;
		m_failMessage = "CorrectionMatrix == NULL";
		return;
	}

	for(int iU= 0; iU<m_nDetU; iU++)
	{
		for (int iV= 0; iV<m_nDetV; iV++)
		{
			*(pSino + m_nDetU * iV + iU) *= *(m_dCorrectionMatix + m_nDetU * iV + iU);
		}
	}
}
