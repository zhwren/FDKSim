#include "ct.h"

Projections::Projections(Scanner* pScanner, double* pData)
{
	m_fail = false;

	m_pscanner = pScanner;
	m_pData = pData;
}

Projections::Projections(Projections& rProj)
{
	m_pscanner = rProj.scanner();
	m_pData = rProj.prjData();
}

Projections::~Projections (void)
{
}