/*****************************************************************************
** File IDENTIFICATION
**
**     Name:                   filter.cpp
**     Purpose:                Routines for signal-procesing filters
**     Progammer:              Kevin Rosenberg
**     Date Started:           Aug 1984
**
**  This is part of the CTSim program
**  Copyright (c) 1983-2009 Kevin Rosenberg
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License (version 2) as
**  published by the Free Software Foundation.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/
#include "ct.h"

int SignalFilter::N_INTEGRAL=500;  //static member
const int SignalFilter::FILTER_INVALID = -1 ;
const int SignalFilter::FILTER_ABS_BANDLIMIT = 0;       // filter times |x|
const int SignalFilter::FILTER_ABS_G_HAMMING = 1;
const int SignalFilter::FILTER_ABS_HANNING = 2;
const int SignalFilter::FILTER_ABS_COSINE = 3;
const int SignalFilter::DOMAIN_INVALID = -1;
const int SignalFilter::DOMAIN_FREQUENCY = 0;
const int SignalFilter::DOMAIN_SPATIAL = 1;

SignalFilter::SignalFilter (const int idFilter, double dFilterMinimum, double dFilterMaximum, int nFilterPoints, double dBandwidth, double dFilterParam, const int idDomain)
{
  init (idFilter, dFilterMinimum, dFilterMaximum, nFilterPoints, dBandwidth, dFilterParam, idDomain);
}

void SignalFilter::init (const int idFilter, double dFilterMinimum, double dFilterMaximum, int nFilterPoints, double dBandwidth, double dFilterParam, const int idDomain)
{
  m_idFilter = idFilter;
  m_idDomain = idDomain;
  m_nFilterPoints = nFilterPoints;
  m_dFilterParam = dFilterParam;
  m_dBandwidth = dBandwidth;
  m_dFilterMin = dFilterMinimum;
  m_dFilterMax = dFilterMaximum;
  m_dFilterInc = (m_dFilterMax - m_dFilterMin) / (m_nFilterPoints - 1);
  m_adFilter = new double [m_nFilterPoints];
  if (m_idDomain == DOMAIN_FREQUENCY) createFrequencyFilter (m_adFilter);
  else if (m_idDomain == DOMAIN_SPATIAL) createSpatialFilter (m_adFilter);
}

SignalFilter::~SignalFilter ()
{
  delete [] m_adFilter;
}

void SignalFilter::createFrequencyFilter (double* adFilter)
{
  double x;
  int i;
  for (x = m_dFilterMin, i = 0; i < m_nFilterPoints; x += m_dFilterInc, i++)
    adFilter[i] = frequencyResponse(x);
}


void SignalFilter::createSpatialFilter (double* adFilter)
{
    double x = m_dFilterMin;

    for (int i = 0; i < m_nFilterPoints; i++, x += m_dFilterInc) 
	{
      if (haveAnalyticSpatial(m_idFilter))
        m_adFilter[i] = spatialResponseAnalytic (x);
      else
        m_adFilter[i] = spatialResponseCalc (x);
    }
}




double SignalFilter::response (double x)
{
  double response = 0;

  if (m_idDomain == DOMAIN_SPATIAL)
    response = spatialResponse (m_idFilter, m_dBandwidth, x, m_dFilterParam);
  else if (m_idDomain == DOMAIN_FREQUENCY)
    response = frequencyResponse (m_idFilter, m_dBandwidth, x, m_dFilterParam);

  return (response);
}


double SignalFilter::spatialResponse (int filterID, double bw, double x, double param)
{
  if (haveAnalyticSpatial(filterID))
    return spatialResponseAnalytic (filterID, bw, x, param);
  else
    return spatialResponseCalc (filterID, bw, x, param, N_INTEGRAL);
}

void SignalFilter::copyFilterData (double* pdFilter, const int iStart, const int nPoints)
{
  int iFirst = clamp (iStart, 0, m_nFilterPoints - 1);
  int iLast = clamp (iFirst + nPoints - 1, 0, m_nFilterPoints - 1);

  for (int i = iFirst; i <= iLast; i++)
    pdFilter[i - iFirst] = m_adFilter[i];
}

/* NAME
*   filter_spatial_response_calc        Calculate filter by discrete inverse fourier
*                                       transform of filters's frequency
*                                       response
*
* SYNOPSIS
*   y = filter_spatial_response_calc (filt_type, x, m_bw, param, n)
*   double y                    Filter's response in spatial domain
*   int filt_type               Type of filter (definitions in ct.h)
*   double x                    Spatial position to evaluate filter
*   double m_bw                 Bandwidth of window
*   double param                General parameter for various filters
*   int n                       Number of points to calculate integrations
*/

double SignalFilter::spatialResponseCalc (double x)
{
  return (spatialResponseCalc (m_idFilter, m_dBandwidth, x, m_dFilterParam, N_INTEGRAL));
}


double SignalFilter::spatialResponseCalc (int filterID, double bw, double x, double param, int n)
{
  double zmin, zmax;

  double zinc = (zmax - zmin) / (n - 1);

  double z = zmin;
  double* q = new double [n];

  for (int i = 0; i < n; i++, z += zinc)
    q[i] = frequencyResponse (filterID, bw, z, param) * cos (TWOPI * z * x);

  double y = 2 * integrateSimpson (zmin, zmax, q, n);
  delete q;

  return (y);
}

double SignalFilter::frequencyResponse (double u)
{
  return frequencyResponse (m_idFilter, m_dBandwidth, u, m_dFilterParam);
}


double SignalFilter::frequencyResponse (int filterID, double bw, double u, double param)
{
  double q;
  double au = fabs (u);
  double abw = fabs (bw);

  switch (filterID) 
  {
  case FILTER_ABS_BANDLIMIT:
    if (au >= (abw / 2) + F_EPSILON)
      q = 0.;
    else
      q = au;
    break;
  case FILTER_ABS_COSINE:
    if (au >= (abw / 2) + F_EPSILON)
      q = 0;
    else if (au < F_EPSILON)
      q = 1;
	else
      q = au * cos(PI * au / abw);
    break;
  case FILTER_ABS_HANNING:
    param = 0.5;
    // follow through to ABS_G_HAMMING
  case FILTER_ABS_G_HAMMING:
    if (au >= (abw / 2) + F_EPSILON)
      q = 0;
    else
      q = au * (param + (1 - param) * cos(TWOPI * au / abw));
    break;
  default:
    q = 0;
    //sys_error (ERR_WARNING, "Frequency response for filter %d not implemented [filter_frequency_response]", filterID);
    break;
  }

  return (q);
}



/* NAME
*   filter_spatial_response_analytic                    Calculate filter by analytic inverse fourier
*                               transform of filters's frequency
*                               response
*
* SYNOPSIS
*   y = filter_spatial_response_analytic (filt_type, x, m_bw, param)
*   double y                    Filter's response in spatial domain
*   int filt_type               Type of filter (definitions in ct.h)
*   double x                    Spatial position to evaluate filter
*   double m_bw                 Bandwidth of window
*   double param                General parameter for various filters
*/

double
SignalFilter::spatialResponseAnalytic (double x)
{
  return spatialResponseAnalytic (m_idFilter, m_dBandwidth, x, m_dFilterParam);
}

const bool
SignalFilter::haveAnalyticSpatial (int filterID)
{
  bool haveAnalytic = false;

  switch (filterID) 
  {
	  case FILTER_ABS_BANDLIMIT:
	  case FILTER_ABS_COSINE:
	  case FILTER_ABS_G_HAMMING:
	  case FILTER_ABS_HANNING:
		haveAnalytic = true;
		break;
	  default:
		break;
  }

  return (haveAnalytic);
}

double SignalFilter::spatialResponseAnalytic (int filterID, double bw, double x, double param)
{
  double q, temp;
  double u = TWOPI * x;
  double w = bw / 2;
  double b = PI / bw;
  double b2 = TWOPI / bw;

  switch (filterID) 
  {
	  case FILTER_ABS_BANDLIMIT:
		q = 2 * integral_abscos (u, w);
		break;
	  case FILTER_ABS_COSINE:
		q = integral_abscos(b-u,w) + integral_abscos(b+u,w);
		break;
	  case FILTER_ABS_HANNING:
		param = 0.5;
		// follow through to ABS_G_HAMMING
	  case FILTER_ABS_G_HAMMING:
		q = 2 * param * integral_abscos(u,w) +
		  (1-param)*(integral_abscos(u-b2,w)+integral_abscos(u+b2,w));
		break;
	  default:
		//sys_error (ERR_WARNING, "Analytic filter type %d not implemented [filter_spatial_response_analytic]", filterID);
		q = 0;
		break;
  }

  return (q);
}



// Functions that are inline in filter.h


//  sinc                        Return sin(x)/x function
//   v = sinc (x, mult)
// Calculates sin(x * mult) / x;

//  integral_abscos     Returns integral of u*cos(u)
//
//   q = integral_abscos (u, w)
//   double q                   Integral value
//   double u                   Integration variable
//   double w                   Upper integration boundary
// Returns the value of integral of u*cos(u)*dV for V = 0 to w



