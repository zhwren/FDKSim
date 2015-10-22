/*****************************************************************************
** FILE IDENTIFICATION
**
**   Name:         backprojectors.cpp         Classes for backprojection
**   Programmer:   Kevin Rosenberg
**   Date Started: June 2000
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
#include "interpolator.h"

const int Backprojector::INTERP_INVALID = -1;
const int Backprojector::INTERP_NEAREST = 0;
const int Backprojector::INTERP_LINEAR = 1;
const int Backprojector::INTERP_CUBIC = 2;
const int Backprojector::INTERP_FREQ_PREINTERPOLATION = 3;
#if HAVE_BSPLINE_INTERP
const int Backprojector::INTERP_BSPLINE = 4;
const int Backprojector::INTERP_1BSPLINE = 5;
const int Backprojector::INTERP_2BSPLINE = 6;
const int Backprojector::INTERP_3BSPLINE = 7;
#endif


Backprojector::Backprojector (Projections* proj, double* im, int interpID, const int interpFactor, Volume vol, const ReconstructionROI* pROI)
{
  m_pBackprojectImplem = NULL;
  m_fail = false;

  init(proj, im, interpID, interpFactor, vol, pROI);
}

void Backprojector::BackprojectView (const double* const viewData, const double viewAngle)
{
  if (m_pBackprojectImplem != NULL)
    m_pBackprojectImplem->BackprojectView (viewData, viewAngle);
}

void Backprojector::PostProcessing()
{
  if (m_pBackprojectImplem != NULL)
    m_pBackprojectImplem->PostProcessing();
}

Backprojector::~Backprojector ()
{
	if (m_pBackprojectImplem != NULL)
	{
		delete m_pBackprojectImplem;
		m_pBackprojectImplem = NULL;
	}
}

// FUNCTION IDENTIFICATION
//     Backproject* projector = selectBackprojector (...)
//
// PURPOSE
//     Selects a backprojector based on BackprojType
//     and initializes the backprojector

bool Backprojector::init(Projections* proj, double* im, int interpID, const int interpFactor, Volume vol, const ReconstructionROI* pROI)
{
  m_idInterpolation = interpID;

  /*if (m_idInterpolation == INTERP_INVALID) 
  {
  m_fail = true;
  m_failMessage = "Invalid interpolation method ";
  }*/

  if (proj->geometry() == Scanner::GEOMETRY_EQUILINEAR)
    m_pBackprojectImplem = static_cast<Backproject*>(new BackprojectEquilinear(proj, im, m_idInterpolation, interpFactor, vol, pROI));
  else if (proj->geometry() == Scanner::GEOMETRY_EQUIANGULAR)
    m_pBackprojectImplem = static_cast<Backproject*>(new BackprojectEquiangular(proj, im, m_idInterpolation, interpFactor, vol, pROI));

  return true;
}


// CLASS IDENTICATION
//   Backproject
//
// PURPOSE
//   Pure virtual base class for all backprojectors.
Backproject::Backproject (Projections* proj, double *im, int interpType, const int interpFactor,
                          Volume vol,const ReconstructionROI* pROI)
{
	nDetU = proj->nDetU();
	nDetV = proj->nDetV();
	detUInc = proj->detIncU();
	detVInc = proj->detIncV();

	m_dFocalLength = proj->focalLength();
	m_dSourceDetectorLength = proj->sourceDetectorLength();

	rotInc = proj->rotInc();

	nx = vol.nx;
	ny = vol.ny;
	nz = vol.nz;

	xInc = (pROI->m_dXMax - pROI->m_dXMin)/nx;
	yInc = (pROI->m_dYMax - pROI->m_dYMin)/ny;
	zInc = (pROI->m_dZMax - pROI->m_dZMin)/nz;

	xMin = pROI->m_dXMin;
	yMin = pROI->m_dYMin;
	zMin = pROI->m_dZMin;

	xMax = pROI->m_dXMax;
	yMax = pROI->m_dYMax;
	zMax = pROI->m_dZMax;

	this->interpType = interpType;
	this->interpFactor = interpFactor;
	this->im = im;

	phi = new double[nx*ny];
	r = new double[nx*ny];

	for(int ix =0; ix <nx; ix++)
	{
		for(int iy=0; iy<ny; iy++)
		{
			double x = xMin + (ix+0.5)*xInc;
			double y = yMin + (iy+0.5)*yInc;

			*(r + nx*iy + ix) = sqrt (x * x + y * y);
			*(phi + nx*iy + ix) = atan2 (y, x);
		}
	}
}

Backproject::~Backproject ()
{
	if(phi!=NULL)
		delete[] phi;

	if(r != NULL)
		delete[] r;
}

void Backproject::PostProcessing()
{
  m_bPostProcessingDone = true;
}



void BackprojectEquiangular::BackprojectView (const double* const filteredProj, const double view_angle)
{
  double beta = view_angle;

  for(int iz = 0; iz < nz; iz ++)
  {
	  for (int ix = 0; ix < nx; ix++) 
	  {
		for (int iy = 0; iy < ny; iy++) 
		{
			double x = xMin + (ix+0.5)*xInc;
			double y = yMin + (iy+0.5)*yInc;
			double z = zMin + (iz+0.5)*zInc;

			double dAngleDiff = beta - *(phi + nx * iy + ix);
			double rcos_t = *(r + nx * iy + ix) * cos (dAngleDiff);
			double rsin_t = *(r + nx * iy + ix) * sin (dAngleDiff);

			/*double rcos_t = x * sin(beta) + y * cos(beta);
			double rsin_t = - x * cos(beta) + y * sin(beta);*/

			double dFLPlusSin = m_dFocalLength + rsin_t;
			double gamma =  atan (rcos_t / dFLPlusSin);

			double dPosU = gamma / detUInc;  // position along detector
			double dPosV = z*(rcos_t / dFLPlusSin) /detVInc ;

			double dL2 = dFLPlusSin * dFLPlusSin + (rcos_t * rcos_t);

			if (interpType == Backprojector::INTERP_NEAREST) 
			{
				int iDetPosU = nearest<int>(dPosU + (nDetU-1)/2.0);  // calc index in the filtered raysum vector
				int iDetPosV = nearest<int>(dPosV + (nDetV-1)/2.0);

				if ((iDetPosU >= 0 && iDetPosU < nDetU) && (iDetPosV >= 0 && iDetPosV < nDetV))
				{
					*(im + iz*nx*ny + iy*nx + ix) += (*(filteredProj + iDetPosU + nDetU*iDetPosV))*(abs(rotInc) *(detUInc)*m_dFocalLength/ dL2);
				}
			} 
			else if (interpType == Backprojector::INTERP_LINEAR) 
			{
				double dPosFloorU = floor(dPosU);
				double dPosFloorV = floor(dPosV);

				int iDetPosU = static_cast<int>(dPosFloorU + (nDetU-1)/2.0);
				int iDetPosV = static_cast<int>(dPosFloorV + (nDetV-1)/2.0);

				double fracU = dPosU - dPosFloorU; // fraction distance from det
				double fracV = dPosV - dPosFloorV; // fraction distance from det
        
				if ((iDetPosU >= 0 && iDetPosU < nDetU-1) 
					&& (iDetPosV >= 0 && iDetPosV < nDetV-1));
				{
					double p1 = (1-fracU)* (*(filteredProj + iDetPosU + nDetU*iDetPosV)) + fracU * (*(filteredProj + iDetPosU+1 + nDetU*iDetPosV));
					double p2 = (1-fracU)* (*(filteredProj + iDetPosU + (nDetU+1)*iDetPosV))  + fracU * (*(filteredProj + iDetPosU+1 + (nDetU+1)*iDetPosV));

					*(im + iz*nx*ny + iy*nx + ix) += ((1-fracV)* p1+ fracV* p2) *(abs(rotInc) *(detUInc)*m_dFocalLength/ dL2);
				}
			} 
			else if (interpType == Backprojector::INTERP_CUBIC) 
			{
				//double d = iDetCenter + dPos;           // position along detector
    //    
				//if (d >= 0 && d < nDet)
				//  im[iz*nx*ny + iy*nx + ix] += pCubicInterp->interpolate (d) / dL2;
		  }
		}   // end for y
	  }     // end for x
  }         //end for z
}

void BackprojectEquilinear::BackprojectView (const double* const filteredProj, const double view_angle)
{
  double beta = view_angle;

  for(int iz = 0; iz < nz; iz++)
  {
	  for (int ix = 0; ix < nx; ix++) 
	  {
		for (int iy = 0; iy < ny; iy++) 
		{
			double x = xMin + (ix+0.5)*xInc;
			double y = yMin + (iy+0.5)*yInc;
			double z = zMin + (iz+0.5)*zInc;

			/*double dAngleDiff = beta - *(phi + nx * iy + ix);
			double rcos_t = *(r + nx * iy + ix) * cos (dAngleDiff);
			double rsin_t = *(r + nx * iy + ix) * sin (dAngleDiff);*/

			double rcos_t = x * sin(beta) + y * cos(beta);
			double rsin_t = - x * cos(beta) + y * sin(beta);

			double dU = (m_dFocalLength + rsin_t) / m_dFocalLength;
			double dDetPosU =  rcos_t / dU;
			double dDetPosV = z / dU;

			// Scale for imaginary detector that passes through origin of phantom, see Kak-Slaney Figure 3.22.
			dDetPosU *= m_dSourceDetectorLength / m_dFocalLength;
			dDetPosV *= m_dSourceDetectorLength / m_dFocalLength;

			double dPosU = dDetPosU / detUInc;  // position along detector array
			double dPosV = dDetPosV / detVInc;

			if (interpType == Backprojector::INTERP_NEAREST) 
			{
				int iDetPosU = nearest<int>(dPosU + (nDetU-1)/2.0);  // calc index in the filtered raysum vector
				int iDetPosV = nearest<int>(dPosV+ (nDetV-1)/2.0);

				if (iDetPosU >= 0 && iDetPosU < nDetU 
					&& iDetPosV >=0 && iDetPosV < nDetV)
				{
					*(im + iz*nx*ny + iy*nx + ix) += *(filteredProj + iDetPosU + nDetU*iDetPosV) / (dU * dU) * abs(rotInc) *(detUInc);
				}
			} 
			else if (interpType == Backprojector::INTERP_LINEAR) 
			{
				double dPosFloorU = floor (dPosU);
				double dPosFloorV = floor (dPosV);

				int iDetPosU = nearest<int>(dPosFloorU + (nDetU-1)/2.0);  // calc index in the filtered raysum vector
				int iDetPosV = nearest<int>(dPosFloorV+ (nDetV-1)/2.0);

				double fracU = dPosU - iDetPosU; 
				double fracV = dPosV - iDetPosV; 

				if (iDetPosU >= 0 && iDetPosU < nDetU - 1 && iDetPosV >= 0 && iDetPosV < nDetV - 1)
				{
					double p1 = (1-fracU)* (*(filteredProj + iDetPosU + nDetU*iDetPosV)) + fracU * (*(filteredProj + iDetPosU+1 + nDetU*iDetPosV));
					double p2 = (1-fracU)* (*(filteredProj + iDetPosU + (nDetU+1)*iDetPosV))  + fracU * (*(filteredProj + iDetPosU+1 + (nDetU+1)*iDetPosV));

					*(im + iz*nx*ny + iy*nx + ix) += ((1-fracV)*p1 + fracV* p2) / (dU * dU) * abs(rotInc) *(detUInc);
				}
			} 
			else if (interpType == Backprojector::INTERP_CUBIC) 
			{
				//double d = iDetCenter + dPos;           // position along detector

				//if (d >= 0 && d < nDet)
				//  im[iz*nx*ny + iy*nx + ix] += pCubicInterp->interpolate (d) / (dU * dU);
			}
		}   // end for y
	}     // end for x
  }       //end for z
}

