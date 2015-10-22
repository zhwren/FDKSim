/*****************************************************************************
** FILE IDENTIFICATION
**
**   Name:          scanner.cpp
**   Purpose:       Classes for CT scanner
**   Programmer:    Kevin Rosenberg
**   Date Started:  1984
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

const int Scanner::GEOMETRY_EQUIANGULAR = 1;
const int Scanner::GEOMETRY_EQUILINEAR = 2;

Scanner::Scanner (int geometryId, int nDetU, int nDetV, double DetUSize, double DetVSize,
	int nView, double dStartView, const double rotAngleInc, 
	double FocalLength, double sourceDetectorLength)
{
	m_idGeometry = geometryId;
	m_nDetU = nDetU;				/* Number of detectors in array */
	m_nDetV = nDetV;				/*Number of detectors in rows*/
	m_dDetUSize = DetUSize;			/*Size or Angle of detector element in array*/
	m_dDetVSize = DetVSize;			/*Size of detector element in rows*/

	m_dFocalLength = FocalLength;						// Focal Length, distance from source to center
	m_dSourceDetectorLength = sourceDetectorLength;		// Distance from source to detectors
	m_dCenterDetectorLength = sourceDetectorLength - FocalLength;		// Distance from center to detectors

	m_nRotView = nView;         /* Number of rotated views */

	m_startView = dStartView;				//Start View angle (norm 0)
	m_rotInc = rotAngleInc;					// Increment in rotation angle between views
	m_rotLen = (m_nRotView + 1)*m_rotInc;	// Rotation angle length in radians (norm 2PI)
}

Scanner::~Scanner (void)
{
}

