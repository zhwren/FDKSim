/*****************************************************************************
** FILE IDENTIFICATION
**
**   Name:          scanner.h
**   Purpose:       Scanner header file
**   Programmer:    Kevin Rosenberg
**   Date Started:  July 1, 1984
**
**  This is part of the CTSim program
**  Copyright (C) 1983-2009 Kevin Rosenberg
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

#ifndef SCANNER_H
#define SCANNER_H

// Projections are collected along an array of ndet detectors.  The data
// for these detectors is stored in the class DetectorArray

class Scanner
{
 public:
  static const int GEOMETRY_EQUILINEAR;
  static const int GEOMETRY_EQUIANGULAR;
  Scanner (int, int, int, double, double, int, double, const double, double, double);
  ~Scanner();

  void setNView(int);
  void setOffsetView(int);

  bool fail() {return m_fail;}
  string failMessage() {return m_failMessage;}

  int geometry() const {return m_idGeometry;}
  unsigned int nDetU() const {return m_nDetU;}
  unsigned int nDetV() const {return m_nDetV;}
  double detIncU() const {return m_dDetUSize;}
  double detIncV() const {return m_dDetVSize;}

  unsigned int nView() const {return m_nRotView;}

  double startView() const {return m_startView;}
  double rotInc() const {return m_rotInc;}
  double rotLen() const {return m_rotLen;}

  double pitchLen() const {return m_pitchLen;}
  double pitchInc() const {return m_pitchInc;}

  double focalLength() const {return m_dFocalLength;}
  double sourceDetectorLength() const {return m_dSourceDetectorLength;}
  double centerDetectorLength() const {return m_dCenterDetectorLength;}


 private:
  bool m_fail;
  string m_failMessage;
  int m_idGeometry;

  unsigned int m_nDetU;         /* Number of detectors in array */
  unsigned int m_nDetV;			/*Number of detectors in rows*/
  double m_dDetUSize;			/*Size or Angle of detector element in array*/
  double m_dDetVSize;			/*Size of detector element in rows*/

  double m_dFocalLength;        // Focal Length, distance from source to center
  double m_dSourceDetectorLength; // Distance from source to detectors
  double m_dCenterDetectorLength; // Distance from center to detectors

  unsigned int m_nRotView;         /* Number of rotated views */

  double m_startView;			//Start View angle (norm 0)
  double m_rotInc;              // Increment in rotation angle between views
  double m_rotLen;              // Rotation angle length in radians (norm 2PI)

  double m_pitchLen;			//Z axis length
  double m_pitchInc;			// Increment in pith between views
};


#endif
