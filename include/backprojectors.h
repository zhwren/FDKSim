/*****************************************************************************
** FILE IDENTIFICATION
**
**      Name:         backproject.h
**      Purpose:      Backprojection classes
**      Programmer:   Kevin Rosenberg
**      Date Started: June 2000
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


#ifndef __BACKPROJECTORS_H
#define __BACKPROJECTORS_H

#undef HAVE_BSPLINE_INTERP

#include "ct.h"


class Backproject;
class ImageFile;
class Projections;
struct ReconstructionROI;

class Backprojector
{
 public:
  static const int INTERP_INVALID;
  static const int INTERP_NEAREST;
  static const int INTERP_LINEAR;
  static const int INTERP_CUBIC;
  static const int INTERP_FREQ_PREINTERPOLATION;

#if HAVE_BSPLINE_INTERP
  static const int INTERP_BSPLINE;
  static const int INTERP_1BSPLINE;
  static const int INTERP_2BSPLINE;
  static const int INTERP_3BSPLINE;
#endif

  Backprojector (Projections* proj, double* im, int interpID, const int interpFactor, Volume vol, const ReconstructionROI* pROI);

  ~Backprojector ();

  void BackprojectView (const double* const viewData, const double viewAngle);
  void PostProcessing();

  bool fail() const {return m_fail;}
  const string& failMessage() const {return m_failMessage;}

 private:

  bool m_fail;
  string m_failMessage;

  int m_idInterpolation;
  Backproject* m_pBackprojectImplem;

  bool init(Projections* proj, double* im, int interpID, const int interpFactor, Volume vol, const ReconstructionROI* pROI);
};


class Backproject
{
 public:
    Backproject (Projections* proj, double* im, int interpID, const int interpFactor, Volume vol, const ReconstructionROI* pROI);

    virtual ~Backproject ();

    virtual void BackprojectView (const double* const viewData, const double viewAngle) = 0;
    virtual void PostProcessing (); // call after backprojecting all views

 protected:
    //void ScaleImageByRotIncrement ();

	double *im;

    int interpType;
	int interpFactor;

	int nViews;
	double rotInc;

    double rotScale;

    int iDetCenter;             // index refering to Z=0 projection
    int nDetU, nDetV;
	double detUInc, detVInc;
   
	int nx, ny, nz;
    double xInc, yInc, zInc;	// size of cells
	double xMin, yMin, zMin;
	double xMax, yMax, zMax;

    int m_interpFactor;

    double m_dFocalLength;
    double m_dSourceDetectorLength;

    bool m_bPostProcessingDone;

	double *r;
	double *phi;

 private:
    Backproject (const Backproject& rhs);
    Backproject& operator= (const Backproject& rhs);
};


class BackprojectEquilinear: public Backproject
{
 public:
  BackprojectEquilinear (Projections* proj, double* im, int interpID, const int interpFactor, Volume vol, const ReconstructionROI* pROI)
	  : Backproject(proj, im, interpID, interpFactor, vol, pROI)
  {}

  void BackprojectView (const double* const t, const double view_angle);

  ~BackprojectEquilinear()
   {}
};

class BackprojectEquiangular: public Backproject
{
 public:
  BackprojectEquiangular (Projections* proj, double* im, int interpID, const int interpFactor, Volume vol, const ReconstructionROI* pROI)
  : Backproject(proj, im, interpID, interpFactor, vol, pROI)
  {}

  void BackprojectView (const double* const t, const double view_angle);

  ~BackprojectEquiangular()
  {}
};


#endif
