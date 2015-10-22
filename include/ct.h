/*****************************************************************************
** FILE IDENTIFICATION
**
**   Name:          ct.h
**   Purpose:       Master header file for CTSim
**   Programmer:    Kevin Rosenberg
**   Date Started:  Aug 1984
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

#ifndef CT_H
#define CT_H


#include <complex>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <cstring>
#include <cstddef>
#include <cstdarg>
#include <cstdlib>

using namespace std;

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <memory>


#define HAVE_FFTW

#ifdef HAVE_FFTW
#include <fftw3.h>
#define HAVE_FFT 1
#endif


struct ReconstructionROI 
{
	double m_dXMin;
	double m_dYMin;
	double m_dZMin;
	double m_dXMax;
	double m_dYMax;
	double m_dZMax;
};

struct Point
{
	double x;
	double y;
	double z;
};

struct Volume
{
	int nx;
	int ny;
	int nz;
};

#include "array2d.h"

#include "ctsupport.h"

#include "scanner.h"
#include "Projections.h"
#include "backprojectors.h"
#include "filter.h"
#include "fourier.h"
#include "procsignal.h"
#include "reconstruct.h"
#include "Correction.h"


#endif

