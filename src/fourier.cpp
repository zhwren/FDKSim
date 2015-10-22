///*****************************************************************************
//** FILE IDENTIFICATION
//**
//**   Name:          fourier.cpp
//**   Purpose:       Fourier transform functions
//**   Programmer:    Kevin Rosenberg
//**   Date Started:  Dec 2000
//**
//**  This is part of the CTSim program
//**  Copyright (c) 1983-2009 Kevin Rosenberg
//**
//**  This program is free software; you can redistribute it and/or modify
//**  it under the terms of the GNU General Public License (version 2) as
//**  published by the Free Software Foundation.
//**
//**  This program is distributed in the hope that it will be useful,
//**  but WITHOUT ANY WARRANTY; without even the implied warranty of
//**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//**  GNU General Public License for more details.
//**
//**  You should have received a copy of the GNU General Public License
//**  along with this program; if not, write to the Free Software
//**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//******************************************************************************/

#include "ct.h"


#ifdef HAVE_FFTW
void Fourier::shuffleNaturalToFourierOrder (fftw_complex* pVector, const int n)
{
  fftw_complex* pTemp = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
  int i;

  if (isOdd(n)) { // Odd
    int iHalfN = (n - 1) / 2;

    pTemp[0][0] = pVector[iHalfN][0];
    pTemp[0][1] = pVector[iHalfN][1];
    for (i = 0; i < iHalfN; i++) {
      pTemp[i + 1][0] = pVector[i + 1 + iHalfN][0];
      pTemp[i + 1][1] = pVector[i + 1 + iHalfN][1];
    }
    for (i = 0; i < iHalfN; i++) {
      pTemp[i + iHalfN + 1][0] = pVector[i][0];
      pTemp[i + iHalfN + 1][1] = pVector[i][1];
    }
  } else {     // Even
    int iHalfN = n / 2;
    pTemp[0][0] = pVector[iHalfN][0];
    pTemp[0][1] = pVector[iHalfN][1];
    for (i = 0; i < iHalfN - 1; i++) {
      pTemp[i + 1][0] = pVector[i + iHalfN + 1][0];
      pTemp[i + 1][1] = pVector[i + iHalfN + 1][1];
    }
    for (i = 0; i < iHalfN; i++) {
      pTemp[i + iHalfN][0] = pVector[i][0];
      pTemp[i + iHalfN][1] = pVector[i][1];
    }
  }

  for (i = 0; i < n; i++) {
    pVector[i][0] = pTemp[i][0];
    pVector[i][1] = pTemp[i][1];
  }
  fftw_free(pTemp);
}

void Fourier::shuffleFourierToNaturalOrder (fftw_complex* pVector, const int n)
{
  fftw_complex* pTemp = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * n));
  int i;
  if (isOdd(n)) { // Odd
    int iHalfN = (n - 1) / 2;

    pTemp[iHalfN][0] = pVector[0][0];
    pTemp[iHalfN][1] = pVector[0][1];
    for (i = 0; i < iHalfN; i++) {
      pTemp[i + 1 + iHalfN][0] = pVector[i + 1][0];
      pTemp[i + 1 + iHalfN][1] = pVector[i + 1][1];
    }
    for (i = 0; i < iHalfN; i++) {
      pTemp[i][0] = pVector[i + iHalfN + 1][0];
      pTemp[i][1] = pVector[i + iHalfN + 1][1];
    }
  } else {     // Even
    int iHalfN = n / 2;
    pTemp[iHalfN][0] = pVector[0][0];
    pTemp[iHalfN][1] = pVector[0][1];
    for (i = 0; i < iHalfN; i++) {
      pTemp[i][0] = pVector[i + iHalfN][0];
      pTemp[i][1] = pVector[i + iHalfN][1];
    }
    for (i = 0; i < iHalfN - 1; i++) {
      pTemp[i + iHalfN + 1][0] = pVector[i+1][0];
      pTemp[i + iHalfN + 1][1] = pVector[i+1][1];
    }
  }

  for (i = 0; i < n; i++) {
    pVector[i][0] = pTemp[i][0];
    pVector[i][1] = pTemp[i][1];
  }

  fftw_free(pTemp);
}
#endif

