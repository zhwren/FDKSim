/*****************************************************************************
** FILE IDENTIFICATION
**
**   Name:         reconstruct.h          Header file for Reconstruction class
**   Programmer:   Kevin Rosenberg
**   Date Started: Aug 84
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

#ifndef __RECONSTRUCT_H
#define __RECONSTRUCT_H


class Backprojector;
class ProcessSignal;
class Correction;

class Reconstructor
{
  private:
    Projections* m_pProjectionDatas;
    ProcessSignal* m_pProcessSignal;
    Backprojector* m_pBackprojector;
    Correction* m_pCorrection;
    int m_nFilteredProjections;
    const bool m_bRebinToParallel;
    bool m_bFail;
    std::string m_strFailMessage;
  public:
    Reconstructor (Projections*, double*, int, double, int, const int, int, 
		 int,int, Volume, ReconstructionROI*, bool bRebinToParallel = false);
    ~Reconstructor();
    bool fail() const {return m_bFail;}
    const std::string& failMessage() const {return m_strFailMessage;}
    void reconstructAllViews ();
    void reconstructView (int iStartView = 0, int iViewCount = -1);
    void postProcessing ();
};

#endif
