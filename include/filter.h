/*****************************************************************************
** FILE IDENTIFICATION
**
**      Name:         filter.h
**      Purpose:      Signal filter header file
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

#ifndef FILTER_H
#define FILTER_H


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include <string>
using std::string;

class SignalFilter 
{
 public:
    static const int FILTER_INVALID;
    static const int FILTER_ABS_BANDLIMIT;      // filter times |x|
    static const int FILTER_ABS_COSINE;
    static const int FILTER_ABS_G_HAMMING;
    static const int FILTER_ABS_HANNING;
    static const int DOMAIN_INVALID;
    static const int DOMAIN_FREQUENCY;
    static const int DOMAIN_SPATIAL;
    SignalFilter(const int, double, double, int, double, double, const int);
    ~SignalFilter (void);

    double* getFilter (void) const { return m_adFilter; }
    const int idFilter(void) const { return m_idFilter;}
    const int idDomain(void) const { return m_idDomain;}
    bool fail() { return m_fail; }
    string failMessage() { return m_failMessage; }

    int getNFilterPoints() { return m_nFilterPoints; }
    double getFilterMin() { return m_dFilterMin; }
    double getFilterMax() { return m_dFilterMax; }
    double getFilterIncrement() { return m_dFilterInc; }

    void copyFilterData(double*, const int, const int);
    double response(double x);
    static double spatialResponse(int, double, double, double);
    static double frequencyResponse(int, double, double, double);
    static double spatialResponseAnalytic(int, double, double, double);
    static double spatialResponseCalc(int, double, double, double, int);
    static void setNumIntegral(int nIntegral) { N_INTEGRAL = nIntegral; }
    static double sinc(double x) { return (fabs(x)>F_EPSILON?(sin(x)/x):1.0); }
    static double sinc(double x, double mult) { return (fabs(x)>F_EPSILON?(sin(x*mult)/x):1.0); }
 private:
    int m_nFilterPoints;
    double m_dBandwidth;
    double m_dFilterParam;
    double m_dFilterInc;
    double m_dFilterMin;
    double m_dFilterMax;
    double* m_adFilter;

    int m_idFilter;
    int m_idDomain;
    static int N_INTEGRAL;
    bool m_fail;
    string m_failMessage;
    static const bool haveAnalyticSpatial(const int filterID);
    void init(const int, double, double, int, double, double, const int);
    void createFrequencyFilter(double*);
    void createSpatialFilter(double*);
    double spatialResponseCalc(double);
    double spatialResponseAnalytic(double);
    double frequencyResponse(double);
    static double integral_abscos (double u, double w)
    { return (fabs(u)>F_EPSILON?(cos(u*w)-1)/(u*u)+w/u*sin(u*w):(w*w/2)); }
};

#endif
