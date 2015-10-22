/*****************************************************************************
** FILE IDENTIFICATION
**
**      File Name:      ctsupport.h
**      Author:         Kevin Rosenberg
**      Purpose:        Header file for CT support library
**      Date Started:   Dec. 83
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

#ifndef CTSUPPORT_H
#define CTSUPPORT_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#define STR_MAX_LEN 255
#define STR_SIZE    STR_MAX_LEN+1

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <string>
#include <vector>
#include <algorithm>


#ifndef TRUE
#define TRUE    1
#endif
#ifndef FALSE
#define FALSE   0
#endif
#define OK      TRUE

/*----------------------------------------------------------------------*/

#define SHOW(var, fmt)  { cerr << "var = " << var << endl; }

/*----------------------------------------------------------------------*/

#define NEWLINE '\n'
#define TAB     '\t'
#define EOS     '\0'
#define BLANK   ' '



/*----------------------------------------------------------------------*/

#if !defined(O_BINARY) && !defined(MSVC)
#define O_BINARY (0)
#endif

#ifndef S_IWRITE
#define S_IWRITE S_IWUSR
#endif

/*----------------------------------------------------------------------*/

#if defined(MSVC) || ! defined(SIZEOF_INT)
   #define SIZEOF_INT 4
   #define SIZEOF_LONG 4
   #define SIZEOF_SHORT 2
   #define SIZEOF_FLOAT 4
   #define SIZEOF_DOUBLE 8
#endif

typedef signed char kint8;
typedef unsigned char kuint8;

#if SIZEOF_INT == 4
    typedef int kint32;
    typedef unsigned int kuint32;
#elif SIZEOF_LONG == 4
    typedef long int kint32;
    typedef unsigned int kuint32;
#endif

#if SIZEOF_SHORT == 2
    typedef short int kint16;
    typedef unsigned short int kuint16;
#elif SIZEOF_INT == 2
    typedef int kint16;
    typedef unsigned int kuint16;
#endif

#if SIZEOF_FLOAT == 4
    typedef float kfloat32;
#endif
#if SIZEOF_DOUBLE == 8
    typedef double kfloat64;
#endif


inline const char* fileBasename (const char* const filename)
{
  const char* p = strrchr (filename, '/');
  return (p ? p + 1 : filename);
}


/* strfuncs.cpp */
char* str_skip_head(const char* str, const char* const charlist);
char* str_skip_head(const char* str, char* charlist);
char *str_lower(char *s);
char *str_wrm_tail(char *str);
char *str_rm_tail(char *str, const char* const charlist);
char *str_upper(char *str);

/* syserror.cpp */
void sys_error(int severity, const char *msg, ...);
void sys_verror (std::string& strOutput, int severity, const char *msg, va_list arg);
void sys_error_level(int severity);
extern unsigned long int g_lSysErrorMaxCount;

// Math Section

#include <cmath>

#define PI      3.14159265358979323846
#define HALFPI  1.57079632679489661923  /* PI divided by 2 */
#define QUARTPI 0.78539816339744830962  /* PI divided by 4 */
#define I_PI    0.31830988618379067154  /* Inverse of PI */
#define I_PID2  0.63661977236758134308  /* Inverse of PID2 */

#define TWOPI   6.28318530717958647692
#define SQRT2   1.414213562373095049

#define F_EPSILON       1.0E-6
#define D_EPSILON       1.0E-10

#define ASSUMEDZERO  1E-10


inline double convertDegreesToRadians (double x)
{ 
	return (x * (PI/180.)); 
}

inline double convertRadiansToDegrees (double x)
{ 
	return (x*(180./PI)); 
}

template<class T>
inline T nearest (double x)
{ 
	return (x > 0 ? static_cast<T>(x+0.5) : static_cast<T>(x-0.5)); 
}

template<class T>
inline T maxValue (T x, T y)
{ 
	return (x > y ? x : y); 
}

inline bool isEven (int n)
{
	return (n % 2) == 0; 
}

inline bool isOdd (int n)
{ return (n % 2) != 0; }



inline int imax (int a, int b)
{ 
	return (a >= b ? a : b); 
}


inline double dmax (double a, double b)
{ 
	return (a >= b ? a : b); 
}

template<class T>
inline T clamp (T value, T lowerBounds, T upperBounds)
{ 
	return (value >= upperBounds ? upperBounds : (value <= lowerBounds ? lowerBounds : value )); 
}

template<class T>
inline T lineLength (T x1, T y1, T x2, T y2)
{ 
	return static_cast<T>( sqrt ((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) ); 
}

template<class T>
inline void minmax_array (const T* array, const int n, T& min, T& max)
{
  max = min = array[0];

  for (int i = 1; i < n; i++)
    if (array[i] < min)
      min = array[i];
    else if (array[i] > max)
      max = array[i];
}


//////////////////////////////////////////////////////////////
// FUNTION DECLARATIONS
//////////////////////////////////////////////////////////////

// mathfuncs.cpp
double normalizeAngle (double theta);
double integrateSimpson (const double xmin, const double xmax, const double *y, const int np);
void vectorNumericStatistics (std::vector<double> vec, const int nPoints, double& min, double& max, double& mean, double& mode, double& median, double& stddev);


#endif  /* #ifndef CTSUPPORT_H */
