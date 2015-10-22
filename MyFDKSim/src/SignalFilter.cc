 /***********************************************************
 *                         _ooOoo_                          *
 *                        o8888888o                         *
 *                        88" . "88                         *
 *                        (| -_- |)                         *
 *                         O\ = /O                          *
 *                     ____/`---'\____                      *
 *                   .   ' \\| |// `.                       *
 *                    / \\||| : |||// \                     *
 *                  / _||||| -:- |||||- \                   *
 *                    | | \\\ - /// | |                     *
 *                  | \_| ''\---/'' | |                     *
 *                   \ .-\__ `-` ___/-. /                   *
 *                ___`. .' /--.--\ `. . __                  *
 *             ."" '< `.____<|>_/___.' >'"".                *
 *            | | : `- \`.;` _ /`;.`/ - ` : | |             *
 *              \ \ `-. \_ __\ /__ _/ .-` / /               *
 *      ======`-.____`-.___\_____/___.-`____.-'======       *
 *                         `=---='                          *
 *                                                          *
 *      .............................................       *
 *             Buddha bless me, No bug forever              *
 *                                                          *
 ************************************************************
 * Author       : ZhuHaiWen                                 *
 * Email        : zhuhw@ihep.ac.cn/zhwren0211@whu.edu.cn    *
 * Last modified: 2015-10-22 09:20:1445476809
 * Filename     : SignalFilter.cc
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#include "SignalFilter.hh"

const int SignalFilter::DOMAIN_FREQUENCY = 1;
const int SignalFilter::DOMAIN_SPATIAL = 2;
const int SignalFilter::FILTER_ABS_COSINE = 1;
const int SignalFilter::N_INTEGRAL = 500;

SignalFilter::SignalFilter(int idFilter, double dFilterMinimum, double dFilterMaximum,
    int nFilterPoints, double dBandwidth, double dFilterParam, int idDomain)
{
  m_FilterID = idFilter;
  m_DomainID = idDomain;
  m_FilterPoints = nFilterPoints;
  m_FilterParams = dFilterParam;
  m_Bandwidth = dBandwidth;
  m_FilterMin = dFilterMinimum;
  m_FilterMax = dFilterMaximum;
  m_FilterIncrement = (m_FilterMax-m_FilterMin)/(m_FilterPoints-1);
  
  m_Filter = new double[m_FilterPoints];
  if(m_DomainID == DOMAIN_FREQUENCY) CreateFrequencyFilter();
  else if(m_DomainID == DOMAIN_SPATIAL) CreateSpatialFilter();
}

SignalFilter::~SignalFilter()
{}

void SignalFilter::CreateFrequencyFilter()
{
  double x;
  int i;
  for(x=m_FilterMin,i=0;i<m_FilterPoints;i++,x+=m_FilterIncrement)
    m_Filter[i] = FrequencyResponse(x);
}

double SignalFilter::FrequencyResponse(double u)
{
  double q;
  double au = fabs(u);
  double abw = fabs(m_Bandwidth);
  
  switch(m_FilterID)
  {
    case FILTER_ABS_COSINE:
      if(au>=(abw/2)+F_EPSILON) q = 0;
      else if(au<F_EPSILON) q = 1;
      else q = au*cos(PI*au/abw);
      break;
  }
  return q;
}

void SignalFilter::CreateSpatialFilter()
{
  double x = m_FilterMin;
  for(int i=0; i<m_FilterPoints; i++,x+=m_FilterIncrement)
  {
    if( HaveAnalyticSpatial() )
      m_Filter[i] = SpatialResponseAnalytic(x);
    else 
      m_Filter[i] = SpatialResponseCalc(x);
  }
}

bool SignalFilter::HaveAnalyticSpatial()
{
  bool haveAnalytic = false;
  switch(m_FilterID)
  {
    case FILTER_ABS_COSINE:
      haveAnalytic = true;
      break;
    default:
      break;
  }
  return haveAnalytic;
}

double SignalFilter::SpatialResponseAnalytic(double x)
{
  double q,temp;
  double u = 2*PI*x;
  double w = m_Bandwidth/2;
  double b = PI/m_Bandwidth;
  double b2 = 2*PI/m_Bandwidth;
  switch(m_FilterID)
  {
    case FILTER_ABS_COSINE:
      q = integral_abscos(b-u,w) + integral_abscos(b+u,w);
      break;
    default:
      q = 0;
      break;
  }
  return q;
}

double SignalFilter::SpatialResponseCalc(double x)
{
  double zmin,zmax;
  double zinc = (zmax-zmin)/(N_INTEGRAL-1);
  double z = zmin;
  double *q = new double[N_INTEGRAL];
  for(int i=0; i<N_INTEGRAL; i++,z+=zinc)
    q[i] = FrequencyResponse(z)*cos(2*PI*z*x);
  double y = 2*integrateSimpson(zmin,zmax,q,N_INTEGRAL);
  delete q;
  return y;
}

double SignalFilter::integrateSimpson(double xmin,double xmax,double* y,int np)
{
  if (np < 2)
    return (0.);
  else if (np == 2)
    return ((xmax - xmin) * (y[0] + y[1]) / 2);

  double area = 0;
  int nDiv = (np - 1) / 2;  // number of divisions
  double width = (xmax - xmin) / (double) (np - 1);     // width of cells

  for (int i = 1; i <= nDiv; i++) {
    int xr = 2 * i;
    int xl = xr - 2;       // 2 * (i - 1) == 2 * i - 2 == xr - 2
    int xm = xr - 1;       // (xl+xr)/2 == (xr+xr-2)/2 == (2*xr-2)/2 = xr-1

    area += (width / 3.0) * (y[xl] + 4.0 * y[xm] + y[xr]);
  }

  if ((np & 1) == 0)            /* do last trapazoid */
    area += width * (y[np-2] + y[np-1]) / 2;

  return (area);
}
