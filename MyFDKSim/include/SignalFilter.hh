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
 * Last modified: 2015-10-22 09:17:1445476623
 * Filename     : SignalFilter.hh
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#ifndef SignalFilter_h
#define SignalFilter_h 1

#ifndef F_EPSILON
#define F_EPSILON 1E-6
#endif
#ifndef PI
#define PI 3.14159265
#endif

#include <math.h>

class SignalFilter
{
  private:
    int m_FilterID;
    int m_DomainID;
    int m_FilterPoints;
    double m_FilterParams;
    double m_Bandwidth;
    double m_FilterMin;
    double m_FilterMax;
    double m_FilterIncrement;
    double* m_Filter;
  public:
    SignalFilter(int, double, double, int, double, double, int);
    ~SignalFilter();
  public:
    int GetFilterID() { return m_FilterID; }
    int GetDomainID() { return m_DomainID; }
    int GetSignalPoints() { return m_FilterPoints; }
    double GetBandwidth() { return m_Bandwidth; }
    
  public:
    static const int DOMAIN_FREQUENCY;
    static const int DOMAIN_SPATIAL;
    static const int FILTER_ABS_COSINE;
    static const int N_INTEGRAL;
  private:
    void CreateFrequencyFilter();
    void CreateSpatialFilter();
    bool HaveAnalyticSpatial();
    double FrequencyResponse(double);
    double SpatialResponseAnalytic(double);
    double SpatialResponseCalc(double);
    double integrateSimpson(double,double,double*,int);
    double integral_abscos(double u,double w)
    { return (fabs(u)>F_EPSILON?(cos(u*w)-1)/(u*u)+w/u*sin(u*w):(w*w/2)); }
};

#endif
