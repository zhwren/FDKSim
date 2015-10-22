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
 * Last modified: 2015-10-21 15:46:1445413580
 * Filename     : reconstruct.cpp
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#include "ct.h"

Reconstructor::Reconstructor (Projections* rProj, double *pReconData,int filterID, double filt_param,
    int filterMethodID, const int zeropad, int filterGenerationID, int interpID,int interpFactor, Volume vol, 
    ReconstructionROI* pROI, bool bRebinToParallel)
:m_bRebinToParallel(bRebinToParallel)
{
  m_pProjectionDatas = new Projections(rProj->scanner(), rProj->prjData());
  m_pCorrection = new Correction(m_pProjectionDatas->geometry(),m_pProjectionDatas->nDetU(),m_pProjectionDatas->nDetV(), 
  	m_pProjectionDatas->detIncU(), m_pProjectionDatas->detIncV(), m_pProjectionDatas->focalLength(), m_pProjectionDatas->sourceDetectorLength());
  
  double filterBW = 1. / m_pProjectionDatas->detIncU();
  
  m_pProcessSignal = new ProcessSignal (filterID, filterMethodID, filterBW, m_pProjectionDatas->detIncU(),
        m_pProjectionDatas->nDetU(), filt_param, SignalFilter::DOMAIN_SPATIAL, filterGenerationID, zeropad, interpFactor,
        m_pProjectionDatas->geometry(), m_pProjectionDatas->focalLength(), m_pProjectionDatas->sourceDetectorLength());
  m_pBackprojector = new Backprojector (m_pProjectionDatas, pReconData, interpID, interpFactor, vol, pROI);
}

Reconstructor::~Reconstructor ()
{
	if(m_pProjectionDatas != NULL) { delete m_pProjectionDatas; m_pProjectionDatas == NULL;}
	if(m_pBackprojector != NULL) {delete m_pBackprojector; m_pBackprojector = NULL;}
	if(m_pProcessSignal != NULL) {delete m_pProcessSignal; m_pProcessSignal = NULL;}
	if(m_pCorrection != NULL) {delete m_pCorrection; m_pCorrection = NULL;}
}



void Reconstructor::reconstructAllViews ()
{
	reconstructView (0, m_pProjectionDatas->nView());
	postProcessing();
}

void Reconstructor::postProcessing()
{
  m_pBackprojector->PostProcessing();
}

void Reconstructor::reconstructView (int iStartView, int iViewCount)
{
  double* adFilteredProj = new double [m_pProjectionDatas->nDetU() * m_pProjectionDatas->nDetV()];   // filtered projections
  if (iViewCount <= 0)
    iViewCount = m_pProjectionDatas->nView() - iStartView;
  for (int iView = iStartView; iView < (iStartView + iViewCount); iView++)  
  {
    double *pTmp = m_pProjectionDatas->prjData() + (iView * m_pProjectionDatas->nDetU() *m_pProjectionDatas->nDetV());
    m_pCorrection->ProcessSignal(pTmp);
    for(int iv =0; iv < m_pProjectionDatas->nDetV(); iv++)
      m_pProcessSignal->filterSignal(pTmp + iv * m_pProjectionDatas->nDetU(), adFilteredProj + iv * m_pProjectionDatas->nDetU());
    
    m_pBackprojector->BackprojectView(adFilteredProj, iView*m_pProjectionDatas->rotInc() + m_pProjectionDatas->startView());
  }
  
  delete adFilteredProj;
}
