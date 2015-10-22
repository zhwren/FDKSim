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
 * Last modified: 2015-10-21 17:30:1445419802
 * Filename     : Correction.cc
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#include "Correction.hh"
#include "ScannerGeometry.hh"
#include <math.h>

Correction::Correction(ScannerGeometry* scanner)
{
  m_detColumn = scanner->GetDetectorColumnCount();
  m_detRow = scanner->GetDetectorRowCount();
  m_detColumnAngle = scanner->GetDetectorColumnAngleSpacing();
  m_detRowLength = scanner->GetDetectorRowDistanceSpacing();
  m_FocalLength = scanner->GetFocalLength();
  m_SourceToDetectorLength = scanner->GetSourceToDetectorLength();
  m_CorrectionMatix = new double[m_detColumn*m_detRow];
  ProcessMatix();
}

Correction::~Correction()
{}

void Correction::ProcessMatix()
{
  for(int column=0; column<m_detColumn; column++)
  {
    double cosa = cos((column-m_detColumn/2+0.5)*m_detColumnAngle);
    for(int row=0; row<m_detRow; row++)
    {
      double temp = (row-m_detRow/2+0.5)*m_detRowLength;
      double cosc = m_SourceToDetectorLength/sqrt(temp*temp + m_SourceToDetectorLength*m_SourceToDetectorLength);
      *(m_CorrectionMatix+row*m_detColumn+column) = cosa*cosc;
    }
  }
}

void Correction::ProcessSignal(double* pSion)
{
  for(int column=0; column<m_detColumn; column++)
    for(int row=0; row<m_detRow; row++)
      *(pSion+m_detColumn*row+column) *= *(m_CorrectionMatix+m_detColumn*row+column);
}
