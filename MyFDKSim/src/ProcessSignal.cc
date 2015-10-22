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
 * Last modified: 2015-10-22 10:14:1445480049
 * Filename     : ProcessSignal.cc
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#include "ProcessSignal.hh"
#include "SignalFilter.hh"
#include "ScannerGeometry.hh"

const int ProcessSignal::FILTER_METHOD_FFTW = 4;
const int ProcessSignal::FILTER_GENERATION_INVERSE_FOURIER = 1;

ProcessSignal::ProcessSignal(int zeropad, int polationFactor, SignalFilter* filter,
    ScannerGeometry* scanner)
{
  m_FilterID = filter->GetFilterID();
  m_DomainID = filter->GetDomainID();
  m_FilterParam = filter->GetPilterParam();
  m_FilterMethod = FilterMethod;
  m_FilterGeneration = FilterGeneration;

  m_FocalLength = scanner->GetFocalLength();
  m_SourceToDetectorLength = scanner->GetSourceToDetectorLength();
  m_SignalIncrement = scanner->GetDetectorColumnAngleSpacing();
  m_Bandwidth = 1./m_SignalIncrement;
  i
