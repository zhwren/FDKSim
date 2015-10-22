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
 * Last modified: 2015-10-21 16:16:1445415387
 * Filename     : ScannerGeometry.cc
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#include "ScannerGeometry.hh"

ScannerGeometry::ScannerGeometry()
{
  m_DetectorColumnCount = 749;
  m_DetectorRowCount = 1;
  m_DetectorColumnAngleSpacing = 0.1*PI/180;
  m_DetectorRowDistanceSpacing = 1;
  m_FocalLength = 600;
  m_SourceToDetectorLength = 600;
  m_TotalRotateNumber = 360;
  m_StartAngle = -2*PI;
  m_IncrementAngle = -2*PI/m_TotalRotateNumber;
}

ScannerGeometry::~ScannerGeometry()
{}
