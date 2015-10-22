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
 * Last modified: 2015-10-21 16:16:1445415389
 * Filename     : ScannerGeometry.hh
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#ifndef ScannerGeometry_h
#define ScannerGeometry_h 1

#ifndef PI
#define PI 3.14159265
#endif

class ScannerGeometry
{
  private:
    int m_DetectorColumnCount;
    int m_DetectorRowCount;
    int m_TotalRotateNumber;

    double m_DetectorColumnAngleSpacing;
    double m_DetectorRowDistanceSpacing;
    double m_FocalLength;
    double m_SourceToDetectorLength;
    double m_DetectroTocenterLength;

    double m_StartAngle;
    double m_IncrementAngle;
    double m_TotalRotateAngle;

    double m_PitchLength;
    double m_IncrementPitch;
  public:
    ScannerGeometry();
    ~ScannerGeometry();
  public:
    void SetDetectorColumnCount(int i) { m_DetectorColumnCount = i; }
    void SetDetectorRowColumn(int i) { m_DetectorRowCount = i; }
    void SetTotalRotateNumber(int i) { m_TotalRotateNumber = i; }
    void SetStartAngle(double i) { m_StartAngle = i; }
    void SetIncrementAngle(double i) { m_IncrementAngle = i; }
    void SetTotalRotateAngle(double i) { m_TotalRotateAngle = i; }
    void SetPitchLength(double i) { m_PitchLength = i; }
    void SetIncrementPitch(double i) { m_IncrementPitch = i; }
  public:
    int GetDetectorColumnCount() { return m_DetectorColumnCount; }
    int GetDetectorRowCount() { return m_DetectorRowCount; }
    int GetTotalRotateNumber() { return m_TotalRotateNumber; }
    double GetIncrementAngle() { return m_IncrementAngle; }
    double GetDetectorColumnAngleSpacing() { return m_DetectorColumnAngleSpacing; }
    double GetDetectorRowDistanceSpacing() { return m_DetectorRowDistanceSpacing; }
    double GetStartAngle() { return m_StartAngle; }
    double GetPitchLength() { return m_PitchLength; }
    double GetIncrementPitch() { return m_IncrementPitch; }
    double GetFocalLength() { return m_FocalLength; }
    double GetSourceToDetectorLength() { return m_SourceToDetectorLength; }
};

#endif
