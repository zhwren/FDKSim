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
 * Last modified: 2015-10-21 16:55:1445417737
 * Filename     : Projections.hh
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#ifndef Projections_h
#define Projections_h 1

class ScannerGeometry;
class Projections
{
  private:
    ScannerGeometry* m_Scanner;
    double* m_Data;
  public:
    Projections(ScannerGeometry*,double*);
    Projections(Projections&);
    ~Projections();
  public:
    ScannerGeometry* GetScannerGeometry() { return m_Scanner; }
    double* GetProjectionData() { return m_Data; }
};

#endif
