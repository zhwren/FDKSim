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
 * Last modified: 2015-10-22 10:10:1445479823
 * Filename     : ProcessSignal.hh
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#ifndef ProcessSignal_h
#define ProcessSignal_h 1

#include <complex>
#include <fftw3.h>

class SignalFilter;
class ScannerGeometry;

class ProcessSignal
{
  private:
  public:
    static const int FILTER_METHOD_FFTW;
    static const int FILTER_GENERATION_INVERSE_FOURIER;
    enum { FORWARD = -1,BACKWARD = 1 };
  public:
    ProcessSignal(int,int,SignalFilter*,ScannerGeometry*);
    ~ProcessSignal();
};
