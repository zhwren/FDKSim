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
 * Last modified: 2015-10-21 16:15:1445415335
 * Filename     : main.cc
 * Phone Number : 18625272373                               *
 * Discription  :                                           *
 ***********************************************************/
#include "Projections.hh"
#include "SomeClasses.hh"
#include "ScannerGeometry.hh"

#include <iostream>
#include <cstring>
#include <stdio.h>
using namespace std;

int main()
{
  ScannerGeometry *scanner = new ScannerGeometry();

  int column = scanner->GetDetectorColumnCount();
  int row = scanner->GetDetectorRowCount();
  int nviews = scanner->GetTotalRotateNumber();

  double *pSino = new double[column*row*nviews];
  memset(pSino, 0, sizeof(double)*column*row*nviews);

  FILE *fp = fopen("../data.txt","rb");
  fread(pSino, sizeof(double), column*row*nviews, fp);
  fclose(fp);

  Projections* aProjection = new Projections(scanner, pSino);

  Volume* vol = new Volume(512, 512, 1);
  double* pReconsData = new double[vol->nx()*vol->ny()*vol->nz()];
  memset(pReconsData, 0, sizeof(double)*vol->nx()*vol->ny()*vol->nz());


  return 1;
}
