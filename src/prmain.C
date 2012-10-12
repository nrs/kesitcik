

#include "pr1.h"
#include "prgl.h"
#include "section.h"
#include "glutMaster.h"
#include "glutWindow.h"
#include <fstream>
GlutMaster  * glutMaster;
MeshWindow  * firstWindow = 0;

using namespace std;


int main(int argc, char** argv) {
  


  supermesh s;
  // triple_rein_rect(s,300,500, 
  //                  3,8,35, 
  //                  3,8,465,
  //                  2,8,250);
//  double_rein_rect(s,300,500, 5,8,80,4,8,420);
//  single_rein_rect(s,300,500,5,8,80);
  duck_section(s);

//  lolfun(s,-0.00350,497.363, -5.39057e-05 ,"lol.vtk");


  interaction(s,200, "interaction.txt");


#if 0
  char cucu[1024];
  Real slope=-1e-6;
  Real bb=0;

  for (unsigned int i=0 ; i<50; i++, slope+= -1e-7){
    sprintf(cucu, "a%03d.vtk",i);
    lolfun(s,slope,250,cucu);


  }
#endif

#if 0
  for (unsigned int i=0 ; i<500; i++, bb+= 1){
    Real M,F;
    // sprintf(cucu, "a%03d.vtk",i);
    // lolfun(m,aa,bb,cucu);
    normal_force(s,0,bb,slope);
//    normal_and_moment(m,F,M,0,bb,slope);
//    moment(s,0,bb,slope);
  }
#endif
#if 0
  unsigned int niter = 0;
  Real x[3] = {0.000,-0.00001,-0.00005}; Real fx[2],pivot = 3.5e-5;
  while (true){
    cout << niter++ << " " << x[2] <<endl;

    x[0] = x[1]; x[1] = x[2];
    fx[0]  = normal_force(s,pivot, 250, x[0]);
    fx[1]  = normal_force(s,pivot, 250, x[1]);
    x[2] = x[1] - fx[1] * (x[1] - x[0])/ (fx[1]-fx[0]);

    if (fabs(fx[1] - 0) < 1e-4 || niter > 10) break;
  }
#endif


  glutMaster   = new GlutMaster();    
  
  firstWindow  = new MeshWindow(glutMaster,
                                700, 700,    // height, width
                                0, 0,    // initPosition (x,y)
                                "Poot"); // title
  firstWindow->glsupermesh = &s;

  firstWindow->init();


  glutMaster->CallGlutMainLoop();

}


