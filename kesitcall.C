

#include "pr1.h"
#include "prgl.h"
#include "section.h"
#include "glutMaster.h"
#include "glutWindow.h"
#include <fstream>
GlutMaster  * glutMaster;
MeshWindow  * firstWindow = 0;

using namespace std;
void gl_window_supermesh(supermesh &s);


int main(int argc, char** argv) {
  


//  supermesh s;
  // triple_rein_rect(s,300,500, 
  //                  3,8,35, 
  //                  3,8,465,
  //                  2,8,250);
//  double_rein_rect(s,300,500, 5,8,80,4,8,420);
//  single_rein_rect(s,300,500,5,8,80);
//  duck_section(s);



//  interaction(s,500, "dck");
  char cucu[1024];
  for (unsigned int i = 4; i<10; i++){
    supermesh *s = new supermesh;
    single_rein_rect(*s,300,500,5,i,80);
    sprintf(cucu,"l%02d",i);
    interaction(*s,500, cucu);
    delete s;
  }



//  gl_window_supermesh(*s);

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

  return 0;

}


void gl_window_supermesh(supermesh &s){
  glutMaster   = new GlutMaster();    
  
  firstWindow  = new MeshWindow(glutMaster,
                                700, 700,    // height, width
                                0, 0,    // initPosition (x,y)
                                "Poot"); // title
  firstWindow->glsupermesh = &s;

  firstWindow->init();


  glutMaster->CallGlutMainLoop();

}
