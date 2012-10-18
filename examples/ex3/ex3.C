

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
  supermesh s;
  // triple_rein_rect(s,300,500, 
  //                  3,8,35, 
  //                  3,8,465,
  //                  2,8,250);
//  double_rein_rect(s,300,500, 5,8,80,4,8,420);
//  single_rein_rect(s,300,500,5,8,80);
  
  
//  gl_window_supermesh(s);
  
#if 1
  char cucu[1024];
  for (unsigned int i = 4; i<10; i++){
    supermesh *s = new supermesh;
//    single_rein_rect(*s,300,500,5,i,80);
    triple_rein_rect(*s,300,500, 
                     3,i,35, 
                     3,i,465,
                     2,i,250);
    sprintf(cucu,"l%02d",i);
    interaction(*s,500, cucu);
    delete s;
  }
#endif


#if 1
  for (unsigned int i = 4; i<10; i++){
    supermesh *s = new supermesh;
    single_rein_rect(*s,300,500,5,i,80);
    sprintf(cucu,"m%02d",i);
    interaction(*s,500, cucu);
    delete s;
  }
#endif

//  gl_window_supermesh(*s);




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


