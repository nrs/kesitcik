

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

void my_single_rein(supermesh &s, Real width, Real height, 
                    unsigned int nrein, Real rad, Real dprime);


int main(int argc, char** argv) {

  
#if 0
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

  supermesh s;
  my_single_rein(s,300,500,3,8,80);

/*
  m_vs_phi(s,200,0., "l01",true);
  m_vs_phi(s,200,-200e3, "l02",true);
  m_vs_phi(s,200,-400e3, "l03",true);
  m_vs_phi(s,200,-600e3, "l04",true);
  m_vs_phi(s,200,-800e3, "l05",true);
  m_vs_phi(s,200,-1000e3, "l06",true);
*/


  interaction(s,500, "lol",true);


#if 0
  for (unsigned int i = 4; i<10; i++){
    supermesh *s = new supermesh;
    single_rein_rect(*s,300,500,5,i,80);
    sprintf(cucu,"m%02d",i);
    interaction(*s,500, cucu);
    delete s;
  }
#endif

//  gl_window_supermesh(s);




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


