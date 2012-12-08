

#include "pr1.h"
#include "prgl.h"
#include "section.h"
#include "glutMaster.h"
#include "glutWindow.h"
#include <fstream>
GlutMaster  * glutMaster;
MeshWindow  * firstWindow = 0;

using namespace std;

// a concrete class that allows cracking: \eps_ultimate>\eps_cracking
class concrete2: public material{
public:
  Real sig(Real eps);
  concrete2(Real a, Real b, Real c, Real d, Real e, Real f);
  concrete2();
private:
  void init(Real f_ck, Real eps_0, Real eps_cr1, Real eps_cr2, 
            Real eps_cu1, Real eps_cu2);
};


class trilinear: public material{
public:
  Real sig(Real eps);
  trilinear(Real a, Real b, Real c, Real d, Real e);
  trilinear();
private:
  void init(Real E, Real E_sh, Real eps_y, Real eps_sh, Real eps_u);
};


struct reinf
{
  Real x[2];
  Real d;
  reinf(Real a, Real b, Real dia){x[0] = a; x[1] = b; d=dia;}
};

void tsect_rein(supermesh &s, Real width, Real height, Real a, Real b, 
                vector<reinf> &r);

void gl_window_supermesh(supermesh &s);


int main(int argc, char** argv) {
  // reinforcement diameter
  Real diameter = 28;
  // create the reinforcement array
  vector<reinf> myrein;
  myrein.push_back(reinf(-110,50,diameter));
  myrein.push_back(reinf(110 ,50,diameter));
  myrein.push_back(reinf(0   ,50,diameter));
  myrein.push_back(reinf(110 ,113,diameter));
  myrein.push_back(reinf(-110,113,diameter));
  
  supermesh s;
  // create the t-section
  tsect_rein(s,300,600,150,1800,myrein);
  // readjust the center of gravity as in the question
  s.cg = node(0,75);

  m_vs_phi(s,1000,0., "lol",true);
  interaction(s,500, "lol",true);

  gl_window_supermesh(s);

  return 0;

}


void gl_window_supermesh(supermesh &s){
  glutMaster   = new GlutMaster();    
  firstWindow  = new MeshWindow(glutMaster, 700, 700, 0, 0, "tsect"); 
  firstWindow->glsupermesh = &s;
  firstWindow->init();
  glutMaster->CallGlutMainLoop();
}

void tsect_rein(supermesh &s, Real width, Real height, Real a, Real b, 
                vector<reinf> &r)
{
  const int nsides=12;

  Real side_clearance = width/3;

  s.meshes = list<mesh> (2);

  rectangular_section(s.meshes.front(),-1*width/2,width/2,height-a,0);
  mesh dummy;
  rectangular_section(dummy,-1*b/2,b/2,height,height-a);
  s.meshes.front().add(dummy);

  for (vector<reinf>::iterator it=r.begin(); it!=r.end(); it++){
    mesh lol;
    circular_section(lol,it->d/2,it->x,nsides);
    s.meshes.back().add(lol);
  }


  s.meshes.front().subtract(s.meshes.back());

// Triangulation  
  s.meshes.front().triangulate_mesh(800);//50
  s.meshes.back().triangulate_mesh(20);//5

// Materials
  s.meshes.front().mat = 
    new concrete2(25*0.85,
                  2e-3,
                  -0.0035, 
                  0.00015,
                  -0.020,
                  0.080);
  s.meshes.back().mat = new trilinear;

// Print info
  s.meshes.front().print_info(); cout<<endl;
  s.meshes.back().print_info();

// Color
  s.meshes.front().color_mesh(0,0,4);
  s.meshes.back().color_mesh(0,0,1);
  s.calc_cg();


}


///////////////////////
// Concrete2 Follows //
///////////////////////

Real concrete2::sig(Real eps){
  if (eps < par[2] || eps > par[3]){
    return 0.;
  }
  Real result;
  if (eps < 0){
    result = -1*par[0] * ((2*-1*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
  }else{
    result = eps*par[6];
  }
  return result;
}

concrete2::concrete2(Real f_ck, Real eps_0, Real eps_cr1, Real eps_cr2, 
                     Real eps_cu1, Real eps_cu2)
{
  init(f_ck, eps_0, eps_cr1, eps_cr2, eps_cu1, eps_cu2);

}

concrete2::concrete2(){
  init(16,2e-3,-0.0035, 0.00015,-0.0035,0.00015);
}

void concrete2::init(Real f_ck, Real eps_0, Real eps_cr1, Real eps_cr2, 
                     Real eps_cu1, Real eps_cu2){
  par = std::vector<Real> (7);
  par[0]=(Real) fabs(f_ck);
  par[1]=(Real) eps_0;
  par[2]=(Real) -1*fabs(eps_cr1); // compresssion cracking
  par[3]=(Real)    fabs(eps_cr2); // tension cracking
  par[4] = (Real) -1*fabs(eps_cu1); // compression ultimate
  par[5] = (Real)    fabs(eps_cu2); // tension ultimate
  epsult[0] = par[4];
  epsult[1] = par[5];
  description="Concrete 2";
  epsulttol[0] = epsult[0]-TOLERANCE; 
  epsulttol[1] = epsult[1]+TOLERANCE;
  par[6]=0.35*sqrt(f_ck)/eps_cr2; // tension slope
}

///////////////////////
// Steel2 Follows    //
///////////////////////


Real trilinear::sig(Real eps){
  if (eps < epsulttol[0] || eps > epsulttol[1])
    return 0.;
  Real result;
  Real epsp = fabs(eps);
  Real pos = eps > 0 ? 1. : -1.;

  if (epsp < par[2]){
    result = pos*epsp*par[0];

  }else if (epsp > par[2] && epsp < par[3]){
    result = pos*par[2]*par[0];
        
  }else if (epsp > par[3] && epsp < par[4]){
    result = pos*(par[1]*(epsp-par[3]) + par[2]*par[0]);
  }else{ result = 0;}
  return result;
}


trilinear::trilinear(Real E, Real E_sh, Real eps_y, Real eps_sh, Real eps_u){
  init(E, E_sh, eps_y, eps_sh, eps_u);
}

trilinear::trilinear(){
  init(2e5, 2e5*0.005, 0.001825, 0.01, 0.15);
}

void trilinear::init(Real E, Real E_sh, Real eps_y, Real eps_sh, Real eps_u){
  par = std::vector<Real> (5);
  par[0]=(Real) E;
  par[1]=(Real) E_sh;
//    par[2]=(Real) eps_cu;
  par[2] = (Real) fabs(eps_y);
  par[3] = (Real) fabs(eps_sh); 
  par[4] = (Real) fabs(eps_u); 

  epsult[0] = (Real) -1*fabs(eps_u); 
  epsult[1] = (Real) fabs(eps_u); 
  description = "Trilinear steel";
  epsulttol[0] = epsult[0]-TOLERANCE; 
  epsulttol[1] = epsult[1]+TOLERANCE;
}
