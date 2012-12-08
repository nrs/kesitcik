#include "pr1.h"

// a concrete class that does not allow cracking: \eps_ultimate=\eps_cracking
  class concrete1: public material{
  public:
    Real sig(Real eps);
    concrete1(Real a, Real b, Real c, Real d);
    concrete1();
  private:
    void init(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2);
  };

  class steel1: public material{
  public:
    Real sig(Real eps);
    steel1(Real a, Real b, Real c, Real d, Real e, Real f);
    steel1();
  private:
    void init(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2,
              Real eps_u1, Real eps_u2);  
  };


void lateral_reinforcements(mesh &m, unsigned int n, Real y, Real width, 
                                  Real radius, unsigned int nsides)
{
  std::list<mesh> rein(n);
  m.clear();
//  mesh result;
  std::list<mesh>::iterator mesh_it ;
  unsigned int i;
  Real center[2];

  i=0;
  for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++,i++){
//    center[0]=-1*width/2 + (i+1)*width/(n+1) ;
    center[0]=-1*width/2 + (i)*width/(n-1);
    center[1]=y;
    circular_section(*mesh_it,radius,center,nsides);

  }

  for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++){
    m.add(*mesh_it);
  }

//  return mesh(result);
}


void single_rein_rect(supermesh &s, Real width, Real height, 
                 unsigned int nrein, Real rad, Real dprime)
{
  const int nsides=12;

  Real side_clearance = width/3;

  s.meshes = list<mesh> (2);

  rectangular_section(s.meshes.front(),-1*width/2,width/2,height,0);
  lateral_reinforcements(s.meshes.back(), nrein, dprime, width-side_clearance,rad,nsides);

  s.meshes.front().subtract(s.meshes.back());

// Triangulation  
  s.meshes.front().triangulate_mesh(500);//50
  s.meshes.back().triangulate_mesh(50);//5

// Materials
  // s.meshes.front().mat=auto_ptr<material>(new concrete1);
  // s.meshes.back().mat=auto_ptr<material>(new steel1);
  s.meshes.front().mat = new concrete1;
  s.meshes.back().mat = new steel1;

// Print info
  s.meshes.front().print_info(); cout<<endl;
  s.meshes.back().print_info();

// Color
  s.meshes.front().color_mesh(0,0,4);
  s.meshes.back().color_mesh(0,0,1);
  s.calc_cg();


}



///////////////////////
// Concrete1 Follows //
///////////////////////

Real concrete1::sig(Real eps){
  if (eps < epsulttol[0] || eps > epsulttol[1]){
    return 0.;
  }
  Real result;
  if (eps < 0){
    result = -1*par[0] * ((2*-1*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
  }else{
    result = 1.6 * ((2*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
  }
  return result;
}
concrete1::concrete1(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2){
  init(f_ck, eps_0, eps_cu1, eps_cu2);
}

concrete1::concrete1(){
  init(16,2e-3,-0.003, 0.003);

//    init(20,2e-3,-3.5e-3, 1e-3);
}

void concrete1::init(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2){
  par = std::vector<Real> (4);
  par[0]=(Real) f_ck;
  par[1]=(Real) eps_0;
//    par[2]=(Real) eps_cu;
  par[2] = (Real) -1*fabs(eps_cu1); // compression ultimate
  par[3] = (Real)    fabs(eps_cu2); // tension ultimate
  epsult[0] = par[2];
  epsult[1] = par[3];
  description="Concrete 1";
  epsulttol[0] = epsult[0]-TOLERANCE; 
  epsulttol[1] = epsult[1]+TOLERANCE;
}


///////////////////////
// Steel1 Follows    //
///////////////////////


Real steel1::sig(Real eps){
//    if (eps > par[2]) throw 0;
  if (eps < epsulttol[0] || eps > epsulttol[1])
    return 0.;
  Real result;
  if (eps < 0){
    if (eps > par[2])
      result = eps*par[0];
    else
      result = par[2]*par[0];

  }else{
    if (eps < par[3])
      result = eps*par[0];
    else
      result = par[3]*par[0];
  }
  return result;
}
steel1::steel1(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2,
               Real eps_u1, Real eps_u2){
  init(f_ck, eps_0, eps_cu1, eps_cu2, eps_u1, eps_u2);
}

steel1::steel1(){
//    init(200000,0.0025,-3.6e-3, 1.1e-3);
  init(200000,0.0025,-0.001825, 0.001825,-0.35,0.35);

}

void steel1::init(Real f_ck, Real eps_0, Real eps_sy1, Real eps_sy2,
                  Real eps_u1, Real eps_u2){
  par = std::vector<Real> (4);
  par[0]=(Real) f_ck;
  par[1]=(Real) eps_0;
//    par[2]=(Real) eps_cu;
  par[2] = (Real) -1*fabs(eps_sy1); // pos yield
  par[3] = (Real)    fabs(eps_sy2); // neg yield
  epsult[0] = eps_u1; // pos ultimate
  epsult[1] = eps_u2; // neg ultimate
  description = "Steel 1";
  epsulttol[0] = epsult[0]-TOLERANCE; 
  epsulttol[1] = epsult[1]+TOLERANCE;
}
