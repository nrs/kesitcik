

#include "pr1.h"
#include "prgl.h"
#include "section.h"
#include "glutMaster.h"
#include "glutWindow.h"
#include <fstream>
GlutMaster  * glutMaster;
MeshWindow  * firstWindow = 0;

using namespace std;


void lolfun( supermesh &s, Real aa, Real bb, Real slope, char *lolchar );
vector<Real> lol2(mesh &m, Real a, Real b);


int main(int argc, char** argv) {
  


/*  
  mesh m1,m2;

  rectangular_section(m1,-150,150,500,0);
  Real center[] = {80,80}; circular_section(m2,8,center,12);

  m1.subtract(m2);

  m1.print_info();
  m1.triangulate_mesh(0.1);
*/
  supermesh s;
  triple_rein_rect(s,300,500, 
                   3,8,35, 
                   3,8,465,
                   2,8,250);
//  double_rein_rect(s,300,500, 5,8,80,4,8,420);
//  single_rein_rect(s,300,500,5,8,80);

  char cucu[1024];
  Real slope=-1e-6;
  Real bb=0;
//  lolfun(s,-0.00350,497.363, -5.39057e-05 ,"lol.vtk");

//  interaction(s,-3e6,3e6,200);
  interaction(s,-8e8,8e8,200);


#if 0
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







void lolfun( supermesh &s, Real aa, Real bb, Real slope, char *lolchar ){

  list<mesh>::iterator mesh_it;
  std::list<node>::iterator node_it;
  std::list<triangle>::iterator tri_it;

  vector<vector<Real> > luli;
  // for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++)
  //   mesh_it->init_node2tri();


//  std::cout << s.meshes.front().mat->sig(0.00001) << std::endl;
  
#if 0
  std::vector<Real> sig, eps;

  eps = numtk::range(-4e-2,4e-2,1000);
  for (unsigned int i=0; i<eps.size();i++)
    sig.push_back(  s.meshes.front().mat->sig(eps[i]) );

  numtk::output_gp(eps,sig,"lolfile.txt");
  node a(0,2);
//  std::cout << numtk::distance(a,0.3,1) << std::endl;
#endif


  Real F = 0;
  unsigned int i = 0;
//  Real aa=-0.00002;
//  Real bb=250.;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
    vector<Real> epsv;
    vector<Real> sigv;
//    Real dum2;
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++){

      epsv.push_back(  slope* ( (*(*tri_it).cent)(1) - bb) +aa);
//      dum2 = m.mat->sig( a* ( (*(*tri_it).cent)(1) - b) );
//    dum = 0.001;
//    cout << m.mat->description << endl;
//      result.push_back(dum2);
//    cout<< result.back() << endl;
    }
    i=0;
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++, i++){


      sigv.push_back( mesh_it->mat->sig(epsv[i]) );
      F+=tri_it->area*sigv.back();

    }

//    vector<Real> dum=lol2(*mesh_it,-0.00002,cort);


    luli.push_back(project2nodes(*mesh_it,sigv));
//    luli.push_back(project2nodes(*mesh_it,sigv));
  }

  printf(" a = %.2e  b = %.2e  F = %.2e\n",aa,bb,F);




  add_line(node(-200,bb), node(200,bb),
           s.miscnodes, s.misclines);
  
  set_color(s.misclines,3);

  out_vtk1(s.meshes,luli, lolchar);


  // for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
  //   for (node_it = mesh_it->nodes.begin(); node_it!=mesh_it->nodes.end(); node_it++){
      
      
  //   }    
  // }
  

  
//   for (vector<vector<Real> >::iterator vec_it = luli.begin(); 
//        vec_it != luli.end(); vec_it++){
//     for (vector<Real>::iterator it = vec_it->begin();
//          it != vec_it->end(); it++){
// //      cout<< *it << endl;
//     }
//   }


}

std::vector<Real> lol2(mesh &m, Real a, Real b){
  std::list<triangle>::iterator tri_it;
  std::vector<Real> result;
  Real dum;
  for (tri_it = m.triangles.begin(); tri_it!=m.triangles.end(); tri_it++){
    dum = m.mat->sig( a* ( (*(*tri_it).cent)(1) - b) );
//    dum = 0.001;
//    cout << m.mat->description << endl;
    result.push_back(dum);
//    cout<< result.back() << endl;
  }

  //Get the vertical bounding limits of the meshes
  
  node *n[2];
//  nodes.push_back


  return result;
}


