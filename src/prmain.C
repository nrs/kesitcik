

#include "pr1.h"
#include "prgl.h"
#include "glutMaster.h"
#include "glutWindow.h"
#include <fstream>
GlutMaster  * glutMaster;
MeshWindow  * firstWindow = 0;

using namespace std;

void single_rein(supermesh &s, Real width, Real height, unsigned int count, Real rad, Real dprime);
void lolfun( supermesh &s, Real aa, Real bb, char *lolchar );
vector<Real> lol2(mesh &m, Real a, Real b);


vector<Real> project2nodes( mesh &m, std::vector<Real> &vec);

int main(int argc, char** argv) {
  
  glutMaster   = new GlutMaster();    
  
  firstWindow  = new MeshWindow(glutMaster,
                                700, 700,    // height, width
                                0, 0,    // initPosition (x,y)
                                "Poot"); // title

/*  
  mesh m1,m2;

  rectangular_section(m1,-150,150,500,0);
  Real center[] = {80,80}; circular_section(m2,8,center,12);

  m1.subtract(m2);

  m1.print_info();
  m1.triangulate_mesh(0.1);
*/
  supermesh m;
  single_rein(m,300,500,5,8,80);

  char cucu[1024];
  Real aa=-1e-6;
  Real bb=0;

#if 0
  for (unsigned int i=0 ; i<50; i++, aa+= -1e-7){
    sprintf(cucu, "a%03d.vtk",i);
    lolfun(m,aa,250,cucu);


  }
#endif

#if 1
  for (unsigned int i=0 ; i<500; i++, bb+= 1){
    sprintf(cucu, "a%03d.vtk",i);
    lolfun(m,aa,bb,cucu);

 
  }
#endif


  firstWindow->glsupermesh = &m;

  firstWindow->init();


  m.meshes.front().numbernodes();


  glutMaster->CallGlutMainLoop();

}


void single_rein(supermesh &s, Real width, Real height, unsigned int count, Real rad, Real dprime){
  int nsides=12;
//  std::list<mesh> m(1);
  s.meshes.push_back(mesh());

  std::list<mesh> rein(count);
  std::list<mesh>::iterator mesh_it ;
  unsigned int i;
  Real center[2];
  rectangular_section(s.meshes.front(),-1*width/2,width/2,height,0);

  i=0;
  for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++,i++){
    center[0]=-1*width/2 + (i+1)*width/(count+1) ;
    center[1]=dprime;
    circular_section(*mesh_it,rad,center,nsides);

  }


  // for (unsigned int i=0; i<count; i++){
  //   center[0]=-1*width/2 + (i+1)*width/(count+1) ;
  //   center[1]=dprime;
  //   circular_section(rein[i],rad,center,nsides);
  // }
//  Real center[] = {80,80}; circular_section(m2,8,center,12);
  for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++){
    s.meshes.front().subtract(*mesh_it);
  }

  s.meshes.front().triangulate_mesh(500);//50

//  std::list<node>::iterator node_it;
  // for (node_it=m[0].holes.begin(); node_it != m[0].holes.end(); node_it++){
  //   std::cout << *node_it;
  //   std::cout<<std::endl;
  // }

  s.meshes.push_back(mesh());
  for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++){
    s.meshes.back().add(*mesh_it);
  }
  s.meshes.back().triangulate_mesh(50);//5

  s.meshes.front().mat = new constitutive::concrete1;
  s.meshes.back().mat = new constitutive::steel1;
  s.meshes.front().print_info(); 
  cout<<endl;
  s.meshes.back().print_info();

// A dummy mesh to keep drawings etc.
  s.meshes.front().color_mesh(0,0,4);
  s.meshes.back().color_mesh(0,0,1);


}


void lolfun( supermesh &s, Real aa, Real bb, char *lolchar ){

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

      epsv.push_back(  aa* ( (*(*tri_it).cent)(1) - bb) );
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

  printf(" a = %6f  b = %6f  F = %6f\n",aa,bb,F);




  add_line(node(-200,bb), node(200,bb),
           s.miscnodes, s.misclines);
  
  set_color(s.misclines,3);

//  out_vtk1(s.meshes,luli, lolchar);


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



//  for (unsigned int i=0; i<eps.size();i++)  std::cout << eps[i] << std::endl;;
//  }
//  catch (...) {
//    std::cout << "LOL" <<std::endl;
//  }
//  m.back().print_info();


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


std::vector<Real> project2nodes( mesh &m, std::vector<Real> &vec){
  m.init_node2tri();
  m.numbernodes();
  std::vector<Real> result;
  std::list<node>::iterator node_it;
  std::set<triangle*>::iterator tri_it;
  if (vec.size() != m.centroids.size()) {
    std::cerr << "Projection on nodes unsuccessful." << std::endl;
    return result;
  }
//  cout<< vec.back()<<endl;

  Real dum;
  for (node_it = m.nodes.begin(); node_it!=m.nodes.end(); node_it++){
    dum=0;
    for (tri_it = m.node2tri[&(*node_it)].begin(); 
         tri_it != m.node2tri[&(*node_it)].end();
         tri_it++){
      dum+=vec[(*tri_it)->id];
//      cout << vec[(*tri_it)->id] <<endl;
    }
    if (m.node2tri[&(*node_it)].size() != 0){
      dum /= (Real) m.node2tri[&(*node_it)].size();
    }
//    dum+=node2tri[&(*node_it)];
    result.push_back(dum);
//    cout<< result.back() << endl;
  }
  return result;

}
