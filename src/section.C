#include "section.h"

void single_rein_rect(supermesh &s, Real width, Real height, 
                 unsigned int nrein, Real rad, Real dprime)
{
  int nsides=12;


  s.meshes = list<mesh> (2);

  rectangular_section(s.meshes.front(),-1*width/2,width/2,height,0);
  lateral_reinforcements(s.meshes.back(), nrein, dprime, width,rad,nsides);

  s.meshes.front().subtract(s.meshes.back());

// Triangulation  
  s.meshes.front().triangulate_mesh(500);//50
  s.meshes.back().triangulate_mesh(50);//5

// Materials
  s.meshes.front().mat = new constitutive::concrete1;
  s.meshes.back().mat = new constitutive::steel1;

// Print info
  s.meshes.front().print_info(); cout<<endl;
  s.meshes.back().print_info();

// Color
  s.meshes.front().color_mesh(0,0,4);
  s.meshes.back().color_mesh(0,0,1);


}


void double_rein_rect(supermesh &s, Real width, Real height,
                 unsigned int nrein, Real rad, Real dprime, Real dprpr)
{
  int nsides=12;

  mesh dummy;
  s.meshes = list<mesh> (2);

  rectangular_section(s.meshes.front(),-1*width/2,width/2,height,0);
  lateral_reinforcements(s.meshes.back(), nrein, dprime, width,rad,nsides);
  lateral_reinforcements(dummy, nrein, height - dprpr, width,rad,nsides);
  
  s.meshes.back().add ( dummy );

  s.meshes.front().subtract(s.meshes.back());

// Triangulation  
  s.meshes.front().triangulate_mesh(500);//50
  s.meshes.back().triangulate_mesh(50);//5

// Materials
  s.meshes.front().mat = new constitutive::concrete1;
  s.meshes.back().mat = new constitutive::steel1;

// Print info
  s.meshes.front().print_info(); cout<<endl;
  s.meshes.back().print_info();

// Color
  s.meshes.front().color_mesh(0,0,4);
  s.meshes.back().color_mesh(0,0,1);


}

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
    center[0]=-1*width/2 + (i+1)*width/(n+1) ;
    center[1]=y;
    circular_section(*mesh_it,radius,center,nsides);

  }

  for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++){
    m.add(*mesh_it);
  }

//  return mesh(result);
}




// void single_rein(supermesh &s, Real width, Real height, unsigned int nrein, Real rad, Real dprime){
//   int nsides=12;
// //  std::list<mesh> m(1);
//   s.meshes.push_back(mesh());

//   std::list<mesh> rein(nrein);
//   std::list<mesh>::iterator mesh_it ;
//    unsigned int i;
//    Real center[2];
//   rectangular_section(s.meshes.front(),-1*width/2,width/2,height,0);

//   i=0;
//   for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++,i++){
//     center[0]=-1*width/2 + (i+1)*width/(nrein+1) ;
//     center[1]=dprime;
//     circular_section(*mesh_it,rad,center,nsides);

//   }


//   // for (unsigned int i=0; i<count; i++){
//   //   center[0]=-1*width/2 + (i+1)*width/(count+1) ;
//   //   center[1]=dprime;
//   //   circular_section(rein[i],rad,center,nsides);
//   // }
// //  Real center[] = {80,80}; circular_section(m2,8,center,12);


  

//   for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++){
//     s.meshes.front().subtract(*mesh_it);
//   }

//   s.meshes.front().triangulate_mesh(500);//50



//   s.meshes.push_back(mesh());
//   for (mesh_it = rein.begin(); mesh_it != rein.end(); mesh_it++){
//     s.meshes.back().add(*mesh_it);
//   }
//   s.meshes.back().triangulate_mesh(50);//5

//   s.meshes.front().mat = new constitutive::concrete1;
//   s.meshes.back().mat = new constitutive::steel1;


//   s.meshes.front().print_info(); cout<<endl;
//   s.meshes.back().print_info();

// // Color
//   s.meshes.front().color_mesh(0,0,4);
//   s.meshes.back().color_mesh(0,0,1);


// }
