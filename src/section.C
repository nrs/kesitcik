#include "section.h"
#include "interval.h"

using namespace constitutive;
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


void double_rein_rect(supermesh &s, Real width, Real height,
                      unsigned int nrein1, Real rad1, Real dprime1,
                      unsigned int nrein2, Real rad2, Real dprime2)
// void double_rein_rect(supermesh &s, Real width, Real height,
//                  unsigned int nrein, Real rad, Real dprime, Real dprpr)
{
  const int nsides=12;
  Real side_clearance = width/3;
  mesh dummy;
  s.meshes = list<mesh> (2);

  rectangular_section(s.meshes.front(),-1*width/2,width/2,height,0);
  lateral_reinforcements(s.meshes.back(), nrein1, dprime1, 
                         width-side_clearance,rad1,nsides);
  lateral_reinforcements(dummy, nrein2, dprime2, 
                         width-side_clearance,rad2,nsides);
  
  s.meshes.back().add ( dummy );

  s.meshes.front().subtract(s.meshes.back());

// Triangulation  
  s.meshes.front().triangulate_mesh(500);//50
  s.meshes.back().triangulate_mesh(50);//5

// Materials
  // s.meshes.front().mat = concrete1();
  // s.meshes.back().mat = steel1();
//  s.meshes.front().mat=auto_ptr<material>(new concrete1);
//  s.meshes.back().mat=auto_ptr<material>(new steel1);
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
void triple_rein_rect(supermesh &s, Real width, Real height,
                      unsigned int nrein1, Real rad1, Real dprime1,
                      unsigned int nrein2, Real rad2, Real dprime2,
                      unsigned int nrein3, Real rad3, Real dprime3)
{
  const int nsides=12;
  Real side_clearance = width/3;
  mesh dummy;
  s.meshes = list<mesh> (2);

  rectangular_section(s.meshes.front(),-1*width/2,width/2,height,0);
  lateral_reinforcements(s.meshes.back(), nrein1, dprime1, 
                         width-side_clearance,rad1,nsides);
  lateral_reinforcements(dummy, nrein2, dprime2, 
                         width-side_clearance,rad2,nsides);
  s.meshes.back().add ( dummy );
  lateral_reinforcements(dummy, nrein3, dprime3,
                         width-side_clearance,rad3,nsides);
  s.meshes.back().add ( dummy );

  s.meshes.front().subtract(s.meshes.back());

// Triangulation  
  s.meshes.front().triangulate_mesh(500);//50
  s.meshes.back().triangulate_mesh(50);//5

// Materials
  s.meshes.front().mat = new concrete1;
  s.meshes.back().mat = new steel1;
//  s.meshes.front().mat=auto_ptr<material>(new concrete1);
//  s.meshes.back().mat=auto_ptr<material>(new steel1);
  // s.meshes.front().mat = concrete1();
  // s.meshes.back().mat = steel1();

// Print info
  s.meshes.front().print_info(); cout<<endl;
  s.meshes.back().print_info();

// Color
  s.meshes.front().color_mesh(0,0,4);
  s.meshes.back().color_mesh(0,0,1);

  s.calc_cg();

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



void duck_section(supermesh &s)
{
  mesh beak, eye;

  s.meshes.clear();

  s.meshes.push_back(mesh());
  s.meshes.back().in_vtk_lines("duck/body.vtk");

  eye.in_vtk_lines("duck/eye.vtk");
  node cgeye(0,0);
  for (list<node>::iterator it = eye.nodes.begin(); 
       it != eye.nodes.end(); it++){
    cgeye+=*it;
  }
  cgeye/=eye.nodes.size();

  cout << cgeye << endl;
  eye.innode.push_back(cgeye);
  s.meshes.back().subtract(eye);
  
  beak.in_vtk_lines("duck/beak.vtk");
  s.meshes.push_back(mesh());
  s.meshes.back().add(eye);
  s.meshes.back().add(beak);

  s.meshes.back().triangulate_mesh(50);
  s.meshes.front().triangulate_mesh(500);

  s.meshes.back().color_mesh(0,0,1);
  s.meshes.front().color_mesh(0,0,4);

  s.meshes.front().mat = new concrete1;
  s.meshes.back().mat = new steel1;
//  s.meshes.front().mat=auto_ptr<material>(new concrete1);
//  s.meshes.back().mat=auto_ptr<material>(new steel1);

  // s.meshes.front().mat = concrete1();
  // s.meshes.back().mat = steel1();

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



Real normal_force( supermesh &s, Real a, Real b, Real slope){
  list<mesh>::iterator mesh_it;
  list<triangle>::iterator tri_it;
  Real F = 0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
    Real eps, sig;
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++){
      eps = slope * ( (*tri_it->cent)(1) - b) + a ;
      sig = mesh_it->mat->sig(eps) ;
      F+=tri_it->area *sig;
    }
  }
//  cout << "dist " << (-1*a/slope + b) << endl;
//  printf(" a= %.2e   b= %.2e   s= %.2e   F= %.2e\n",a,b,slope,F);
  return F;
}


// void normal_and_moment( supermesh &s, Real &F, Real &M, Real a, Real b, Real slope){
//   list<mesh>::iterator mesh_it;
//   list<triangle>::iterator tri_it;
//   F = 0, M = 0;
//   for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
//     Real eps, sigma, y, dist,force;
//     for (tri_it = mesh_it->triangles.begin(); 
//          tri_it!=mesh_it->triangles.end(); tri_it++){
//       y = (*tri_it->cent)(1);
//       eps = slope * (y - b) + a ;
//       sigma = mesh_it->mat->sig(eps) ;
//       dist = y - b;
//       force = tri_it->area * sigma;
//       F += force;
//       M += force * dist;
//     }
//   }
//   printf("a= %.2e   b= %.2e   s= %.2e   F= %.2e   M= %.2e\n",a,b,slope,F,M);
// }


Real moment( supermesh &s, Real a, Real b, Real slope){
  list<mesh>::iterator mesh_it;
  list<triangle>::iterator tri_it;
  Real M = 0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
    Real eps, sigma, y, dist;
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++){
      y = (*tri_it->cent)(1);
//      dist = y - (-1*a/slope + b) ;
      dist = y - s.cg(1);
      eps = slope * (y-b) + a;
//      Real eps2 = slope * (dist-b) + a;
//      if (dist>0  ) 
//      cout << "ASDASDASDASDAS " << dist <<" " << y << " "<< (-1*a/slope + b)  <<endl;
      sigma = mesh_it->mat->sig(eps) ;
      M += tri_it->area * sigma * dist;
    }
  }
//  printf("a= %.2e   b= %.2e   s= %.2e   M= %.2e\n",a,b,slope,M);
  return M;
}





namespace numtk{
  inline Real slope(const Point &a, const Point &b){
    return (a(1) - b(1)) / (a(0) - b(0));
  }

}


vector<Real> root_bisect_force(supermesh &s, Real a, Real b,Real c, Real d, Real force ){
  Real lo = min(c,d);
  Real up = max(c,d);
  vector<Real> initial = numtk::range(lo, up,20);
  vector<Real> result;
  const Real mytolerance = 1e-10;
  const unsigned int iterlim = 100; 

  for (unsigned int i = 0; i<initial.size()-1; i++){
//    cout << a <<" "<< b << " " << initial[i]<<" "<<normal_force(s,a,b,initial[i  ]) << endl;
    if ( numtk::sign(normal_force(s,a,b,initial[i  ])-force)!=
         numtk::sign(normal_force(s,a,b,initial[i+1])-force) ){
      
      Real x=initial[i],z,zold,y=initial[i+1];
      z = zold = (x+y)/2;
      unsigned int niter = 0;
      while (true){
        if (niter>iterlim) break;
        if (numtk::sign(normal_force(s,a,b,z) -force)==
            numtk::sign(normal_force(s,a,b,x) -force))
          x=z; else y=z;
        z = (x+y)/2;

        if (fabs((z-zold)/z) < mytolerance){
//          cout << "#iter = " << niter << " " << z << endl;
          result.push_back(z); break;
        }
        zold=z;
        niter++;
      }
    }
  }
  return result;
}


vector<Real> root_bisect_moment(supermesh &s, Real a, Real b,Real c, Real d, Real mom ){
  Real lo = min(c,d);
  Real up = max(c,d);
  vector<Real> initial = numtk::range(lo, up,100);
  vector<Real> result;
  const Real mytolerance = 1e-10;
  const unsigned int iterlim = 100; 
  
//  for (unsigned int i = 0; i<initial.size(); i++)
//    cout << "PRT "<< initial[i] <<" " << moment(s,a,b,initial[i]) << endl;
    
  for (unsigned int i = 0; i<initial.size()-1; i++){
    Real p =  moment(s,a,b,initial[i  ])-mom;
    Real q =  moment(s,a,b,initial[i+1])-mom;
    if ( numtk::sign(p)!= numtk::sign(q) && 
         !(numtk::isinf(p) || numtk::isinf(q)) ){
      
      Real x=initial[i],z,zold,y=initial[i+1];
      z = zold = (x+y)/2;
      unsigned int niter = 0;
      
      while (true){
//        cout << x << " " << y << endl;
        if (niter>iterlim) {
          cout << "No root." <<endl;
          break;
        }
        if (numtk::sign(moment(s,a,b,z) - mom)==
            numtk::sign(moment(s,a,b,x) -mom))
          x=z; else y=z;
        z = (x+y)/2;

        if (fabs((z-zold)/z) < mytolerance){
          cout << "#iter = " << niter << " " << z << endl;
          cout << "LOL " <<  (-1*a/z + b) << endl;

          result.push_back(z); 
//          cout << moment(s,a,b,z) << " " << mom << endl;
          break;
        }
        zold=z;
        niter++;
      }
    }
  }
  return result;
}

void stress_vtk(supermesh &s, Real a, Real b, Real slope, char *lolchar ){

  list<mesh>::iterator mesh_it;
  std::list<node>::iterator node_it;
  std::list<triangle>::iterator tri_it;

  vector<vector<Real> > lulisig, lulieps;

  unsigned int i = 0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
    vector<Real> epsv;
    vector<Real> sigv;

    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++){
      epsv.push_back(  slope* ( (*(*tri_it).cent)(1) - b) +a );
    }
    i=0;
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++, i++){
      sigv.push_back( mesh_it->mat->sig(epsv[i]) );
    }
    lulisig.push_back(project2nodes(*mesh_it,sigv));
    lulieps.push_back(project2nodes(*mesh_it,epsv));
  }
  printf(" a = %.2e  b = %.2e  s = %.2e\n",a,b,slope);

  out_vtk1(s.meshes, lulisig, lulieps, lolchar);

}




void sample_materials(supermesh &s, Real st, Real en, Real div, char *ofilename){
  std::vector<Real> sig, eps;
  list<mesh>::iterator it;
  char fn[1024]; unsigned int i=0;
  eps = numtk::range(st,en,div);
  for (it = s.meshes.begin(); it != s.meshes.end(); it++,i++){
    sig.clear();
    sprintf(fn, "%s_mat%02d.txt",ofilename, i);
    for (unsigned int i=0; i<eps.size();i++)
      sig.push_back(  it->mat->sig(eps[i]) );
    numtk::output_gp(eps,sig,fn);
  }

}



compoundInterval get_slope_range(Point p, supermesh &s)
{
  list<mesh>::iterator mesh_it;

  compoundInterval all;

  for (mesh_it=s.meshes.begin(); mesh_it!=s.meshes.end(); mesh_it++){
    if (mesh_it->mat == NULL){
      cerr<< "Input still has undefined materials." << endl;
      return all;
    }
  }
  if (s.meshes.size()==0){
    cerr<< "No meshes defined in the input." << endl;
    return all;
  }

  unsigned int i;

  std::list<node>::iterator node_it;
  std::list<triangle>::iterator tri_it;

  //pivot points
  vector<vector<Point> > pivots;

  // the intervals must be either all negative or all positive.
  // pivsign holds that information
  vector<vector<bool> >pivsign;

  // allowable and prohibited intervals
  // represented by a line  between two points
  compoundInterval pro;

  Real ex[2],ep[2];

  i=0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++,i++){
    vector<Point> dum1;
    pivots.push_back(vector<Point>());
    pivsign.push_back(vector<bool>());

    ex[0] = ex[1] = (*mesh_it->triangles.begin()->cent)(1);
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!= mesh_it->triangles.end(); tri_it++){
      (*tri_it->cent)(1) < ex[0] ? ex[0]=(*tri_it->cent)(1) : 0;
      (*tri_it->cent)(1) > ex[1] ? ex[1]=(*tri_it->cent)(1) : 0;
    }
    ep[0] = min(mesh_it->mat->epsult[0],mesh_it->mat->epsult[1]) ;
    ep[1] = max(mesh_it->mat->epsult[0],mesh_it->mat->epsult[1]) ;
    
    Point a(ex[0],ep[0]);
    Point b(ex[1],ep[0]);
    Point c(ex[0],ep[1]);
    Point d(ex[1],ep[1]);

    pivots.back().push_back(a);
    pivots.back().push_back(b);
    pivots.back().push_back(c);
    pivots.back().push_back(d);

    pivsign.back().push_back(1);
    pivsign.back().push_back(0);
    pivsign.back().push_back(0);
    pivsign.back().push_back(1);
    
  }


  vector<compoundInterval> dum;
  for (unsigned int j = 0; j < pivots.size(); j++){
    dum.push_back(compoundInterval());
    for (unsigned int m = 0; m < 2; m++){
      Real l,u;
      l = numtk::slope(p,pivots[j][m]);
      u = numtk::slope(p,pivots[j][m+2]);
//      cout <<i<< ": A " << l << " " << u <<endl;
      if ((numtk::isnan(u) || numtk::isinf(u) || 
           numtk::isnan(l) || numtk::isinf(l)))
        continue;
      dum.back()+=interval(l,u);

    }
//    cout << endl;
  }
  if (dum.size() != 0){
    compoundInterval lol = dum[0];
    for (unsigned int q = 0; q < dum.size(); q++){
//          cout << i << " " << q << "  "<< dum[q] << endl;
      lol = lol.intersection(dum[q]);
//          lol +=(dum[q]);
    } 
    all = lol;
  }

  vector<compoundInterval> dum2;

  for (unsigned int j = 0; j < pivots.size(); j++){
    dum2.push_back(compoundInterval());
    for (unsigned int m = 0; m < 2; m++){
      Real l,u;
      l = numtk::slope(p,pivots[j][m*2]);
      u = numtk::slope(p,pivots[j][m*2+1]);

      if (numtk::isnan(u) && numtk::isnan(l)) continue;
      if (numtk::isinf(u) && numtk::isinf(l)) continue;
//      cout <<i<< ": U " << l << " " << u <<endl;
        
      if (numtk::isinf(l)){
        bool dumb = numtk::isposit(pivots[j][m*2](0)-p(0)) *
          !numtk::isposit(u);
        all.subtract_infinity(u,!dumb);
        continue;
      }
      if (numtk::isinf(u)){
        bool dumb = numtk::isposit(pivots[j][m*2+1](0)-p(0))*
          !numtk::isposit(l);
        all.subtract_infinity(l,!dumb);
        continue;
      }
      // if (numtk::isnan(l)){
      //   all.subtract_infinity(0,!pivsign[i][k]);
      //   continue;
      // }
      // if (numtk::isnan(u)){
      //   all.subtract_infinity(0,!pivsign[i][k]);
      //   continue;
      // }
      all-=interval(u,l);
    }
  }
  
  return all;
}




void m_vs_phi(supermesh &s, unsigned int ndiv, Real N,
              char *ofilename, bool plotothers)
{
  vector<vector<Real> > result;
  list<mesh>::iterator mesh_it;
  char dumchar[1024];

//  if (plotothers) sample_materials(s,-5e-3,5e-3,1000,ofilename);

  for (mesh_it=s.meshes.begin(); mesh_it!=s.meshes.end(); mesh_it++){
    if (mesh_it->mat == NULL){
      cerr<< "Input still has undefined materials." << endl;
      return;
    }
  }
  if (s.meshes.size()==0){
    cerr<< "No meshes defined in the input." << endl;
    return;
  }

  unsigned int i;

  std::list<node>::iterator node_it;
  std::list<triangle>::iterator tri_it;

  //pivot points
  vector<vector<Point> > pivots;

  //The vector that holds the solution, namely:
  vector<vector<Real> > pivslopes;

  // the intervals must be either all negative or all positive.
  // pivsign holds that information
  vector<vector<bool> >pivsign;

  // allowable and prohibited intervals
  // represented by a line  between two points
  vector<vector<compoundInterval> > all;
  vector<vector<compoundInterval> > pro;

  // pivot compound intervals:
  vector<vector<compoundInterval> >pivint;


  Real ex[2],ep[2];

  i=0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++,i++){
    vector<Point> dum1;
    pivots.push_back(vector<Point>());
    pivsign.push_back(vector<bool>());

    ex[0] = ex[1] = (*mesh_it->triangles.begin()->cent)(1);
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!= mesh_it->triangles.end(); tri_it++){
      (*tri_it->cent)(1) < ex[0] ? ex[0]=(*tri_it->cent)(1) : 0;
      (*tri_it->cent)(1) > ex[1] ? ex[1]=(*tri_it->cent)(1) : 0;
    }
    ep[0] = min(mesh_it->mat->epsult[0],mesh_it->mat->epsult[1]) ;
    ep[1] = max(mesh_it->mat->epsult[0],mesh_it->mat->epsult[1]) ;
    
    Point a(ex[0],ep[0]);
    Point b(ex[1],ep[0]);
    Point c(ex[0],ep[1]);
    Point d(ex[1],ep[1]);

    pivots.back().push_back(a);
    pivots.back().push_back(b);
    pivots.back().push_back(c);
    pivots.back().push_back(d);

    pivsign.back().push_back(1);
    pivsign.back().push_back(0);
    pivsign.back().push_back(0);
    pivsign.back().push_back(1);

    pivint.push_back(vector<compoundInterval>());
    all.push_back(vector<compoundInterval>());
    pro.push_back(vector<compoundInterval>());

    // expand the vectors
    for (unsigned int k = 0; k < 4; k++){
      pivint.back().push_back(compoundInterval());
      all.back().push_back(compoundInterval());
      pro.back().push_back(compoundInterval());
    }
#if 0
    Real l,u; //lower, upper
    for (unsigned int k = 0; k < 4; k++){
      for (unsigned int m = 0; m < 2; m++){
        l = numtk::slope(pivots.back()[k],pivots.back()[m]);
        u = numtk::slope(pivots.back()[k],pivots.back()[m+2]);
        if (numtk::isnan(u) || numtk::isinf(u) || numtk::isnan(l) || numtk::isinf(l))
          continue;
        all.back()[k]+=interval(l,u);
      }                 
        //   all.back()[k]+=;
        // interval(numtk::slope(pivots.back()[k],pivots.back()[1] ),
        //          numtk::slope(pivots.back()[k],pivots.back()[3] ));
    }

    for (unsigned int k = 0; k < 4; k++){
      pro.back()[k]-=
        interval(numtk::slope(pivots.back()[k],pivots.back()[0] ),
                 numtk::slope(pivots.back()[k],pivots.back()[1] ));
      pro.back()[k]-=
        interval(numtk::slope(pivots.back()[k],pivots.back()[2] ),
                 numtk::slope(pivots.back()[k],pivots.back()[3] ));
    }
#endif

    
  }
#if 0
  for (unsigned int i = 0; i < pivots.size(); i++){
    for (unsigned int k = 0; k < 4; k++){ // loop the points
      vector<compoundInterval> dum;
      for (unsigned int j = 0; j < pivots.size(); j++){
        dum.push_back(compoundInterval());
        for (unsigned int m = 0; m < 2; m++){
          Real l,u;
          l = numtk::slope(pivots[i][k],pivots[j][m]);
          u = numtk::slope(pivots[i][k],pivots[j][m+2]);
          cout <<i<< ": A " << l << " " << u <<endl;
          if ((numtk::isnan(u) || numtk::isinf(u) || 
                numtk::isnan(l) || numtk::isinf(l)))
            continue;
          dum.back()+=interval(l,u);

        }
        cout << endl;
      }
      if (dum.size() != 0){
        compoundInterval lol = dum[0];
        for (unsigned int q = 0; q < dum.size(); q++){
//          cout << i << " " << q << "  "<< dum[q] << endl;
          lol = lol.intersection(dum[q]);
//          lol +=(dum[q]);
        } 
        all[i][k] = lol;
      }else continue;

      vector<compoundInterval> dum2;

      for (unsigned int j = 0; j < pivots.size(); j++){
        dum2.push_back(compoundInterval());
        for (unsigned int m = 0; m < 2; m++){
          Real l,u;
          l = numtk::slope(pivots[i][k],pivots[j][m*2]);
          u = numtk::slope(pivots[i][k],pivots[j][m*2+1]);

          if (numtk::isnan(u) && numtk::isnan(l)) continue;
          if (numtk::isinf(u) && numtk::isinf(l)) continue;
          cout <<i<< ": U " << l << " " << u <<endl;
        
          if (numtk::isinf(l)){
            bool dumb = numtk::isposit(pivots[j][m*2](0)-pivots[i][k](0)) *
              !numtk::isposit(u);
            all[i][k].subtract_infinity(u,!dumb);
            continue;
          }
          if (numtk::isinf(u)){
            bool dumb = numtk::isposit(pivots[j][m*2+1](0)-pivots[i][k](0))*
              !numtk::isposit(l);
            all[i][k].subtract_infinity(l,!dumb);
            continue;
          }
          if (numtk::isnan(l)){
            all[i][k].subtract_infinity(0,!pivsign[i][k]);
            continue;
          }
          if (numtk::isnan(u)){
            all[i][k].subtract_infinity(0,!pivsign[i][k]);
            continue;
          }
        
 
          all[i][k]-=interval(u,l);
          // cout <<i<< ": A " << l << " " << u <<endl;
          // if (!(numtk::isnan(u) || numtk::isinf(u) || 
          //       numtk::isnan(l) || numtk::isinf(l)))
          // dum2.back()+=interval(l,u);
        }
      }
    }
  }


  unsigned int j = 0;
  for (unsigned int i=0; i < all.size(); i++)
    for (unsigned int k=0; k < all[i].size(); k++)
      j+=all[i][k].size();
  if (j==0) {
    cerr<< "No solutions possible for this combination." << endl;
    return;
  }


  for (unsigned int i = 0; i < pivots.size(); i++){
    for (unsigned int j = 0; j < pivots[i].size(); j++){
    cout << i << "("<< all[i][j].size() <<") ("
         <<pivots[i][j](0)<<","<<pivots[i][j](1)<< ") ** " << all[i][j];
    cout << endl;
    }
  }
#endif
  for (unsigned int i = 0; i < pivots.size(); i++){
    for (unsigned int j = 0; j < pivots[i].size(); j++){
    cout << i << "("<< all[i][j].size() <<") ("
         <<pivots[i][j](0)<<","<<pivots[i][j](1)<< ") ** "<< endl;
    }
  }

  ofstream of2; 
  sprintf(dumchar, "%s_mphi.txt",ofilename);
  of2.open(dumchar);
  cout << "Writing " << dumchar << "." << endl;


  for (unsigned int i = 0; i < pivots.size(); i++){
    for (unsigned int j = 0; j < 2; j++){
      
      vector<Real> strains = numtk::range(pivots[i][j](1), pivots[i][j+2](1),ndiv);
      cout << " LOL " << pivots[i][j](1) << "  " <<  pivots[i][j+2](1) << endl;
      for (unsigned int k = 0; k < strains.size(); k++){
//        Point p(pivots[i][j](0),strains[k]);
        compoundInterval inter = get_slope_range( Point(pivots[i][j](0),strains[k]) ,s );
//        cout << pivots[i][j](0)<<" "<<  strains[k]<< " "<< inter << endl;

        for (unsigned int l=0; l<inter.size(); l++){
          vector<Real> root = 
            root_bisect_force(s,strains[k], pivots[i][j](0), inter[l].s, inter[l].e, N);
          for (unsigned int q=0; q<root.size(); q++){
            of2 << root[q] << "  " << 
              moment(s,strains[k], pivots[i][j](0),root[q]) << endl;
          }

//          cout << pivots[i][j](0)<<" "<<  strains[k]<< " "<< inter << endl;

//        cout << get_slope_range( Point(pivots[i][j](0),strains[k]) ,s ) << endl;

        }
        
      }
    }
  }
  of2.close();



}












vector<vector<Real> > interaction(supermesh &s, unsigned int ndiv, 
                                  char *ofilename, bool plotothers)
{
  vector<vector<Real> > result;
  list<mesh>::iterator mesh_it;
  char dumchar[1024];

  if (plotothers) sample_materials(s,-5e-3,5e-3,1000,ofilename);

  // for (mesh_it=s.meshes.begin(); mesh_it!=s.meshes.end(); mesh_it++){
  //   if (mesh_it->mat == NULL){
  //     cerr<< "Input still has undefined materials." << endl;
  //     return result ;
  //   }
  // }
  if (s.meshes.size()==0){
    cerr<< "No meshes defined in the input." << endl;
    return result;
  }

  unsigned int i;

  std::list<node>::iterator node_it;
  std::list<triangle>::iterator tri_it;

  //pivot points
  vector<Point> pivots;

  //The vector that holds the solution, namely:
  vector<vector<Real> > pivslopes;

  // the intervals must be either all negative or all positive.
  // pivsign holds that information
  vector<bool> pivsign;

  // allowable and unallowable intervals
  // represented by a line  between two points
  vector<vector<vector<Point> > > allow;
  vector<vector<vector<Point> > > unall;

  // pivot compound intervals:
  vector<compoundInterval> pivint;


  Real ex[2],ep[2];

  i=0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++,i++){
    vector<Point> dum1;
    ex[0] = ex[1] = (*mesh_it->triangles.begin()->cent)(1);
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!= mesh_it->triangles.end(); tri_it++){
      (*tri_it->cent)(1) < ex[0] ? ex[0]=(*tri_it->cent)(1) : 0;
      (*tri_it->cent)(1) > ex[1] ? ex[1]=(*tri_it->cent)(1) : 0;
    }
    ep[0] = min(mesh_it->mat->epsult[0],mesh_it->mat->epsult[1]) ;
    ep[1] = max(mesh_it->mat->epsult[0],mesh_it->mat->epsult[1]) ;
    
    Point a(ex[0],ep[0]);
    Point b(ex[1],ep[0]);
    Point c(ex[0],ep[1]);
    Point d(ex[1],ep[1]);

    pivots.push_back(a);
    pivots.push_back(b);
    pivots.push_back(c);
    pivots.push_back(d);

    pivsign.push_back(1);
    pivsign.push_back(0);
    pivsign.push_back(0);
    pivsign.push_back(1);

    vector<vector<Point> > dum2;

    dum1.push_back(a);
    dum1.push_back(c);
    dum2.push_back(dum1);
    dum1.clear();

    dum1.push_back(b);
    dum1.push_back(d);
    dum2.push_back(dum1);
    dum1.clear();

    allow.push_back(dum2);
    dum2.clear();


    dum1.push_back(a);
    dum1.push_back(b);
    dum2.push_back(dum1);
    dum1.clear();

    dum1.push_back(c);
    dum1.push_back(d);
    dum2.push_back(dum1);
    dum1.clear();

    unall.push_back(dum2);
    dum2.clear();
    
  }

  // compoundInterval lol;
  // lol._intervals.push_back(interval(0,5));
  // lol  += interval(4,6);
  // cout << "ASDASDASD " << lol << endl;

  for (unsigned int i = 0; i < pivots.size(); i++){
    Real upper,lower;
    pivint.push_back(compoundInterval());
    
    vector<compoundInterval> suballowable;

    for (unsigned int j = 0; j < allow.size(); j++){
      suballowable.push_back(compoundInterval());

      for (unsigned int k = 0; k < allow[j].size(); k++){
        upper = numtk::slope(pivots[i], allow[j][k][0]);
        lower = numtk::slope(pivots[i], allow[j][k][1]);
        cout <<i<< ": A " << lower << " " << upper <<endl;

        if (numtk::isnan(upper) || numtk::isinf(upper) ||
            numtk::isnan(lower) || numtk::isinf(lower))
          continue;
//        cout << "s " << suballowable[j].size() << endl;
        suballowable[j]+=interval(upper,lower);
//        cout << "e " << suballowable[j].size() << endl<<endl;
      }
    }
    
    if (suballowable.size() != 0){
     compoundInterval lol = suballowable[0];
//      vector<interval> lol = suballowable[0];
      for (unsigned int j = 0; j < suballowable.size(); j++){
        lol = lol.intersection(suballowable[j]);
//        cout << "ASDASD " << suballowable[j].size() << endl;
      }
      pivint[i]=lol;
    }else{
      continue;
    }

  // subtraction of unallowable intervals
#if 1
    for (unsigned int j = 0; j < unall.size(); j++){
      for (unsigned int k = 0; k < unall[j].size(); k++){
        
        upper = numtk::slope(pivots[i], unall[j][k][0]);
        lower = numtk::slope(pivots[i], unall[j][k][1]);
        
        if (numtk::isnan(upper) && numtk::isnan(lower)) continue;
        if (numtk::isinf(upper) && numtk::isinf(lower)) continue;
        cout <<i<< ": U " << lower << " " << upper <<endl;
        
        if (numtk::isinf(lower)){
          bool dumb = numtk::isposit(unall[j][k][1](0)-pivots[i](0)) *
            !numtk::isposit(upper);
          pivint[i].subtract_infinity(upper,!dumb);
          continue;
        }
        if (numtk::isinf(upper)){
          bool dumb = numtk::isposit(unall[j][k][0](0)-pivots[i](0))*
            !numtk::isposit(lower);
          pivint[i].subtract_infinity(lower,!dumb);
          continue;
        }
        if (numtk::isnan(lower)){
          pivint[i].subtract_infinity(0,!pivsign[i]);
          continue;
        }
        if (numtk::isnan(upper)){
          pivint[i].subtract_infinity(0,!pivsign[i]);
          continue;
        }
        
        // if (numtk::sign(lower) != numtk::sign(upper)){
        //   bool dumb = numtk::isposit(upper);
        //   subtract_infinity(upper, dumb,pivint[i]);
        //   subtract_infinity(lower, !dumb,pivint[i]);
        //   continue;
        // }
        pivint[i]-=interval(upper,lower);
      }
    }
#endif
    cout << endl;

  }


//  cout << allow.size() << " " << unall.size() << endl;
//   interval dd; 
//   dd.e = 9.13946e-06; dd.s =-0; pivint[0].push_back(dd);
//   dd.e = 0.000308364; dd.s =-0; pivint[0].push_back(dd);
//   dd.e = 0.000271017; dd.s =-0; pivint[0].push_back(dd);
// //  dd.s = ; dd.e =; pivint[0].push_back(dd);
// //  subtract_interval(-1.5,-5,pivint[0]);
//   subtract_interval(0.000271017, 0.000308364,pivint[0]);


  for (unsigned int i = 0; i < pivots.size(); i++){
    cout << i << "("<< pivint[i].size() <<") ("
         <<pivots[i](0)<<","<<pivots[i](1)<< ") ** " << pivint[i];
    cout << endl;
  }

  unsigned int j = 0;
  for (unsigned int i=0; i < pivint.size(); i++)
    j+=pivint[i].size();
  if (j==0) {
    cerr<< "No solutions possible for this combination." << endl;
    return result;
  }


#if 1

  vector<Real> moments;
  vector<Real> forces;

  ofstream of2; 
  sprintf(dumchar, "%s_int.txt",ofilename);
  of2.open(dumchar);
  cout << "Writing " << dumchar << "." << endl;
  for (unsigned int i=0 ; i<pivots.size(); i++){
//    pivslopes.push_back(vector<Real>());
    for (unsigned int k=0 ; k<pivint[i].size(); k++){
      vector<Real> slopes = numtk::range(pivint[i][k].s, pivint[i][k].e,ndiv);

      // vector<Real> lulz = 
      //   root_bisect_force
      //   (s,pivots[i](1),pivots[i](0),pivint[i][k].s,pivint[i][k].e,0);
      // pivslopes[i].insert(pivslopes[i].begin(),lulz.begin(),lulz.end());
      // cout << lulz.size() << endl;
      
      sprintf(dumchar, "%s_b%d%d0",ofilename,i,k);
      if (plotothers) plot_y_vs_sig(s,pivots[i](1), pivots[i](0), pivint[i][k].s, dumchar);
      sprintf(dumchar, "%s_b%d%d1",ofilename,i,k);
      if (plotothers) plot_y_vs_sig(s,pivots[i](1), pivots[i](0), pivint[i][k].e, dumchar);
      
      sprintf(dumchar, "%s_lol%d%d0.vtk",ofilename,i,k);
      if (plotothers) stress_vtk(s, pivots[i](1), pivots[i](0), pivint[i][k].s,dumchar);
      sprintf(dumchar, "%s_lol%d%d1.vtk",ofilename,i,k);
      if (plotothers) stress_vtk(s, pivots[i](1), pivots[i](0), pivint[i][k].e,dumchar);
      // vector<Real> forces;
      // for (unsigned int j=0 ; j<slopes.size(); j++)
      //   forces.push_back( normal_force(s,pivots[i](1),pivots[i](0),slopes[j]) );
      // sprintf(dumchar, "a%02d%02d.txt",i,k);
      // numtk::output_gp(slopes,forces,dumchar);

      for (unsigned int j=0 ; j<slopes.size(); j++){
        moments.push_back( moment(s,pivots[i](1),pivots[i](0),slopes[j])/1e6 );
        forces.push_back( normal_force(s,pivots[i](1),pivots[i](0),slopes[j])/1e3 );
        of2 << moments.back() << " " << forces.back() << " " << slopes[j] << endl;
      }
      of2 << endl;
//      sprintf(dumchar, "a%02d%02d.txt",i,k);
//      numtk::output_gp(slopes,moments,dumchar);

    }
  }
//  numtk::output_gp(moments,forces,"interaction2.txt");
  of2.close();
#endif



#if 0
  ofstream of; of.open("interaction.txt");
  vector<Real> indep = numtk::range(ranst,ranen,ndiv);
  
  for (unsigned int i=0 ; i<pivint.size(); i++){
    for (unsigned int j=0 ; j<pivint[i].size(); j++){
      for (unsigned int k=0 ; k<indep.size(); k++){
#ifndef USE_FORCE
        vector<Real> root = root_bisect_moment
          (s, pivots[i](1), pivots[i](0), pivint[i][j].s, pivint[i][j].e, indep[k]);
#else
        vector<Real> root = root_bisect_force
          (s, pivots[i](1), pivots[i](0), pivint[i][j].s, pivint[i][j].e, indep[k]);
#endif
        cout << "ASDASDASD " << indep[k]<<" " <<root.size() << endl;
        for (unsigned int l=0; l<root.size(); l++){
          cout << "  " << l << ": " << indep[k] << " " << root[l] << endl;
#ifndef USE_FORCE
          Real for1 = normal_force(s,pivots[i](1),pivots[i](0),root[l]);
          of << indep[k]/1e6 << " " << for1/1e3 << endl;
#else
          Real mom1 = moment(s,pivots[i](1),pivots[i](0),root[l]);
          of << mom1/1e6 << " " << indep[k]/1e3 << endl;
#endif
        }
      }
    }
  }

  of.close();
  numtk::exec("gnuplot -p -e \"plot 'interaction.txt' using 1:2 w p \"");
#endif



#if 0
  cout << endl<<endl<<endl << "LOL" << endl;
  Real indep1 = 0;
  for (unsigned int i=0 ; i<pivint.size(); i++){
    for (unsigned int j=0 ; j<pivint[i].size(); j++){
      cout << i << "  " << j << endl; 
#ifndef USE_FORCE
      vector<Real> root = root_bisect_moment
        (s, pivots[i](1), pivots[i](0), pivint[i][j].s, pivint[i][j].e, indep1);
#else
      vector<Real> root = root_bisect_force
        (s, pivots[i](1), pivots[i](0), pivint[i][j].s, pivint[i][j].e, indep1);
#endif
      for (unsigned int l=0; l<root.size(); l++){
//          cout << "  " << l << ": " << indep1 << " " << root[l] << endl;
#ifndef USE_FORCE
        Real for1 = normal_force(s,pivots[i](1),pivots[i](0),root[l]);
        Real mom1 = moment(s,pivots[i](1),pivots[i](0),root[l]);
        cout << "  AD " << indep1 << " " << for1 <<  " " << mom1 << endl;
#else
        Real for1 = normal_force(s,pivots[i](1),pivots[i](0),root[l]);
        Real mom1 = moment(s,pivots[i](1),pivots[i](0),root[l]);
        cout << mom1 << " " << indep1 << endl;
        cout << "  AD " << indep1 << " " << for1 <<  " " << mom1 << endl;
#endif
        
      }
    }
  }
#endif

  
}



void plot_y_vs_sig( supermesh &s, Real a, Real b, Real slope, string lolchar ){

  list<mesh>::iterator mesh_it;
  std::list<node>::iterator node_it;
  std::list<triangle>::iterator tri_it;

  vector<Real> epsv,y,sigv;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){

    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++){
      y.push_back( (*tri_it->cent)(1) );
      epsv.push_back(  slope* ( y.back() - b) + a );
      sigv.push_back( mesh_it->mat->sig(epsv.back()) );
//      cout << y.back() << " " << epsv.back() << " " << sigv.back() << endl;
    }
  }
//  cout << endl;
  string c1;
  c1=lolchar+".txt";
  numtk::output_gp(y,epsv,sigv,c1);

  
  printf(" a = %.2e  b = %.2e  slope = %.2e\n",a,b,slope);



}





#if 0


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




#endif

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





// void union_interval(vector<interval> &v)
// {
//   if (v.size()==0) return;
//   sort(v.begin(), v.end(), interval_sorter);
  
//   vector<interval> r;
//   r.push_back(v[0]);
//   for (unsigned int i = 0; i < v.size(); i++){

//     if (r.back().e < v[i].s)
//       r.push_back(v[i]);
//     else if (r.back().e == v[i].s)
//       r.back().e = v[i].e;
//     if (r.back().e < v[i].e)
//       r.back().e = v[i].e;
//   // if y[-1].end < x.start:
//   //     y.append(x)
//   // elif y[-1].end == x.start:
//   //     y[-1].end = x.end
//   }
//   v=r;

// }

// void add_interval(Real a, Real b, vector<interval> &v)
// {
//   if (numtk::iszero(a-b)) return;

//   interval c;
//   c.s = min(a,b);
//   c.e = max(a,b);
//   if (v.size()==0){
//     v.push_back(c); return;
//   }
//   v.push_back(c);
//   union_interval(v);
//   // cout << c.s << " " << c.e << endl;
//   // for (unsigned int j = 0; j < v.size(); j++)
//   //   cout << "(" << v[j].s << "," <<  v[j].e << "), ";
//   // cout << endl<<endl;
  
// }


// void subtract_interval(Real a, Real b, vector<interval> &v)
// {
//   if (numtk::iszero(a-b)) return;
//   interval c,d;
//   c.s = min(a,b);
//   c.e = max(a,b);

//   sort(v.begin(), v.end(), interval_sorter);

  
//   vector<interval> r;
//   for (unsigned int i = 0; i < v.size(); i++){
//     if (c.e <= v[i].s){
//       r.push_back(v[i]);
//     }else if (c.s <= v[i].s && c.e >= v[i].s && c.e <= v[i].e){
//       d.s = c.e; d.e = v[i].e; r.push_back(d);
//     } else if (c.s <= v[i].s && c.e >= v[i].e){
//       continue;
//     } else if (c.s >= v[i].s && c.e <= v[i].e){
//       d.s=v[i].s; d.e=c.s; r.push_back(d);
//       d.s=c.e; d.e=v[i].e; r.push_back(d);
//     } else if (c.s >= v[i].s && c.s <= v[i].e && c.e >= v[i].e){
//       d.s = v[i].s; d.e = c.s; r.push_back(d);
//     } else if (c.s >= v[i].e){
//       r.push_back(v[i]);
//     }
//   }
//   v.clear();
//   for (unsigned int i = 0; i < r.size(); i++)
//     if (r[i].s!=r[i].e && ~numtk::iszero(r[i].s - r[i].e))
//         v.push_back(r[i]);

// // cout << c.s << " " << c.e << endl;
//   // for (unsigned int j = 0; j < v.size(); j++)
//   //   cout << "(" << v[j].s << "," <<  v[j].e << "), ";
//   // cout << endl<<endl;

// // The algorithm in python

// // for x in s[:]:
// //     if z.end < x.start:
// //         break
// //     elif z.start < x.start and z.end > x.start and z.end < x.end:
// //         x.start=z.end
// //     elif z.start < x.start and z.end > x.end:
// //         s.remove(x)
// //     elif z.start > x.start and z.end < x.end:
// //         s.append(tp(x.start,z.start))
// //         s.append(tp(z.end,x.end))
// //         s.remove(x)
// //     elif z.start > x.start and z.start < x.end and z.end > x.end:
// //         x.end=z.start
// //     elif z.start > x.end:
// //         continue
// }

// vector<interval> intersection_interval(vector<interval> a, vector<interval> b)
// {
//   union_interval(a); union_interval(b);
//   vector<interval> r;
//   interval d;
//   for (unsigned int i = 0; i < a.size(); i++){
//     for (unsigned int j = 0; j < b.size(); j++){
//       if (a[i].e < b[i].s || a[i].s > b[i].e){
//         continue;
//       }else if ( a[i].s <= b[i].s && a[i].e >= b[i].s && a[i].e <= b[i].e){
//         d.s = b[i].s; d.e = a[i].e; r.push_back(d);
//       }else if ( a[i].s <= b[i].s && a[i].e >= b[i].e ){
//         r.push_back(b[i]);
//       }else if ( a[i].s >= b[i].s && a[i].e <= b[i].e ){
//         r.push_back(a[i]);
//       }else if ( a[i].s >= b[i].s && a[i].s <= b[i].e && a[i].e >= b[i].e){
//         d.s = a[i].s; d.e = b[i].e; r.push_back(d);
//       }
//     }
//   }
//   vector<interval> s;
//   for (unsigned int i = 0; i < r.size(); i++)
//     if (r[i].s!=r[i].e && ~numtk::iszero(r[i].s - r[i].e))
//         s.push_back(r[i]);

//   return s;

// }


// void subtract_infinity(Real a, bool side, vector<interval> &v)
// {
//   interval d;

//   // side == 1 for positive, == 0 for negative
  
//   vector<interval> r;
//   if (side == true){
//     for (unsigned int i = 0; i < v.size(); i++){
//       if (a <= v[i].s && a <= v[i].e){
//         continue;
//       }else if (a >= v[i].s && a <= v[i].e){
//         d.s = v[i].s; d.e = a; r.push_back(d);
//       } else if (a >= v[i].s && a >= v[i].e){
//         r.push_back(v[i]);
//       }
//     }
//   }else {
//     for (unsigned int i = 0; i < v.size(); i++){
//       if (a <= v[i].s && a <= v[i].e){
//         r.push_back(v[i]);
//       }else if (a >= v[i].s && a <= v[i].e){
//         d.s = a; d.e = v[i].e; r.push_back(d);
//       } else if (a >= v[i].s && a >= v[i].e){
//         continue;
//       }
//     }
//   }
//   v.clear();
//   for (unsigned int i = 0; i < r.size(); i++)
//     if (r[i].s!=r[i].e && ~numtk::iszero(r[i].s - r[i].e))
//         v.push_back(r[i]);


// }
