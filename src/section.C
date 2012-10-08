#include "section.h"
#include "interval.h"

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
      F+=tri_it->area*sig;
    }
  }
//  printf(" a= %.2e   b= %.2e   s= %.2e   F= %.2e\n",a,b,slope,F);
  return F;
}


void normal_and_moment( supermesh &s, Real F, Real M, Real a, Real b, Real slope){
  list<mesh>::iterator mesh_it;
  list<triangle>::iterator tri_it;
  F = 0, M = 0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
    Real eps, sigma, y, dist;
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++){
      y = (*tri_it->cent)(1);
      eps = slope * (y - b) + a ;
      sigma = mesh_it->mat->sig(eps) ;
      dist = y - b;
      F += tri_it->area * sigma;
      M += tri_it->area * sigma * dist;
    }
  }
  printf("a= %.2e   b= %.2e   s= %.2e   F= %.2e   M= %.2e\n",a,b,slope,F,M);
}


Real moment( supermesh &s, Real a, Real b, Real slope){
  list<mesh>::iterator mesh_it;
  list<triangle>::iterator tri_it;
  Real M = 0;
  for (mesh_it = s.meshes.begin(); mesh_it != s.meshes.end(); mesh_it++){
    Real eps, sigma, y, dist;
    for (tri_it = mesh_it->triangles.begin(); 
         tri_it!=mesh_it->triangles.end(); tri_it++){
      y = (*tri_it->cent)(1);
      eps = slope * (y - b) + a ;
      sigma = mesh_it->mat->sig(eps) ;
      dist = y - b;
      M += tri_it->area * sigma * dist;
    }
  }
  printf("a= %.2e   b= %.2e   s= %.2e   M= %.2e\n",a,b,slope,M);
  return M;
}





namespace numtk{
  inline Real slope(const Point &a, const Point &b){
    return (a(1) - b(1)) / (a(0) - b(0));
  }
}

void interaction(supermesh &s)
{
  list<mesh>::iterator mesh_it;
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
  vector<Point> pivots;
  // allowable and unallowable intervals
  // represented by a line  between two points
  vector<vector<vector<Point> > > allow;
  vector<vector<vector<Point> > > unall;

  // pivot compound intervals:
  vector<compoundInterval> pivint;
  // the intervals must be either all negative or all positive.
  // pivsign holds that information
  vector<bool> pivsign;


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

  compoundInterval lol;
  lol._intervals.push_back(interval(0,5));
  lol  += interval(4,6);
  cout << "ASDASDASD " << lol << endl;




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

  // subtraction of unallowed intervals
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
         <<pivots[i](0)<<","<<pivots[i](1)<< ") ** " ;
    for (unsigned int j = 0; j < pivint[i].size(); j++)
      cout << "(" << pivint[i][j].s << "," <<  pivint[i][j].e << "), ";
    cout << endl;
  }

  unsigned int j = 0;
  for (unsigned int i=0; i < pivint.size(); i++)
    j+=pivint[i].size();
  if (j==0) {
    cerr<< "No solutions possible for this combination." << endl;
    return;
  }


  char cucu[1024],cucu2[1024];
  for (unsigned int lol=0;lol<pivots.size(); lol++){
    sprintf(cucu2, "b%d",lol);
//    plot_y_vs_sig(s,pivots[lol](1),pivots[lol](0),pivint[lol][0].s,cucu2);
  }
//  plot_y_vs_sig(s,pivots[1](1),pivots[1](0),-9.09902e-06,"cc");
  for (unsigned int i=0 ; i<pivots.size(); i++){
    for (unsigned int k=0 ; k<pivint[i].size(); k++){
      vector<Real> slopes = numtk::range(pivint[i][k].s, pivint[i][k].e,100);
      vector<Real> forces;

      sprintf(cucu2, "b%d%d0",i,k);
      plot_y_vs_sig(s,pivots[i](1), pivots[i](0), pivint[i][k].s, cucu2);
      sprintf(cucu2, "b%d%d1",i,k);
      plot_y_vs_sig(s,pivots[i](1), pivots[i](0), pivint[i][k].e, cucu2);
      
      for (unsigned int j=0 ; j<slopes.size(); j++)
        forces.push_back( normal_force(s,pivots[i](1),pivots[i](0),slopes[j]) );
      sprintf(cucu, "a%02d%02d.txt",i,k);
//      numtk::output_gp(slopes,forces,cucu);
      numtk::output_gp(slopes,forces,cucu);

 
    }

  }
  // for (Real r = -0.0035; r < 0 ; r+=0.0001)
  //   cout << r<<" " <<s.meshes.front().mat->sig(r) <<endl;

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
