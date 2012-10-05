
#include "pr1.h"
#include "prgl.h"
//#include "predicates.c"
#include "predicates.h"
#include "triangle.h"
#include <time.h>

//#define INSERTINTOVECTOR(a,b) a.insert(a.begin(),b.begin(),b.end())
namespace numtk{
  std::vector<Real> range(Real st, Real en, unsigned int a){
    std::vector<Real> result;
    for (unsigned int i=0; i < a; i++){
      result.push_back(st+i*(en-st)/a);
    }
    result.push_back(en);
    return result;
  }
  void output_gp(std::vector<Real> &a, std::vector<Real> &b, char *filename){
    cout << "Writing " << filename << "." << endl;
    std::ofstream afile;
    afile.open(filename);
    for (unsigned int i=0; i<a.size();i++)
      afile << a[i] << " " << b[i] << std::endl;
    afile.close();
  }
  Real distance(node a, Real m, Real n){
    return (Real) fabs(-1*m*a(0)+a(1)-n)/sqrt(m*m+1);
  }
}

unsigned int random_range(unsigned int min, unsigned int max) //range : [min, max)
{
  static bool first = true;
  if ( first ) 
  {  
    srand(time(NULL)); //seeding for the first time only!
    first = false;
  }
//  std::cout<< (max)<<std::endl;
  if (max-min != 0){
    return min + rand() % (max - min);
  }else{
    std::cerr<<"Invalid range for random number generation"<<std::endl;
    return 0;
  }
}
Point random_vector(Real rad){
// t = 2*pi*random()
// u = random()+random()
// r = if u>1 then 2-u else u
// [r*cos(t), r*sin(t)]
  static bool first = true;
  if ( first ) 
  {  
    srand(time(NULL)); //seeding for the first time only!
    first = false;
  }
//  std::cout << "ASDAS " <<   random_Real()<< std::endl;
  Real t = 2 * M_PI * random_Real();
  Real u = random_Real() + random_Real();
  Real r = u>1 ? 2-u : u;
  
  return Point(rad*r*cos(t),rad*r*sin(t), 0);

}

Real calculate_triangle_area_2d(Real a[2], Real b[2], Real c[2]){
  Real d[2],e[2], area;
  for (unsigned int i = 0; i<2; i++) {
    d[i] = c[i] - a[i];
    e[i] = b[i] - a[i];
  }
  area = fabs( (d[0]*e[1] - e[0]*d[1])/2 );
  return area;
}

Real calculate_triangle_area_2d(node &aa, node &bb, node &cc){
  Real a[2],b[2],c[2],d[2],e[2], area;
  for (unsigned int i = 0; i<2; i++) {
    a[i]=aa(i); b[i]=bb(i);  c[i]=cc(i);
  }
  for (unsigned int i = 0; i<2; i++) {
    d[i] = c[i] - a[i];
    e[i] = b[i] - a[i];
  }
  area = fabs( (d[0]*e[1] - e[0]*d[1])/2 );
//  std::cout<< "Area " << area<< std::endl;
  return area;
}
bool triangle::is_node_inside(node &n){
  Real a0 = calculate_triangle_area_2d(n,*nodes[0],*nodes[1]);
  Real a1 = calculate_triangle_area_2d(n,*nodes[1],*nodes[2]);
  Real a2 = calculate_triangle_area_2d(n,*nodes[2],*nodes[0]);
  // Don't recalculate triangle's original area.
  if ( fabs(area - (a0 + a1 + a2)) > VOLUME_TOLERANCE ) {
    return false;
  } else {
    return true;
  }
}

bool mesh::is_node_inside(node &n){
  std::list<triangle>::iterator it;
  for (it = triangles.begin(); it != triangles.end(); it++){
    if (it->is_node_inside(n)){
#ifdef DEBUG
      std::cout << "Node is inside."<< std::endl;
#endif
      return true;
    }
  }
  return false;
}


triangle::triangle(node *n0, node *n1, node *n2, mesh &m){
  line *l0, *l1, *l2;
  m.lines.push_front( line(n0,n1 , true) );
  l0=&(*m.lines.begin());
  m.lines.push_front( line(n1,n2 , true) );
  l1=&(*m.lines.begin());
  m.lines.push_front( line(n2,n0 , true) );
  l2=&(*m.lines.begin());
  this->lines[0]=l0;
  this->lines[1]=l1;
  this->lines[2]=l2;
  this->nodes[0]=n0;
  this->nodes[1]=n1;
  this->nodes[2]=n2;
  calc_area();
}


triangle::triangle(node *n0, node *n1, node *n2){
  line *l0, *l1, *l2;
//  m.lines.push_front( line(n0,n1 , true) );
//  l0=&(*m.lines.begin());
//  m.lines.push_front( line(n1,n2 , true) );
//  l1=&(*m.lines.begin());
//  m.lines.push_front( line(n2,n0 , true) );
//  l2=&(*m.lines.begin());
  // this->lines[0]=l0;
  // this->lines[1]=l1;
  // this->lines[2]=l2;
  this->nodes[0]=n0;
  this->nodes[1]=n1;
  this->nodes[2]=n2;
  calc_area();
}



triangle::triangle(line *l1, line *l2, line *l3){
  this->lines[0]=l1;
  this->lines[1]=l2;
  this->lines[2]=l3;
  std::set<node*> n;
  n.insert(l1->nodes[0]); n.insert(l1->nodes[1]);
  n.insert(l2->nodes[0]); n.insert(l2->nodes[1]);
  n.insert(l3->nodes[0]); n.insert(l3->nodes[1]);
  unsigned int i=0;
  for(  std::set<node*>::iterator it = n.begin(); i<3 ; i++,it++ ){ 
    this->nodes[i]=*it;
  }
  calc_area();
//  std::cout << n.size() << std::endl;
}

Real triangle::calc_area(){
  area = calculate_triangle_area_2d( *(nodes[0]) ,*(nodes[1]),*(nodes[2]) );
  return area;
}


bool do_lines_intersect_exclusive(REAL a[2], REAL b[2], REAL c[2], REAL d[2]){
// Exclusive of the endpoints of the lines
// To spell it out: suppose you're looking at two line segments, 
// [AB] and [CD]. The segments intersect if and only if ((A and 
// B are of different sides of [CD]) and (C and D are on different 
// sides of [AB])).

// To see whether two points, P and Q, are on different sides of a 
// line segment [EF], compute two cross products, one for P and one for Q:

// (Fx - Ex)(Py - Fy) - (Fy - Ey)(Px - Fx)

// (Fx - Ex)(Qy - Fy) - (Fy - Ey)(Qx - Fx)

// If the results have the same sign (both positive or both negative) 
// then forget it, the points are on the same side, the segments do not 
// intersect. If one is positive and the other negative, then the points 
// are on opposite sides. 

  static bool first2 = true;
  if (first2) 
  {  
    predicates::exactinit();
    first2 = false;
    return false;
  }

  bool p1, p2, p3, p4;
  REAL r1, r2, r3, r4;
//  p1 = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]) > 0 ? 1 : 0;
//  p2 = (b[0]-a[0])*(d[1]-b[1]) - (b[1]-a[1])*(d[0]-b[0]) > 0 ? 1 : 0;
//  std::cout<<orient2d(a,b,c)<< " " << orient2d(a,b,d) << std::endl;

  r1 = predicates::orient2dexact(a,b,c);
  r2 = predicates::orient2dexact(a,b,d);

  p1 = r1 > 0 ? 1 : 0;
  p2 = r2 > 0 ? 1 : 0;

  if (p1==p2) return false;

//  p3 = (d[0]-c[0])*(a[1]-d[1]) - (d[1]-c[1])*(a[0]-d[0]) > 0 ? 1 : 0;
//  p4 = (d[0]-c[0])*(b[1]-d[1]) - (d[1]-c[1])*(b[0]-d[0]) > 0 ? 1 : 0;
  r3 = predicates::orient2dexact(c,d,a);
  r4 = predicates::orient2dexact(c,d,b);

  p3 = r3 > 0 ? 1 : 0;
  p4 = r4 > 0 ? 1 : 0;

  if (fabs(r1)<TOLERANCE || fabs(r2)<TOLERANCE || 
      fabs(r3)<TOLERANCE || fabs(r4)<TOLERANCE){
    return false;
  }
  return p3!=p4;

}


bool do_lines_intersect_inclusive(REAL a[2], REAL b[2], REAL c[2], REAL d[2]){
// inclusive of the end points of the lines
  static bool first2 = true;
  if (first2) {  
    predicates::exactinit();
    first2 = false;
    return false;
  }

  bool p1, p2, p3, p4;
  REAL r1, r2, r3, r4;

  r1 = predicates::orient2dexact(a,b,c);
  r2 = predicates::orient2dexact(a,b,d);
  p1 = r1 > 0 ? 1 : 0;
  p2 = r2 > 0 ? 1 : 0;
  if (p1==p2) return false;
  r3 = predicates::orient2dexact(c,d,a);
  r4 = predicates::orient2dexact(c,d,b);
  p3 = r3 > 0 ? 1 : 0;
  p4 = r4 > 0 ? 1 : 0;
  return p3!=p4;
}

bool line::intersect_line(line &l) {
  Real a[2],b[2],c[2],d[2];
  for (unsigned int i = 0 ; i<2 ; i++){
    a[i] = (*this->nodes[0])(i);
    b[i] = (*this->nodes[1])(i);
    c[i] = (*l.nodes[0])(i);
    d[i] = (*l.nodes[1])(i);
  }
  return do_lines_intersect_exclusive(a,b,c,d);
}


line::line(node *n1, node *n2, const bool rh) {
  this->nodes[0]=n1;
  this->nodes[1]=n2;
  this->rh = rh;
}




// line::line(node ns[2]) {
//   this->nodes[0]=&ns[0];
//   this->nodes[1]=&ns[1];
//   rh=true;
// }
bool line::operator==(const line &l){
//  std::cout << "ASDASDASDASDASDASDASD";
  if( (this->nodes[0] == this->nodes[0] && this->nodes[1] == this->nodes[1]) ||
      (this->nodes[0] == this->nodes[1] && this->nodes[1] == this->nodes[0]) ){
    return true;
  } else {
    return false;
  }
}

void rectangular_section(mesh &m, Real left, Real right, Real top, Real bottom){
  m.clear();
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;

  node *n[4];
  for (unsigned int i=0; i< 4; i++) m.nodes.push_back(node());
  unsigned int i = 0;
  for (node_it=m.nodes.begin();node_it!=m.nodes.end() ; node_it++,i++) {
    n[i] = &(*node_it);
  }

 (*n[0])(0)=left;  (*n[0])(1)=bottom;
  (*n[1])(0)=right;  (*n[1])(1)=bottom;
  (*n[2])(0)=right;  (*n[2])(1)=top;
  (*n[3])(0)=left;  (*n[3])(1)=top;

  m.lines.push_back(line(n[0],n[1]));
  m.lines.push_back(line(n[1],n[2]));
  m.lines.push_back(line(n[2],n[3]));
  m.lines.push_back(line(n[3],n[0]));


  for (line_it=m.lines.begin(); line_it!=m.lines.end(); line_it++) 
    m.hull.push_back( &(*line_it) );
}


void circular_section(mesh &m, Real rad, Real cen[2], int nsides){
  m.clear();
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
  m.innode.push_back( node(cen[0],cen[1]) );

//  for (unsigned int i=0; i< 4; i++) m.nodes.push_back(node());
  unsigned int i = 0;
  Real th,x,y,dum;
//  int nsides=10;
  for(int j=0; j<nsides  ; j++)
  {
//    if (360%==j) continue;
    dum=360*j/nsides;
    th=M_PI * dum / 180.0;
    x = rad * cos(th)+cen[0];
    y = rad * sin(th)+cen[1];
    m.nodes.push_back(node(x,y));
  }

  std::vector<node*> node_ptrs;
  for (node_it = m.nodes.begin(); node_it!=m.nodes.end(); node_it++)
    node_ptrs.push_back(&(*node_it));

  for (unsigned int i=0; i<node_ptrs.size()-1; i++){
    m.lines.push_back( line(node_ptrs[i], node_ptrs[i+1]) );
  }
  m.lines.push_back( line(node_ptrs.back(), node_ptrs.front()) );
  
// (*n[0])(0)=0;  (*n[0])(1)=0;
//  (*n[1])(0)=width;  (*n[1])(1)=0;
//  (*n[2])(0)=width;  (*n[2])(1)=height;
//  (*n[3])(0)=0;  (*n[3])(1)=height;

  // m.lines.push_back(line(n[0],n[1]));
  // m.lines.push_back(line(n[1],n[2]));
  // m.lines.push_back(line(n[2],n[3]));
  // m.lines.push_back(line(n[3],n[0]));


  for (line_it=m.lines.begin(); line_it!=m.lines.end(); line_it++) 
    m.hull.push_back( &(*line_it) );
}


void mesh::init_centroids(){
  std::list<triangle>::iterator tri_it;
  centroids.clear();

  for (tri_it = triangles.begin(); tri_it != triangles.end(); tri_it++){
    centroids.push_back(node(  
                          ((*tri_it->nodes[0])(0)+
                           (*tri_it->nodes[1])(0)+
                           (*tri_it->nodes[2])(0))/3,
                          ((*tri_it->nodes[0])(1)+
                           (*tri_it->nodes[1])(1)+
                           (*tri_it->nodes[2])(1))/3
                            ));
    tri_it->cent = &centroids.back();
  }

}

void mesh::init_node2tri(){
  std::list<triangle>::iterator tri_it;
  for (tri_it = triangles.begin(); tri_it != triangles.end(); tri_it++){
    for (unsigned int i = 0; i<3; i++){
      node2tri[ &(*tri_it->nodes[0]) ].insert( &(*tri_it) );
      node2tri[ &(*tri_it->nodes[1]) ].insert( &(*tri_it) );
      node2tri[ &(*tri_it->nodes[2]) ].insert( &(*tri_it) );
    }
  }
}


void mesh::generate_seed_triangle(Real avlen, Real probrad){
// Seed triangle
  clear();
  node *n0,*n1,*n2;
  line *l0,*l1,*l2;
  for (unsigned int i=0; i< 3; i++) nodes.push_back(node());

  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
//  std::list<line*>::iterator bine_it;

  node_it=nodes.begin(); 
  n0 = &(*node_it);
  node_it++;
  n1 = &(*node_it);
  //nodes[0](0)=0;nodes[0](1)=0;nodes[0](2)=0;
  (*node_it)(0) = avlen; 
  *node_it += random_vector(probrad);

  *node_it++;
  n2 = &(*node_it);
  (*node_it)(0) = avlen/2;
  (*node_it)(1) = avlen*sqrt(3)/2;
  *node_it += random_vector(probrad);

  triangles.push_front(triangle(n0,n1,n2,*this));

  for (line_it=lines.begin(); line_it!=lines.end(); line_it++) 
    hull.push_back( &(*line_it) );
}

void report(struct triangulateio *io, int markers, int reporttriangles,
            int reportneighbors, int reportsegments, int reportedges,
            int reportnorms)
{
  int i, j;

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }
}


void mesh::triangulate_mesh(Real minarea)
{
  char lolchar[1024];

  struct triangulateio in, out, vorout;
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
  unsigned int i;

  numbernodes();
  /* Define input points. */
  if (nodes.size() == 0){
    std::cerr << "Lol no points." <<std::endl;
    return;
  }
  in.numberofpoints = nodes.size();
  in.numberofsegments  = lines.size();
  in.numberofpointattributes = 0;
  in.numberofregions = 0;

  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
  in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  
  if ( holes.size()>0){
    in.numberofholes = holes.size();
    in.holelist= (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    i=0;
    for (node_it = holes.begin(); node_it != holes.end(); node_it++,i++){
//      std::cout << holes.size() << " " <<(*node_it)(0) << std::endl;
      in.holelist[i*2]  = (*node_it)(0);
      in.holelist[i*2+1]= (*node_it)(1);
    }
  } else{
    in.numberofholes = 0;
    in.holelist= (REAL *) NULL;
  }

  // for (unsigned int i = 0; i< holes.size()*2; i++)
  //   std::cout << in.holelist[i] << std::endl;


  i=0;
  for (line_it=lines.begin(); 
       line_it!=lines.end(); line_it++,i++){
//    std::cout << "ASDASD" << (*line_it).nodes[0]->id <<std::endl;
    in.segmentlist[i*2]  = (*line_it).nodes[0]->id;
    in.segmentlist[i*2+1]= (*line_it).nodes[1]->id;
  }


  i=0;
  for (node_it=nodes.begin(); 
       node_it!=nodes.end(); node_it++,i++){
    in.pointlist[i*2]  = (*node_it)(0);
    in.pointlist[i*2+1]= (*node_it)(1);
  }

  in.pointmarkerlist = (int *) NULL;
  in.regionlist= (REAL *) NULL;
  in.segmentmarkerlist = (int *) NULL;
 
  // in.pointattributelist = (REAL *) malloc(in.numberofpoints *
  //                                         in.numberofpointattributes *
  //                                         sizeof(REAL));

  // in.numberofregions = 1;
  // in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
  // in.regionlist[0] = 0.5;
  // in.regionlist[1] = 5.0;
  // in.regionlist[2] = 7.0;            /* Regional attribute (for whole mesh). */
  // in.regionlist[3] = 0.1;          /* Area constraint that will not be used. */


//  report(&in, 0, 0, 0, 1, 0, 0);

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out' and a voronoi diagram in `vorout'.  */

  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  out.pointattributelist = (REAL *) NULL;
  out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  out.triangleattributelist = (REAL *) NULL;
//  out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  /* Needed only if segments are output (-p or -c) and -P not used: */
  out.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  out.segmentmarkerlist = (int *) NULL;
  out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), preserve the convex hull (c), number everything from  */
  /*   zero (z), assign a regional attribute to each element (A), and  */
  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  /*   neighbor list (n).                                              */

  sprintf(lolchar, "Qpqzea%.4f",(Real) minarea);
//  std::cout << lolchar << std::endl;
  triangulate(lolchar, &in, &out, (struct triangulateio *) NULL);
//  triangulate("pqzena5", &in, &out, NULL);
//  triangulate(triswitches, &in, &out, &vorout);

//  printf("Triangulation:\n\n");
//  report(&out, 0, 1, 1, 1, 1, 0);


  clear();
  for (i=0; i<out.numberofpoints; i++)
    nodes.push_back(node( out.pointlist[i*2],out.pointlist[i*2+1]));
  std::vector<node*> node_ptrs;
  for (node_it = nodes.begin(); node_it!=nodes.end(); node_it++)
    node_ptrs.push_back(&(*node_it));

  for (i=0; i<out.numberofedges; i++)
    lines.push_back(line( 
                      node_ptrs[out.edgelist[i*2]],
                      node_ptrs[out.edgelist[i*2+1]]));
  
  for (i=0; i<out.numberoftriangles; i++)
    triangles.push_back(triangle( 
                          node_ptrs[out.trianglelist[i*3]],
                          node_ptrs[out.trianglelist[i*3+1]],
                          node_ptrs[out.trianglelist[i*3+2]]));

  // Calculate the centroids of the elements.
  init_centroids();

//  printf("Initial Voronoi diagram:\n\n");
//  report(&vorout, 0, 0, 0, 0, 1, 1);

  if (in.holelist!=NULL) free(in.holelist);
  free(in.pointlist);
  // free(in.pointattributelist);
  // free(in.pointmarkerlist);
  // free(in.regionlist);
  free(out.pointlist);
  free(out.pointattributelist);
  free(out.pointmarkerlist);
  free(out.trianglelist);
  free(out.triangleattributelist);
//  free(out.trianglearealist);
//  free(out.neighborlist);
  free(out.segmentlist);
  free(out.segmentmarkerlist);
  free(out.edgelist);
  free(out.edgemarkerlist);


}


void mesh::subtract(mesh &m){
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;

  m.numbernodes();
  unsigned int lastindex=nodes.size();
//  node_it = nodes.end(); node_it--;
//  std::cout << *node_it;
  holes.insert(holes.begin(),m.innode.begin(), m.innode.end());

  nodes.insert(nodes.end(),m.nodes.begin(),m.nodes.end());

  std::vector<node*> node_ptrs;
  for (node_it = nodes.begin(); node_it!=nodes.end(); node_it++)
    node_ptrs.push_back(&(*node_it));

  for (line_it = m.lines.begin(); line_it!=m.lines.end(); line_it++){
    lines.push_back(line( node_ptrs[lastindex+line_it->nodes[0]->id],
                          node_ptrs[lastindex+line_it->nodes[1]->id]));
  }


}


void mesh::add(mesh &m){
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;

  m.numbernodes();
  unsigned int lastindex=nodes.size();

  nodes.insert(nodes.end(),m.nodes.begin(),m.nodes.end());
  innode.insert(innode.begin(),m.innode.begin(), m.innode.end());

  std::vector<node*> node_ptrs;
  for (node_it = nodes.begin(); node_it!=nodes.end(); node_it++)
    node_ptrs.push_back(&(*node_it));

  for (line_it = m.lines.begin(); line_it!=m.lines.end(); line_it++){
    lines.push_back(line( node_ptrs[lastindex+line_it->nodes[0]->id],
                          node_ptrs[lastindex+line_it->nodes[1]->id]));
  }

}


void mesh::numbernodes(){
  unsigned int i=0;
  for (std::list<node>::iterator it=nodes.begin(); it!=nodes.end() ; it++,i++){
    it->id=i;
  }
  i=0;
  for (std::list<line>::iterator it=lines.begin(); it!=lines.end() ; it++,i++){
    it->id=i;
  }
  i=0;
  for (std::list<triangle>::iterator it=triangles.begin(); it!=triangles.end() ; it++,i++){
    it->id=i;
  }
}


mesh::mesh(){
  nnodes=10;
  avlen=1;
  probrad=0.4;
  mat = NULL;
}


void mesh::generate1(unsigned int n, Real a, Real p){
  nnodes=n;
  avlen=a;
  probrad=p;
  generate1();
}


void mesh::generate1(){
  std::cout << "Generating...";
  clock_t start, end;
  start = clock();//time(NULL);

// Form the seed triangle
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
  std::list<line*>::iterator bine_it;

  generate_seed_triangle(avlen,probrad);

//  hull.insert(hull.begin(),lines.begin(),lines.end());
//  std::vector<line>::iterator rl;
#if 1
  for (unsigned int i = 0 ; (nodes.size() < nnodes)&&(i<10*nnodes) ; i++)
#else
  for (unsigned int i = 0 ; (nodes.size() < nnodes) ; i++)
#endif

  {
//    line rl = get_random_hull_elem();
    unsigned int rn=random_range(0,hull.size());
//    random_n< std::vector<line>::iterator > (hull.begin(),hull.end());
    add_node2(rn, probrad, avlen);
//    display();    glutPostRedisplay();
  }

  numbernodes();
  end = clock();//time(NULL);
  std::cout << " in " << (end-start) / (REAL) CLOCKS_PER_SEC << " seconds." << std::endl;

  print_info();

}







void mesh::print_info(){
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
  std::list<line*>::iterator bine_it;

  std::cout<< "# nodes      = " << nodes.size()<<std::endl;
  std::cout<< "# lines      = " << lines.size()<<std::endl;
  std::cout<< "# hull lines = " << hull.size()<<std::endl;
  std::cout<< "# triangles  = " << triangles.size()<<std::endl;

#if 0
  for (node_it=nodes.begin() ; node_it != nodes.end() ; node_it++){
    std::cout << &(*node_it) <<" " << *node_it <<std::endl;
  }
  std::cout << std::endl;
  for (line_it=lines.begin() ; line_it != lines.end() ; line_it++)
    std::cout << &(*line_it)<< " " << line_it->nodes[0]->id << " " 
              << line_it->nodes[1]->id << std::endl;

  for (bine_it=hull.begin() ; bine_it != hull.end() ; bine_it++){
    std::cout << (*bine_it)->id<< " " << (*bine_it) << std::endl;
  }
//    std::cout <<lines[i].nodes[0]->id <<" " << *(lines[i].nodes[0]) <<std::endl;
//    std::cout <<lines[i].nodes[0]->id <<" " << *(lines[i].nodes[1]) <<std::endl;

#endif

}

void mesh::color_nodes(unsigned int c){
  std::list<node>::iterator it;

  for (it= nodes.begin(); it!=nodes.end(); it++) it->color=c;
  
}
void mesh::color_lines(unsigned int c){
  std::list<line>::iterator it;
  for (it= lines.begin(); it!=lines.end(); it++) it->color=c;

}
void mesh::color_triangles(unsigned int c){
  std::list<triangle>::iterator it;
  for (it= triangles.begin(); it!=triangles.end(); it++) it->color=c;

}
void mesh::color_mesh(unsigned int a,unsigned int b,unsigned int c){
  color_nodes(a);
  color_lines(b);
  color_triangles(c);

}


void mesh::clear(){
  nodes.clear(); lines.clear(); hull.clear();triangles.clear();
  holes.clear(); centroids.clear();
}

bool mesh::add_node2(unsigned int ln, const Real probrad, const Real avlen){
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
  std::list<line*>::iterator bine_it;

  bine_it=hull.begin();
  std::advance(bine_it,ln);
  line *l=(*bine_it);//hull[ln];

  node *n0 = l->nodes[0];
  node *n1 = l->nodes[1];
//  Real length = calculate_line_length_2d(*n0,*n1);
  line *l0, *l1;
  Point median, difference, poot, final,finam;
  median = (*n1 + *n0);//(l.nodes[1] + l.nodes[0]);
  for (unsigned int i=0; i<2;i++) median(i)/=2;
//  std::cout<<median;

  difference = *n1 - *n0;
  poot(0)=-1*difference(1); poot(1)= difference(0);
  for (unsigned int i=0; i<2;i++) poot(i)*=sqrt(3)/2;

  Point trolol = random_vector(probrad);
  if (l->rh==true){
    final = median - poot + trolol;
  }else{
    final = median + poot + trolol;
  }

   // final(0)+=0.05;
   // finam(0)+=0.05;

//  std::cout<<final << std::endl;
    

//  node newnode(final(0),final(1));
  node lolnode(final(0),final(1));

  if ( is_node_inside(lolnode) ){
      return false;
  }

// Intersection test    
  line lolline1(&lolnode,n1,true);
  line lolline2(n0,&lolnode,true);

//  line lolline3(n1,&lolnode,true);line lolline4(&lolnode,n0,true);


#if 1
  for (std::list<line*>::iterator it = hull.begin(); 
       it!=hull.end() ; it++){
    if ( *it != l ){
      if ((*it)->intersect_line(lolline1) || (*it)->intersect_line(lolline2)) {
#ifdef DEBUG
        std::cout << "Lines intersect."<< std::endl;
        std::cout << std::endl;
#endif
        return false;
      }
    }
  }
#else
  // for (std::list<line>::iterator it = lines.begin(); 
  //      it!=lines.end() ; it++){
  //   if (&(*it) != l && 
  //       ((*it).intersect_line(lolline1) || (*it).intersect_line(lolline2)) ){
  //     std::cout << "Lines intersect."<< std::endl;
  //     return false;
  //   }
  // }
#endif

#if 1
// Area tests:
  Real area = calculate_triangle_area_2d(*n0, *n1, lolnode) ;
  Real equilateral = 0.43301270189 * avlen*avlen;

  if ( area < 0.2* equilateral||
       area > 2* equilateral){
    return false;
  }

#endif


// Remove from hull
  hull.erase(bine_it);


  nodes.push_front(lolnode);
  node *newnode = &(*nodes.begin());
//  node *newnode = &(*( nodes.end() ));
//  std::cout << &(nodes.back()) << std::endl;


// Add new lines to lines and hull
//  std::vector<line> newlines;
//  std::list<line>::reverse_iterator rline_it;
  lines.push_front(line(newnode,n1,true));
//  lines.push_front(line(newnode,n0,true));
  lines.push_front(line(n0,newnode,true));
  line_it=lines.begin();
//  line_it++;
  l0=& (*line_it);
  line_it++;
  l1=& (*line_it);

  triangles.push_front( triangle(l,l0,l1) );

  hull.push_front( l0 );
  hull.push_front( l1 );
//  std::cout<< "ASDASDAS "<< & lines.at(lines.size()-1)<< " " << newnode <<std::endl;
//  std::cout<< "ASDASDAS "<< & lines.at(lines.size()-2)<< " " << newnode <<std::endl;


//  newnode->id=999;
//  std::cout<<newnode->id << std::endl;


//  hull.erase(std::find(hull.begin(),hull.end(),l));
//  hull.erase(std::remove(hull.begin(), hull.end(), l), hull.end());
  return true;

}




void mesh::out_vtk1(bool ohull){
  FILE *outfile;
//  char vtkfilename[FILENAMESIZE];
  char vtkfilename[] = "out.vtk";
  unsigned int NEL;
  int nnodes = 2;
  int celltype = 3;
  if (ohull)
    NEL = hull.size();
  else
    NEL = lines.size();
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
  std::list<line*>::iterator bine_it;

  unsigned int NN = nodes.size();
/*
  if (ofilename != (char *) NULL && ofilename[0] != '\0') {
  strcpy(vtkfilename, ofilename);
  } else if (b->outfilename[0] != '\0') {
  strcpy(vtkfilename, b->outfilename);
  } else {
  strcpy(vtkfilename, "unnamed");
  }
  strcat(vtkfilename, ".vtk");
*/

  printf("Writing %s.\n", vtkfilename);

  outfile = fopen(vtkfilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("File I/O Error:  Cannot create file %s.\n", vtkfilename);
    return;
  }

  fprintf(outfile, "# vtk DataFile Version 2.0\n");
  fprintf(outfile, "Unstructured Grid\n");
  fprintf(outfile, "ASCII\n"); // BINARY
  fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(outfile, "POINTS %d double\n", NN);

  
  for(node_it=nodes.begin(); node_it!=nodes.end() ; node_it++){
    Real x = (*node_it)(0);//pointloop[0];
    Real y = (*node_it)(1);//pointloop[1];
    Real z = (*node_it)(2);//pointloop[2];

    fprintf(outfile, "%.6e %.6e %.6e\n", x, y, z);
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "CELLS %d %d\n", NEL, NEL*(nnodes+1));
  //NEL rows, each has 1 type id + 4 node id's
 
  if(ohull){
    // for(unsigned int i=0; i<NEL ; i++){
    //   unsigned int n1 = hull[i]->nodes[0]->id;
    //   unsigned int n2 = hull[i]->nodes[1]->id;
    //   fprintf(outfile, "%d  %4d %4d\n", nnodes, n1, n2);
    // }

    for(bine_it = hull.begin(); bine_it != hull.end() ; bine_it++){
      unsigned int n1 = (*bine_it)->nodes[0]->id;
      unsigned int n2 = (*bine_it)->nodes[1]->id;
      fprintf(outfile, "%d  %4d %4d\n", nnodes, n1, n2);
    }

  }else{
//  std::cout << "ASDASDASDASD" << std::endl;

    for(line_it = lines.begin(); line_it != lines.end() ; line_it++){
      unsigned int n1 = line_it->nodes[0]->id;
      unsigned int n2 = line_it->nodes[1]->id;
      fprintf(outfile, "%d  %4d %4d\n", nnodes, n1, n2);
    }

  }
  fprintf(outfile, "\n");
  
  fprintf(outfile, "CELL_TYPES %d\n", NEL);
  for(int tid=0; tid<NEL; tid++){
    fprintf(outfile, "%d\n", celltype);
  }
  fprintf(outfile, "\n");
  
  
  fclose(outfile);
  
}

void out_vtk1(list<mesh> &m, vector<vector<Real> > &r, char *lolchar){
  if (r.size() != m.size()){
    cerr << "out_vtk1: vector sizes do not match." << endl;
    return;
  }

  list<mesh>::iterator mesh_it;
  vector<vector<Real> >::iterator vec_it;
  unsigned int NEL=0;
  unsigned int NN=0;
  
  vector<unsigned int> NNN(1,0);
  vector<unsigned int> NELL(1,0);

  
  for (mesh_it = m.begin(); mesh_it != m.end(); mesh_it++){
    NN+=mesh_it->nodes.size();
    NEL+=mesh_it->triangles.size();
    NNN.push_back(NNN.back()+ mesh_it->nodes.size());
    NELL.push_back(NELL.back()+mesh_it->triangles.size());
  }


  FILE *outfile;
//  char vtkfilename[FILENAMESIZE];
//  char vtkfilename[] = "out.vtk";
  char *vtkfilename = lolchar;
  int nnodes = 3;
  int celltype = 5;

  list<node>::iterator node_it;
  list<triangle>::iterator tri_it;


  printf("Writing %s.\n", vtkfilename);

  outfile = fopen(vtkfilename, "w");
  if (outfile == (FILE *) NULL) {
    printf("File I/O Error:  Cannot create file %s.\n", vtkfilename);
    return;
  }

  fprintf(outfile, "# vtk DataFile Version 2.0\n");
  fprintf(outfile, "Unstructured Grid\n");
  fprintf(outfile, "ASCII\n"); // BINARY
  fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(outfile, "POINTS %d double\n", NN);
  
  for (mesh_it = m.begin(); mesh_it != m.end(); mesh_it++){
    
    for(node_it=mesh_it->nodes.begin(); node_it!=mesh_it->nodes.end() ; node_it++){
      Real x = (*node_it)(0);
      Real y = (*node_it)(1);
      Real z = (*node_it)(2);
      
      fprintf(outfile, "%.6e %.6e %.6e\n", x, y, z);
    }
    
  }
  fprintf(outfile, "\n");

  fprintf(outfile, "CELLS %d %d\n", NEL, NEL*(nnodes+1));

  unsigned i = 0;
  for (mesh_it = m.begin(); mesh_it != m.end(); mesh_it++, i++){
    
    for(tri_it = mesh_it->triangles.begin(); 
        tri_it != mesh_it->triangles.end() ; tri_it++){
      unsigned int n1 = tri_it->nodes[0]->id + NNN[i];
      unsigned int n2 = tri_it->nodes[1]->id + NNN[i];
      unsigned int n3 = tri_it->nodes[2]->id + NNN[i];
      fprintf(outfile, "%d  %4d %4d %4d\n", nnodes, n1, n2, n3);
    }
  }

  fprintf(outfile, "\n");
  
  fprintf(outfile, "CELL_TYPES %d\n", NEL);
  for(int tid=0; tid<NEL; tid++){
    fprintf(outfile, "%d\n", celltype);
  }
  fprintf(outfile, "\nPOINT_DATA %d\n\n",NN);
  fprintf(outfile,"SCALARS alpha float 1\nLOOKUP_TABLE default\n");
  
  for (vec_it = r.begin(); vec_it != r.end(); vec_it++){
    // if (vec_it->size() != NN) {
    //   cerr<< "Skipping scalar set with incompatible size." << endl;
    //   continue;
    // }
    for (vector<Real>::iterator it = vec_it->begin();
         it != vec_it->end(); it++){
      fprintf(outfile, "%6e\n", *it);
    }
  }

  fclose(outfile);

}

#if 0
bool mesh::add_node1(unsigned int ln, const Real probrad, const Real avlen){
  std::list<node>::iterator node_it;
  std::list<line>::iterator line_it;
  std::list<line*>::iterator bine_it;

  bine_it=hull.begin();
  std::advance(bine_it,ln);
  line *l=(*bine_it);//hull[ln];

  node *n0 = l->nodes[0];
  node *n1 = l->nodes[1];
  Real length = calculate_line_length_2d(*n0,*n1);
  line *l0, *l1;
  Point median, difference, poot, final,finam;
  median = (*n1 + *n0);//(l.nodes[1] + l.nodes[0]);
  for (unsigned int i=0; i<2;i++) median(i)/=2;
//  std::cout<<median;

  difference = *n1 - *n0;
  poot(0)=-1*difference(1); poot(1)= difference(0);
  for (unsigned int i=0; i<2;i++) poot(i)*=sqrt(3)/2;

  Point trolol = random_vector(probrad);
//  if (l->rh==true){
  final = median - poot + trolol;
//  }else{
  finam = median + poot + trolol;
//  }

   // final(0)+=0.05;
   // finam(0)+=0.05;

//  std::cout<<final << std::endl;
    
  node lolnode;
//  node newnode(final(0),final(1));
  node lolnode1(final(0),final(1));
  node lolnode2(finam(0),finam(1));
  if ( is_node_inside(lolnode1) ){
    if ( is_node_inside(lolnode2) )
      return false;
    else
//      std::cout << "false" << std::endl;
      lolnode = lolnode2;
  } else {
//    std::cout << "true" << std::endl;
    lolnode = lolnode1;

  }
// Intersection test    
  line lolline1(n1,&lolnode,true),lolline2(&lolnode,n0,true);


  // Real a[2]={-1,0};
  // Real b[2]={1,0};
  // Real c[2]={0,-1};
  // Real d[2]={0,0.001};

  // std::cout<<"ASDASDASDASD  "<< do_lines_intersect(a,b,c,d) << std::endl;
#if 1
  for (std::list<line*>::iterator it = hull.begin(); 
       it!=hull.end() ; it++){
    if ( *it != l &&
      ((*it)->intersect_line(lolline1) || (*it)->intersect_line(lolline2)) ){
#ifdef DEBUG
      std::cout << "Lines intersect."<< std::endl;
#endif
      return false;
    }
  }
#else
  for (std::list<line>::iterator it = lines.begin(); 
       it!=lines.end() ; it++){
    if (&(*it) != l && 
        ((*it).intersect_line(lolline1) || (*it).intersect_line(lolline2)) ){
      std::cout << "Lines intersect."<< std::endl;
      return false;
    }
  }


#endif

// Area tests:
  Real area = calculate_triangle_area_2d(*n0, *n1, lolnode) ;
//  Real equilateral = 0.43301270189 * length*length;
  
  Real equilateral = 0.43301270189 * avlen*avlen;

  if ( area < 0.2* equilateral||
       area > 2* equilateral){
    return false;
  }


// Remove from hull
  hull.erase(bine_it);


  nodes.push_front(lolnode);
  node *newnode = &(*nodes.begin());
//  node *newnode = &(*( nodes.end() ));
//  std::cout << &(nodes.back()) << std::endl;


// Add new lines to lines and hull
//  std::vector<line> newlines;
//  std::list<line>::reverse_iterator rline_it;
  lines.push_front(line(n1,newnode,true));
  lines.push_front(line(newnode,n0,true));
  line_it=lines.begin();
//  line_it++;
  l0=& (*line_it);
  line_it++;
  l1=& (*line_it);

  triangles.push_front( triangle(l,l0,l1) );

  hull.push_front( l0 );
  hull.push_front( l1 );
//  std::cout<< "ASDASDAS "<< & lines.at(lines.size()-1)<< " " << newnode <<std::endl;
//  std::cout<< "ASDASDAS "<< & lines.at(lines.size()-2)<< " " << newnode <<std::endl;


//  newnode->id=999;
//  std::cout<<newnode->id << std::endl;


//  hull.erase(std::find(hull.begin(),hull.end(),l));
//  hull.erase(std::remove(hull.begin(), hull.end(), l), hull.end());
  return true;

}

#endif
