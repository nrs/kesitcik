#ifndef __PR1_H__
#define __PR1_H__

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <queue>
#include <algorithm>
#include <utility>
#include <iterator>
#include <memory>
#include <limits>
#include <stdio.h>
#include <time.h>
#include "point.h"
#include "constitutive.h"
#include <cmath>
#include <cstring>

#define VOLUME_TOLERANCE 1e-4

class node;
class line;
class triangle;
class mesh;

using namespace std;

unsigned int random_range(unsigned int min, unsigned int max);
Point random_vector(Real rad);

Real calculate_triangle_area_2d(node &aa, node &bb, node &cc);
Real calculate_triangle_area_2d(Real a[2], Real b[2], Real c[2]);


namespace numtk{
  std::vector<Real> range(Real st, Real en, unsigned int a);
  std::string exec(char* cmd);

  void output_gp(std::vector<Real> &a, std::vector<Real> &b, string filename);
  void output_gp(std::vector<Real> &a, 
                 std::vector<Real> &b, std::vector<Real> &c, string filename);

  Real distance(node a, Real m, Real n);

  template<typename T> inline bool isnan(T value)
  {
    value=fabs(value);
    return value != value;
  }

// requires #include <limits>
  template<typename T> inline bool isinf(T value)
  {
    return std::numeric_limits<T>::has_infinity &&
      fabs(value) == std::numeric_limits<T>::infinity();
  }

  template<typename T> inline bool iszero(T value)
  {
    return (fabs(value) < TOLERANCE);
  }

  template <typename T> inline int sign(T val) {
    return (T(0) < val) - (val < T(0));
  }
  template <typename T> inline bool isposit(T val) {
    return val >= 0;
  }
}

inline Real random_Real(){
  return ((Real)std::rand()/(Real)RAND_MAX);

}
class node : public Point
{
public:
  unsigned int id;
  node(const Real x, const Real y);
  node();
  unsigned int color;

};


inline node::node(){
  node(0,0);
}


inline node::node(const Real x, const Real y):Point(x,y,0){
  // _coords[0]=x;
  // _coords[1]=y;
  // _coords[2]=0;
  
}

class line 
{
public:
  node *nodes[2];
//  line(node *n1, node *n2);
  line(node *n1, node *n2, const bool rh=true);
  bool intersect_line(line &l);
//  line(node ns[2]);
  bool operator==(const line &l);
  bool rh;
  unsigned int id;
  unsigned int color;
};



class triangle
{
public:
  node *nodes[3];
  line *lines[3];
  node *cent;
  triangle(line *l1, line *l2, line *l3);
  triangle(node *n1, node *n2, node *n3, mesh &m); // don't use this generally
  triangle(node *n1, node *n2, node *n3); 

  bool is_node_inside(node &n);
  bool operator==(const line &l);
  unsigned int id;
  Real area;
  Real calc_area();
  unsigned int color;

};



inline void set_color(list<line> &l, unsigned int c){
  for (list<line>::iterator it=l.begin() ; it != l.end() ; it++)
    it->color = c;
}

inline void set_color(list<node> &l, unsigned int c){
  for (list<node>::iterator it=l.begin() ; it != l.end() ; it++)
    it->color = c;
}

inline void set_color(list<triangle> &l, unsigned int c){
  for (list<triangle>::iterator it=l.begin() ; it != l.end() ; it++)
    it->color = c;
}



inline void add_line(node a, node b, list<node> &ln, list<line> &ll){
  node *c[2];
  ln.push_back(a); c[0] = &ln.back();
  ln.push_back(b); c[1] = &ln.back();
  ll.push_back(line(c[0],c[1]));

}

using constitutive::material;
class mesh
{
public:
  
  // unsigned int nnodes;
  // Real avlen, probrad;

  std::list<node> nodes;
  std::list<node> centroids;

  std::list<line> lines;
  std::list<line*> hull;
  std::list<triangle> triangles;

  std::map<node*, std::set<triangle*> > node2tri;
  void init_node2tri();

//  auto_ptr<material> mat(new material());
  material *mat;

  bool in_vtk_lines(char *infilename);

  list<node> holes;
  list<node> innode;

  // void generate1(unsigned int n, Real a, Real p);
  // void generate1();

  void numbernodes();
  line* get_random_hull_elem();
  // bool add_node1(unsigned int ln, const Real probrad, const Real avlen);
  // bool add_node2(unsigned int ln, const Real probrad, const Real avlen);

  void print_info();
  void out_vtk1(bool ohull=false);

  // void generate_seed_triangle(Real avlen, Real probrad);
  bool is_node_inside(node &n);
  bool intersect_line(line &l);
  void clear();
  void triangulate_mesh(Real minarea);

  void init_centroids();
  void subtract(mesh &m);
  void add(mesh &m);

  Real area;
  Real calc_area();
  Real recalc_area();
// center of gravity
  node cg;
  node calc_cg();

  void color_nodes(unsigned int c);
  void color_lines(unsigned int c);
  void color_triangles(unsigned int c);
  void color_mesh(unsigned int a,unsigned int b,unsigned int c);

  // mesh(unsigned int n, Real a, Real p);
  mesh(){mat = NULL;}
  ~mesh(){
    if (mat != NULL)
      delete mat;
  }

};

bool do_lines_intersect_exclusive(Real a[2], Real b[2], Real c[2], Real d[2]);
bool do_lines_intersect_inclusive(Real a[2], Real b[2], Real c[2], Real d[2]);

/*
template <typename I> I random_element(I begin, I end)
{
    const unsigned long n = std::distance(begin, end);
    const unsigned long divisor = (RAND_MAX + 1) / n;
    I result=begin;
    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);
    std::advance(result, k);
    return result;//*std::advance(begin, k);
}
*/


class supermesh
{
public:
  list<mesh> meshes;

  list<mesh> miscmeshes;
  list<node> miscnodes;
  list<line> misclines;
  list<triangle> misctris;

// center of gravity
  node cg;
  node calc_cg();

  Real area;
  Real calc_area();

};


template <class InputIterator> InputIterator 
random_n(InputIterator first, InputIterator last) {
   typename std::iterator_traits<InputIterator>::difference_type distance = 
        std::distance(first, last);
   InputIterator result = first;
   const unsigned long divisor = (RAND_MAX + 1) /distance;

   std::cout<<divisor<<std::endl;
   if (distance > 1) {
      // Uses std::rand() naively.  Should replace with more uniform solution. 
      std::advance( result, std::rand() % distance );
   }
   return result;
}

inline Real calculate_line_length_2d(node &aa, node &bb){
  Real a[2],b[2];
  for (unsigned int i = 0; i<2; i++) {
    a[i]=aa(i); b[i]=bb(i); 
  }
  return sqrt( (a[0]-b[0]) * (a[0]-b[0]) + (a[1]-b[1]) * (a[1]-b[1]) );
}

void rectangular_section(mesh &m, Real left, Real right, Real top, Real bottom);
void circular_section(mesh &m, Real rad, Real cen[2],int nsides);
void out_vtk1(list<mesh> &m, vector<vector<Real> > &r, char *lolchar);

void out_vtk1(list<mesh> &m, vector<vector<Real> > &r1, 
              vector<vector<Real> > &r2, char *lolchar);



// template <class I> I&
// listelem(unsigned int i) {
//    typename std::iterator_traits<I>::difference_type distance = 
//         std::distance(first, last);
//    I result = first;
//    const unsigned long divisor = (RAND_MAX + 1) /distance;

//    std::cout<<divisor<<std::endl;
//    if (distance > 1) {
//       // Uses std::rand() naively.  Should replace with more uniform solution. 
//       std::advance( result, std::rand() % distance );
//    }
//    return result;
// }

/*
unsigned int get_random_from_range(unsigned int beg, unsigned int end) {
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(beg,end); // define the range
    return (distr(eng));
//    for(int n=0; n<40; ++n)
//        std::cout << distr(eng) << ' '; // generate numbers

}
*/





// inline line* mesh::get_random_hull_elem(){
//   unsigned int ind = random_range(0,hull.size());
// //  std::cout<<hull.size()<<std::endl;
//   return (hull[ ind ]);
// }


#endif
// Local Variables:
// mode: c++
// End:
