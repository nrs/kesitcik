

#ifndef __SECTION_H__ 
#define __SECTION_H__

#include "pr1.h"
#include <vector>
#include <ostream>

void single_rein_rect(supermesh &s, Real width, Real height, 
                 unsigned int nrein, Real rad, Real dprime);

void double_rein_rect(supermesh &s, Real width, Real height,
                      unsigned int nrein1, Real rad1, Real dprime1,
                      unsigned int nrein2, Real rad2, Real dprime2);

void triple_rein_rect(supermesh &s, Real width, Real height,
                      unsigned int nrein1, Real rad1, Real dprime1,
                      unsigned int nrein2, Real rad2, Real dprime2,
                      unsigned int nrein3, Real rad3, Real dprime3);

void lateral_reinforcements(mesh &m, unsigned int n, 
                            Real y, Real width, Real radius,
                            unsigned int nsides);


Real normal_force( supermesh &s, Real a, Real b, Real slope);
void normal_and_moment( supermesh &s, Real &F, Real &M, Real a, Real b, Real slope);
Real moment( supermesh &s, Real a, Real b, Real slope);

vector<Real> project2nodes( mesh &m, std::vector<Real> &vec);

//void interaction(supermesh &s);
vector<vector<Real> > interaction(supermesh &s, Real ranst, Real ranen, 
                                  unsigned int ndiv);

void plot_y_vs_sig( supermesh &s, Real a, Real b, Real slope, string lolchar);


#endif
// Local Variables:
// mode: c++
// End:
