

#ifndef __SECTION_H__ 
#define __SECTION_H__

#include "pr1.h"
#include <vector>
#include <ostream>



void sample_materials(supermesh &s, Real st, Real en, Real div, char *ofilename);



Real normal_force( supermesh &s, Real a, Real b, Real slope);
void normal_and_moment( supermesh &s, Real &F, Real &M, Real a, Real b, Real slope);
Real moment( supermesh &s, Real a, Real b, Real slope);

vector<Real> project2nodes( mesh &m, std::vector<Real> &vec);

//void interaction(supermesh &s);
vector<vector<Real> > interaction(supermesh &s, unsigned int ndiv, 
                                  char *ofilename, bool plotothers=false);

void m_vs_phi(supermesh &s, unsigned int ndiv, Real N, 
              char *ofilename, bool plotothers=false);

void plot_y_vs_sig( supermesh &s, Real a, Real b, Real slope, string lolchar);


#endif
// Local Variables:
// mode: c++
// End:
