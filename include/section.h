

#ifndef __SECTION_H__ 
#define __SECTION_H__

#include "pr1.h"



void single_rein_rect(supermesh &s, Real width, Real height, 
                 unsigned int nrein, Real rad, Real dprime);

void double_rein_rect(supermesh &s, Real width, Real height,
                 unsigned int nrein, Real rad, Real dprime, Real dprpr);

void lateral_reinforcements(mesh &m, unsigned int n, 
                            Real y, Real width, Real radius,
                            unsigned int nsides);


Real normal_force( supermesh &s, Real a, Real b, Real slope);
void normal_and_moment( supermesh &s, Real F, Real M, Real a, Real b, Real slope);
Real moment( supermesh &s, Real a, Real b, Real slope);

vector<Real> project2nodes( mesh &m, std::vector<Real> &vec);

void interaction(supermesh &s);

void interaction_1(supermesh &s);
void interaction_2(supermesh &s);

void plot_y_vs_sig( supermesh &s, Real a, Real b, Real slope, string lolchar);

#endif
// Local Variables:
// mode: c++
// End:
