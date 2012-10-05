

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

#endif
// Local Variables:
// mode: c++
// End:
