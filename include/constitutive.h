#ifndef __CONSTITUTIVE_H__
#define __CONSTITUTIVE_H__

#include "point.h"
#include <vector>
#include <iostream>
#include <ostream>
#include <string>

using std::ostream;

namespace constitutive{

  class material{
  public:
    virtual Real sig(Real eps){return 0;}
    Real tanmod_f(Real eps, Real h){return (sig(eps+h)-sig(eps))/h;}
    Real tanmod_b(Real eps, Real h){return (sig(eps)-sig(eps-h))/h;}
    Real tanmod_c(Real eps, Real h){return (sig(eps+h)-sig(eps-h))/2/h;}
    std::vector<Real> par;
    Real epsult[2];
    Real epsulttol[2];
    std::string description;
//    std::vector<Real> eps_ult;
    void print(std::ostream& os = std::cout) const;
//    friend ostream& operator<< (ostream &out, material &mat);
    friend std::ostream& operator << (std::ostream& os, const material& t){
      t.print(os);return os;
    }
  };




}
#endif
// Local Variables:
// mode: c++
// End:
