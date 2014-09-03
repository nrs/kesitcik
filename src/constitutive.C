#include "kesitcik/constitutive.h"

namespace constitutive{

  void material::print(std::ostream& os) const{
    os<< "mat_param: (";
    for (std::vector<Real>::const_iterator it = par.begin(); it != par.end()-1; it++){
      os<< *it << ", ";
    }
    os << par.back() << ")";
  }


}
