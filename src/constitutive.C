#include "constitutive.h"

namespace constitutive{
  void material::print(std::ostream& os) const{
    os<< "mat_param: (";
    for (std::vector<Real>::const_iterator it = par.begin(); it != par.end()-1; it++){
      os<< *it << ", ";
    }
    os << par.back() << ")";
  }



//////////////////////
// Concrete Follows //
//////////////////////


  Real concrete1::sig(Real eps){
//    if (eps > par[2]) throw 0;
    if (eps < par[2] || eps > par[3]) 
      return 0.;
    Real result;
    if (eps < 0){
      result = -1*par[0] * ((2*-1*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
    }else{
      result = 1 * ((2*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
    }
    return result;
  }
  concrete1::concrete1(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2){
    init(f_ck, eps_0, eps_cu1, eps_cu2);
  }

  concrete1::concrete1(){
    init(20,2e-3,-3.5e-3, 1e-3);
//    concrete1(20,2e-3,3.5e-3);
  }

  void concrete1::init(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2){
    par = std::vector<Real> (4);
    par[0]=(Real) f_ck;
    par[1]=(Real) eps_0;
//    par[2]=(Real) eps_cu;
    par[2] = (Real) -1*fabs(eps_cu1); // compression ultimate
    par[3] = (Real)    fabs(eps_cu2); // tension ultimate
    epsult[0] = par[2];
    epsult[1] = par[3];
    description="Concrete 1";
  }


//////////////////////
// Steel Follows    //
//////////////////////


  Real steel1::sig(Real eps){
//    if (eps > par[2]) throw 0;
    if (eps < par[2] || eps > par[3]) 
      return 0.;
    Real result;
    if (eps < 0){
      result = -1* 1 * ((2*-1*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
    }else{
      if (eps < par[1])
        result = eps*par[0];
      else
        result = par[1]*par[0];
    }
    return result;
  }
  steel1::steel1(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2){
    init(f_ck, eps_0, eps_cu1, eps_cu2);
  }

  steel1::steel1(){
    init(200000,0.0025,-3.5e-3, 0.01875);
//    steel1(20,2e-3,3.5e-3);
  }

  void steel1::init(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2){
    par = std::vector<Real> (4);
    par[0]=(Real) f_ck;
    par[1]=(Real) eps_0;
//    par[2]=(Real) eps_cu;
    par[2] = (Real) -1*fabs(eps_cu1); // compression ultimate
    par[3] = (Real)    fabs(eps_cu2); // tension ultimate
    epsult[0] = par[2];
    epsult[1] = par[3];
    description = "Steel 1";
  }

}
