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
    if (eps < epsulttol[0] || eps > epsulttol[1]){
      return 0.;
    }
    Real result;
    if (eps < 0){
      result = -1*par[0] * ((2*-1*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
    }else{
      result = 1.6 * ((2*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
    }
    return result;
  }
  concrete1::concrete1(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2){
    init(f_ck, eps_0, eps_cu1, eps_cu2);
  }

  concrete1::concrete1(){
    init(16,2e-3,-0.003, 0.003);

//    init(20,2e-3,-3.5e-3, 1e-3);
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
    epsulttol[0] = epsult[0]-TOLERANCE; 
    epsulttol[1] = epsult[1]+TOLERANCE;
  }


//////////////////////
// Steel Follows    //
//////////////////////


  Real steel1::sig(Real eps){
//    if (eps > par[2]) throw 0;
    if (eps < epsulttol[0] || eps > epsulttol[1])
      return 0.;
    Real result;
    if (eps < 0){
      if (eps > par[2])
        result = eps*par[0];
      else
        result = par[2]*par[0];

    }else{
      if (eps < par[3])
        result = eps*par[0];
      else
        result = par[3]*par[0];
    }
    return result;
  }
  steel1::steel1(Real f_ck, Real eps_0, Real eps_cu1, Real eps_cu2,
                 Real eps_u1, Real eps_u2){
    init(f_ck, eps_0, eps_cu1, eps_cu2, eps_u1, eps_u2);
  }

  steel1::steel1(){
//    init(200000,0.0025,-3.6e-3, 1.1e-3);
    init(200000,0.0025,-0.001825, 0.001825,-0.35,0.35);

  }

  void steel1::init(Real f_ck, Real eps_0, Real eps_sy1, Real eps_sy2,
                    Real eps_u1, Real eps_u2){
    par = std::vector<Real> (4);
    par[0]=(Real) f_ck;
    par[1]=(Real) eps_0;
//    par[2]=(Real) eps_cu;
    par[2] = (Real) -1*fabs(eps_sy1); // pos yield
    par[3] = (Real)    fabs(eps_sy2); // neg yield
    epsult[0] = eps_u1; // pos ultimate
    epsult[1] = eps_u2; // neg ultimate
    description = "Steel 1";
    epsulttol[0] = epsult[0]-TOLERANCE; 
    epsulttol[1] = epsult[1]+TOLERANCE;
  }

}
