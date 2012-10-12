

#ifndef __INTERVAL_H__ 
#define __INTERVAL_H__

#include <vector>
#include <ostream>
#include "pr1.h"



class interval{
public:
  Real s,e;
  interval(Real a, Real b);
  interval(){return;}
};

inline interval::interval(Real a, Real b){
  s=min(a,b); e=max(a,b);
}

inline bool interval_sorter(interval const& lhs, interval const& rhs) {
  return lhs.s < rhs.s;
}

class compoundInterval{
public:
  void add(Real a, Real b);
  void subtract(Real a, Real b);
  
  compoundInterval operator+ (interval c);
  compoundInterval operator- (interval c);
  compoundInterval& operator+= (interval c);
  compoundInterval& operator-= (interval c);

  compoundInterval operator+ (compoundInterval c);
  compoundInterval operator- (compoundInterval c);
  compoundInterval& operator+= (compoundInterval c);
  compoundInterval& operator-= (compoundInterval c);


  compoundInterval intersection (compoundInterval c);
  void subtract_infinity(Real a, bool side);

  inline interval& operator[] (unsigned int i){
    return _intervals[i];
  }

  inline void clear(){ _intervals.clear();}
  inline unsigned int size() { return _intervals.size();}

  compoundInterval& operator= (const compoundInterval& c);  
  void print(std::ostream& os = std::cout) const{
    for (vector<interval>::const_iterator i=_intervals.begin(); 
         i!=_intervals.end(); i++){
      os << "(" << i->s << ", " << i->e << ") ";
    }
  }
  friend std::ostream& operator << (std::ostream& os, const compoundInterval& t){
    t.print(os);return os;
  }

  void yunion();
  vector<interval> _intervals;
};


#endif
// Local Variables:
// mode: c++
// End:
