

#include "interval.h"

compoundInterval compoundInterval::operator+ (compoundInterval c){
  compoundInterval v=*this;
  for (unsigned int i = 0; i< c.size(); i++)
    v+=c[i];
  return v;
}
compoundInterval compoundInterval::operator- (compoundInterval c){
  compoundInterval v=*this;
  for (unsigned int i = 0; i< c.size(); i++)
    v-=c[i];
  return v;
  
}
compoundInterval& compoundInterval::operator+= (compoundInterval c){
  for (unsigned int i = 0; i< c.size(); i++)
    *this+=c[i];
  return *this;
}
compoundInterval& compoundInterval::operator-= (compoundInterval c){
  for (unsigned int i = 0; i< c.size(); i++)
    *this-=c[i];
  return *this;
}


compoundInterval compoundInterval::operator+ (interval c){
  if (numtk::iszero(c.e - c.s)) return *this;
  compoundInterval v=*this;
  if (_intervals.size()==0){
    v._intervals.push_back(c); 
    return v;
  }else{
    v._intervals.push_back(c);
    v.yunion();
    return v;
  }

}
compoundInterval compoundInterval::operator- (interval c){
  if (numtk::iszero(c.s - c.e)) return *this;
  interval d;

  sort(_intervals.begin(), _intervals.end(), interval_sorter);

  
  vector<interval> r;
  for (unsigned int i = 0; i < _intervals.size(); i++){
    if (c.e <= _intervals[i].s){
      r.push_back(_intervals[i]);
    }else if (c.s <= _intervals[i].s && c.e >= _intervals[i].s && c.e <= _intervals[i].e){
      d.s = c.e; d.e = _intervals[i].e; r.push_back(d);
    } else if (c.s <= _intervals[i].s && c.e >= _intervals[i].e){
      continue;
    } else if (c.s >= _intervals[i].s && c.e <= _intervals[i].e){
      d.s=_intervals[i].s; d.e=c.s; r.push_back(d);
      d.s=c.e; d.e=_intervals[i].e; r.push_back(d);
    } else if (c.s >= _intervals[i].s && c.s <= _intervals[i].e && c.e >= _intervals[i].e){
      d.s = _intervals[i].s; d.e = c.s; r.push_back(d);
    } else if (c.s >= _intervals[i].e){
      r.push_back(_intervals[i]);
    }
  }

  compoundInterval result;

  for (unsigned int i = 0; i < r.size(); i++)
    if (r[i].s!=r[i].e && ~numtk::iszero(r[i].s - r[i].e))
        result._intervals.push_back(r[i]);

  return result;
}

void compoundInterval::yunion(){
  if (_intervals.size()==0) return;
  sort(_intervals.begin(), _intervals.end(), interval_sorter);
  
  vector<interval> r;
  r.push_back(_intervals[0]);
  for (unsigned int i = 0; i < _intervals.size(); i++){

    if (r.back().e < _intervals[i].s)
      r.push_back(_intervals[i]);
    else if (r.back().e == _intervals[i].s)
      r.back().e = _intervals[i].e;
    if (r.back().e < _intervals[i].e)
      r.back().e = _intervals[i].e;
  }
  _intervals=r;
}

compoundInterval& compoundInterval::operator+= (interval c){
  if (numtk::iszero(c.e - c.s)) return *this;

  if (_intervals.size()==0){
    _intervals.push_back(c); 
    return *this;
  }else{
    _intervals.push_back(c);
    yunion();
    return *this;
  }

}
compoundInterval& compoundInterval::operator-= (interval c){
  if (numtk::iszero(c.s - c.e)) return *this;
  interval d;

  sort(_intervals.begin(), _intervals.end(), interval_sorter);

  
  vector<interval> r;
  for (unsigned int i = 0; i < _intervals.size(); i++){
    if (c.e <= _intervals[i].s){
      r.push_back(_intervals[i]);
    }else if (c.s <= _intervals[i].s && c.e >= _intervals[i].s && c.e <= _intervals[i].e){
      d.s = c.e; d.e = _intervals[i].e; r.push_back(d);
    } else if (c.s <= _intervals[i].s && c.e >= _intervals[i].e){
      continue;
    } else if (c.s >= _intervals[i].s && c.e <= _intervals[i].e){
      d.s=_intervals[i].s; d.e=c.s; r.push_back(d);
      d.s=c.e; d.e=_intervals[i].e; r.push_back(d);
    } else if (c.s >= _intervals[i].s && c.s <= _intervals[i].e && c.e >= _intervals[i].e){
      d.s = _intervals[i].s; d.e = c.s; r.push_back(d);
    } else if (c.s >= _intervals[i].e){
      r.push_back(_intervals[i]);
    }
  }

  _intervals.clear();

  for (unsigned int i = 0; i < r.size(); i++)
    if (r[i].s!=r[i].e && ~numtk::iszero(r[i].s - r[i].e))
        _intervals.push_back(r[i]);

  return *this;
}

compoundInterval& compoundInterval::operator= (const compoundInterval& c){
  this->_intervals = c._intervals;
  return *this;
}


compoundInterval compoundInterval::intersection (compoundInterval c)
{
  this->yunion(); c.yunion();
  compoundInterval r;
  interval d;
  for (unsigned int i = 0; i < this->size(); i++){
    for (unsigned int j = 0; j < c.size(); j++){
      if ((*this)[i].e < c[i].s || (*this)[i].s > c[i].e){
        continue;
      }else if ( (*this)[i].s <= c[i].s && (*this)[i].e >= c[i].s && (*this)[i].e <= c[i].e){
        d.s = c[i].s; d.e = (*this)[i].e; r._intervals.push_back(d);
      }else if ( (*this)[i].s <= c[i].s && (*this)[i].e >= c[i].e ){
        r._intervals.push_back(c[i]);
      }else if ( (*this)[i].s >= c[i].s && (*this)[i].e <= c[i].e ){
        r._intervals.push_back((*this)[i]);
      }else if ( (*this)[i].s >= c[i].s && (*this)[i].s <= c[i].e && (*this)[i].e >= c[i].e){
        d.s = (*this)[i].s; d.e = c[i].e; r._intervals.push_back(d);
      }
    }
  }
  compoundInterval s;
  for (unsigned int i = 0; i < r.size(); i++)
    if (r[i].s!=r[i].e && ~numtk::iszero(r[i].s - r[i].e))
        s._intervals.push_back(r[i]);

  return s;

}


void compoundInterval::subtract_infinity(Real a, bool side)
{
  interval d;

  // side == 1 for positive, == 0 for negative
  
  compoundInterval r;
  if (side == true){
    for (unsigned int i = 0; i < this->size(); i++){
      if (a <= (*this)[i].s && a <= (*this)[i].e){
        continue;
      }else if (a >= (*this)[i].s && a <= (*this)[i].e){
        d.s = (*this)[i].s; d.e = a; r._intervals.push_back(d);
      } else if (a >= (*this)[i].s && a >= (*this)[i].e){
        r._intervals.push_back((*this)[i]);
      }
    }
  }else {
    for (unsigned int i = 0; i < this->size(); i++){
      if (a <= (*this)[i].s && a <= (*this)[i].e){
        r._intervals.push_back((*this)[i]);
      }else if (a >= (*this)[i].s && a <= (*this)[i].e){
        d.s = a; d.e = (*this)[i].e; r._intervals.push_back(d);
      } else if (a >= (*this)[i].s && a >= (*this)[i].e){
        continue;
      }
    }
  }
  this->clear();
  for (unsigned int i = 0; i < r.size(); i++)
    if (r[i].s!=r[i].e && ~numtk::iszero(r[i].s - r[i].e))
        this->_intervals.push_back(r[i]);


}

