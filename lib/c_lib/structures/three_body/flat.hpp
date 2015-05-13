#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__FLAT_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__FLAT_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/fct.hpp>
#include <meson_deca/lib/c_lib/complex.hpp>

namespace resonances {

  // Non-resonant (flat) resonance, 3-body decay
  struct flat_3 : resonances::resonance_base_3
  {
    // Constructor
    flat_3(particle _P, particle _a, particle _b, particle _c) : 
      resonance_base_3(_P, _a, _b, _c) {};

  
    // Returns 1 if we are within Dalitz plot bounds, 0 else.
    template <typename T>
    std::vector<T>
    value(const T& m2_ab, const T& m2_bc) {

      std::vector<T> res(2, 0.0);
      if (fct::valid(m2_ab, m2_bc, 
		     this->P, this->a, this->b, this->c) == true) {
	res[0] = 1.0;
      }
      return res;
    }
  };

}

#endif
