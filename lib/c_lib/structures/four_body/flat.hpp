#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__FLAT_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__FLAT_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/fct.hpp>
#include <meson_deca/lib/c_lib/complex.hpp>
#include <meson_deca/lib/c_lib/structures/four_body/base.hpp> // base class

namespace resonances {

  // Non-resonant (flat) resonance, 5-dimensional 4 body decay.
  struct flat_4 : resonances::resonance_base_4
  {
    const particle R_1;
    const particle R_2;

    // Constructor
    flat_4(particle _P, particle _R_1, particle _R_2,
	   particle _a, particle _b, particle _c, particle _d) : 
      resonance_base_4(_P, _a, _b, _c, _d), R_1(_R_1), R_2(_R_2) {};

    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename 
		boost::math::tools::promote_args<T0, T1, T2, T3, T4>::type>
    value(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
	  const T3 &m2_34, const T4& m2_13)
    {

      typedef typename 
	boost::math::tools::promote_args<T0, T1, T2, T3, T4>::type T_res;

      T_res m2_123;
      m2_123 = m2_12 + m2_13 + m2_23 - this->a.m2 - this->b.m2 - this->c.m2;

      std::vector<T_res> res(2, 0.0);

      if (fct::valid_5d(m2_12, m2_14, m2_23, m2_34, m2_13, m2_123,
			this->P, this->R_1, this->R_2,
			this->a, this->b, this->c, this->d) == true)
	{  
	  res[0] = 1.0;
	}		

      return res;
    }
  };

}

#endif
