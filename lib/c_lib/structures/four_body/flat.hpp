#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__FLAT_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__FLAT_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/fct.hpp>
#include <meson_deca/lib/c_lib/complex.hpp>

namespace resonances {
  struct flat // Non-resonant (flat) resonance, 3-body decay
  {    
    double dummy_member; // I don't know, but without members
    // that does not compile. C++ is hard.

    flat(double dummy) : dummy_member(dummy) {}; // Constructor

    template <typename T>
    std::vector<T>
    value(const T& m2_ab, const T& m2_bc,
          const particle &p, const particle &a, 
          const particle &b, const particle &c) {

      std::vector<T> res(2, 0.0);
      if (fct::valid(m2_ab, m2_bc, p, a, b, c) == true) {
	res[0] = 1.0;
      }
      return res;
    }
  };


  struct flat_5d // Non-resonant (flat) resonance, 5-dimensional 4 body decay.
  {    
    double dummy_member; // I don't know, but without members
    // that does not compile. C++ is hard.

    flat_5d(double dummy) : dummy_member(dummy) {}; // Constructor

    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T3, T4>::type>
    value(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
	  const T3 &m2_34, const T4& m2_13,
	  const particle &P, const particle &a, 
	  const particle &b, const particle &c, const particle &d)
    {

      typedef typename boost::math::tools::promote_args<T0, T1, T2, T3, T4>::type T_res;

      T_res m2_123;
      m2_123 = m2_12 + m2_13 + m2_23 - a.m2 - b.m2 - c.m2;

      std::vector<T_res> res(2, 0.0);
      if (fct::valid_5d(m2_12, m2_14, m2_23, m2_34, m2_13, m2_123,
			P, a, b, c, d) == true) {
	res[0] = 1.0;
      }
      return res;
    }
  };
}

#endif
