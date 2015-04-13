#ifndef MESON_DECA__LIB__C_LIB__FCT__ZEMACH_HPP
#define MESON_DECA__LIB__C_LIB__FCT__ZEMACH_HPP

#include <vector>
#include <cmath>
#include <complex>

#include <meson_deca/lib/c_lib/fct.hpp>

namespace fct {

  template <typename T0, typename T1, typename T2>
  typename boost::math::tools::promote_args<T0,T1,T2>::type
  zemach(int J, const T0& m2_ab, const T1& m2_bc, const T2& m2_R,
         const particle &a, const particle &b, const particle &c) {

    if (J == 0) return 1;
    if (J == 1) {
      return m2_R + a.m2 + b.m2 + c.m2 - m2_ab - 2 * m2_bc - 
	(m2_R - c.m2) * (a.m2 - b.m2) / m2_ab;
    }

    if (J == 2) {
      return fct::zemach(1, m2_ab, m2_bc, m2_R, a, b, c) -
	(m2_ab - 2.*m2_R - 2.*c.m2 + (m2_R - c.m2) * (m2_R - c.m2) / m2_ab) *
	(m2_ab - 2.*a.m2 - 2.*b.m2 + (a.m2 - b.m2) * (a.m2 - b.m2) / m2_ab) /
        3.0;
    }

    return 0;
  }

}

#endif
