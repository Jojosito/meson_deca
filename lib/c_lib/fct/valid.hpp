#ifndef MESON_DECA__LIB__C_LIB__FCT__VALID_HPP
#define MESON_DECA__LIB__C_LIB__FCT__VALID_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/structures/struct_particles.hpp> 
// class particle

namespace fct {

  template <typename T>
  inline
  bool valid(const T &m2_ab, const T &m2_bc, 
             const particle &p, const particle &a, 
             const particle &b, const particle &c)
  {
    if ( (m2_ab < (a.m2 + b.m2 + 2. * sqrt(a.m2 * b.m2))) ||
         (m2_ab > (p.m2 + c.m2 - 2. * sqrt(p.m2 * c.m2)))   ) {
      return false;
    }

   typedef typename boost::math::tools::promote_args<T, double>::type T_res;

    T_res E_b = (m2_ab - a.m2 + b.m2) / 2. / sqrt(m2_ab);
    T_res E_c = (p.m2 - m2_ab - c.m2) / 2. / sqrt(m2_ab);
    T_res P_b = sqrt(E_b * E_b - b.m2);
    T_res P_c = sqrt(E_c * E_c - c.m2);

    if ((fabs(m2_bc - b.m2 - c.m2 - 2. * E_b * E_c)) <= (2. * P_b * P_c)){
      return true;
    }

    return false;
  }

}
#endif
