#ifndef MESON_DECA__LIB__C_LIB__FCT__VALID_HPP
#define MESON_DECA__LIB__C_LIB__FCT__VALID_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/structures/struct_particles.hpp> 
// class particle

/*
 * Check whether we are in the energetically allowed region of the decay.
 * 
 * FUNCTIONS
 * valid(m2_ab, m2_bc, p, a, b, c) - Check m2_ab, m2_bc for p->abc decay
 * valid(m2_ab, m2_bc, m2_p, m2_a, m2_b, m2_c) - As above, other input type
 * valid_5d() // WRONG. To be corrected.
 * 
 */

namespace fct {

  /**
   * bool valid(m2_ab, m2_bc, p, a, b, c)
   *
   * Determines whether we are in an energetically allowed
   * phase space region of the decay p-> a + b + c.
   *
   * m2_ab, m2_bc are the invariant square masses of a and b, b and c 
   * (scalars). p, a, b, c are decay particles, usually specified in
   * c_lib/structures/particles.hpp .
   *
   * Since we do not need the momenta of p, a, b, c, we can use their
   * masses as arguments - the overloaded function may be found below.
   * 
   */
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


  /*
   * Overloaded input arguments to allow particle masses.
   */
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline
  bool valid(const T0 &m2_ab, const T0 &m2_bc, 
             const T1 &m2_p, const T2 &m2_a, const T3 &m2_b, const T4 &m2_c)
  {
    if ( (m2_ab < (m2_a + m2_b + 2. * sqrt(m2_a * m2_b))) ||
         (m2_ab > (m2_p + m2_c - 2. * sqrt(m2_p * m2_c)))   ) {
      return false;
    }

    typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4,double>::type T_res;

    T_res E_b = (m2_ab - m2_a + m2_b) / 2. / sqrt(m2_ab);
    T_res E_c = (m2_p - m2_ab - m2_c) / 2. / sqrt(m2_ab);
    T_res P_b = sqrt(E_b * E_b - m2_b);
    T_res P_c = sqrt(E_c * E_c - m2_c);

    /* verbose debugging
    std::cout << "Valid range: from " 
	      << (E_b + E_c) * (E_b + E_c) - (P_b + P_c)
	      << " to "
	      << (E_b + E_c) * (E_b + E_c) - (P_b - P_c)
	      << "\n";
    */

    if ((fabs(m2_bc - m2_b - m2_c - 2. * E_b * E_c)) <= (2. * P_b * P_c)){
      return true;
    }

    return false;
  }


  /**
   * bool valid_5d(m2_12, ..., m2_123, p, a, b, c)
   *
   * Determines whether we are in an energetically allowed
   * phase space region of the decay P -> R_1 d -> R_2 c d -> a b c d.
   *
   * The variable m2_123 is passed, although it can also be derived 
   * from the invariant masses before that. We pass it to save one
   * computation line; maybe it does not really matter (sigh).
   */
  template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
  bool valid_5d(const T0 &m2_12, const T1 &m2_14, const T2 &m2_23,
		const T3 &m2_34, const T4& m2_13, const T5 &m2_123,
		const particle &P, const particle &R_1, const particle &R_2,
		const particle &a, const particle &b, 
		const particle &c, const particle &d)
  {
    typedef typename boost::math::tools::promote_args<T0,T2,T4>::type T_123;
    typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

    // 3 body decay  P -> R_2 c d
    if (fct::valid(m2_123, m2_34, P.m2, R_2.m2, c.m2, d.m2) == false) {
      /*std::cout << "Called valid(" << m2_123
		<< ", " << m2_34
		<< ", " << P.m2
		<< ", " << m2_12
		<< ", " << c.m2
		<< ", " << d.m2
		<< ").\n";*/
      return false;
    }

    // 3 body decay R_1 -> a b c
    if (fct::valid(m2_12, m2_23, R_1.m2, a.m2, b.m2, c.m2) == false) {
      // verbose debugging
      //std::cout << "Call 2.";
      return false;
    }

    return true;
  }

}
#endif
