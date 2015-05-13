#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__D_R1D_R2cd_abcd_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__D_R1D_R2cd_abcd_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/fct.hpp> // Breit-Wigner, Blatt-Weisskopf, etc.
#include <meson_deca/lib/c_lib/complex.hpp> // Complex numbers

#include <meson_deca/lib/c_lib/structures/four_body/base.hpp> // base class


namespace resonances {

  // Amplitude function for the 4 particle decay
  //   P -> R_1 d -> R_2 c d -> a b c d
  //
  // Model-dependent description is glued from
  // two Breit-Wigner 3-body decays
  //   P   -> R_2 c d
  //   R_1 -> a   b c
  //
  // For detailed description, see
  //   https://www.overleaf.com/2630893ckqcxt#/6945226/
  // or feel free to contact me at arseniy.tsipenyuk@gmail.com.
  struct P_R1d_R2cd_abcd : resonances::resonance_base_4
  {

    const int l_1; // Orbital angular momentum between R_1 and d
    const int l_2; // Orbital angular momentum between R_2 and c
    const int l_3; // Orbital angular momentum between a and b
    const particle R_1; // First decay resonance (e.g. a_1)
    const particle R_2; // 2nd order decay resonance (e.g. rho_0)
    const double W_R_1, W_R_2; // Width of the 1st, 2nd resonance
 

    // Default constructor
    P__R1d__R2cd__abcd__bw_decay(particle _P, particle _a, particle _b, 
				 particle _c, particle _d,
				 int _l_1, int _l_2, int _l_3,  
				 particle _R_1, particle _R_2,
				 double _W_R_1, double _W_R_2) : 
      decay_base(_P,_a, _b, _c, _d), 
      l_1(_l_1), l_2(_l_2), l_3(_l_3),
      R_1(_R_1), R_2(_R_2), W_R_1(_W_R_1), W_R_2(_W_R_2) {};


    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABCD (not symmetrized)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    // m2_12 is the invariant square mass of particles a and b.
    // Analogously, m2_34 is i.sq.m. of c and d, m2_23 - of b and c, etc.
    value(const T0& m2_12, const T1& m2_14, const T2& m2_23,
	  const T3& m2_34, const T4& m2_13) {

      typedef typename boost::math::tools::promote_args<T0,T2,T4>::type T_123;
      typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

      T0 m_12;
      m_12 = sqrt(m2_12);

      // Invariant square mass of particles a,b,c together
      T_123 m_123, m2_123;
      m2_123 = m2_12 + m2_13 + m2_23 - a.m2 - b.m2 - c.m2;
      m_123 = sqrt(m2_123);
      
      // Form factor P -> R_1 d
      T_123 F_P = fct::blatt_weisskopf(this->l_1, this->P.r2, this->P.m2, 
				 m_123, this->d.m) / 
	fct::blatt_weisskopf(this->l_1, this->P.r2, this->P.m2,
			     this->R_1.m, this->d.m);

      // Form factor R_1 -> R_2 c
      T_123 F_R_1 = fct::blatt_weisskopf(this->l_2, this->R_1.r2, m2_123,
				 m_12, this->c.m) / 
	fct::blatt_weisskopf(this->l_2, this->R_1.r2, this->R_1.m2,
			     m_12, this->c.m);

      // Form factor R_2 -> a b
      T0 F_R_2 = fct::blatt_weisskopf(this->l_3, this->R_2.r2, m2_12,
				   this->a.m, this->b.m) / 
	fct::blatt_weisskopf(this->l_3, this->R_2.r2, this->R_2.m2,
			     this->a.m, this->b.m);

      // Dynamical (Breit-Wigner) form factor of the first resonance
      T_123 width_R_1 = fct::breit_wigner::relativistic_width(this->R_1.m,
        W_R_1, this->l_2, this->R_1.r, m2_123, m2_12, c.m2);
      std::vector<T_123> T_R_1 = fct::breit_wigner::value(this->R_1.m, m2_123,width_R_1);

      // Dynamical (Breit-Wigner) form factor of the 2nd resonance
      T_123 width_R_2 = fct::breit_wigner::relativistic_width(this->R_2.m,
        W_R_2, this->l_3, this->R_2.r, m2_12, a.m2, b.m2);
      std::vector<T_123> T_R_2 = fct::breit_wigner::value(this->R_2.m, m2_12,width_R_2);

      // Zemach tensors
      // Calculate transformed variables z2, cos2_theta for 
      // the decay D-> R_1 d -> R_2 c d
      T_123 p2_c = fct::breakup_momentum::p2(m2_123,
					     m_12, this->c.m);
      T_123 p2_d = fct::breakup_momentum::p2(this->P.m2, 
					      m_123, this->d.m);

      T_123 E_c = sqrt(this->c.m2 + p2_c);
      T_123 E_d = sqrt(this->d.m2 + p2_d);

      T_res p_c_dot_p_d = (-0.5) * (m2_34 - c.m2 - d.m2 - 2.0 * E_c * E_d);
      T_res cos2_theta = p_c_dot_p_d * p_c_dot_p_d / p2_c / p2_d;

      T_123 s = m2_123 + d.m2 + 2.0 * m_123 * E_d;
      T_123 z2 = p2_d / s;

      T_res Z_1 = fct::zemach(this->P.J, this->R_1.J, l_1, z2, cos2_theta);


      // Calculate transformed variables z2, cos2_theta for 
      // the decay R_1 -> R_2 c -> a b c
      // p2_c stays the same as above
      //p2_c = fct::breakup_momentum::p2(m2_123, m_12, this->c.m);
      T_123 p2_b = fct::breakup_momentum::p2(m2_12, this->a.m, this->b.m);

      // E_c stays the same as above
      // E_c = sqrt(this->c.m2 + p2_c);
      T_123 E_b = sqrt(this->b.m2 + p2_b);

      T_res p_b_dot_p_c = (-0.5) * (m2_23 - b.m2 - c.m2 - 2.0 * E_b * E_c);
      // The variables below are re-used.
      cos2_theta = p_b_dot_p_c * p_b_dot_p_c / p2_b / p2_c;

      s = m2_12 + c.m2 + 2.0 * m_12 * E_c;
      z2 = p2_c / s;

      T_res Z_2 = fct::zemach(this->R_1.J, this->R_2.J, l_2, z2, cos2_theta);

      // Combine the factors to the decay amplitude
      std::vector<T_res> A(2);
      A = complex::scalar::mult(F_P * F_R_1 * Z_1 * F_R_2  * Z_2,
				complex::scalar::mult(T_R_1,T_R_2));

      return A;
    }

  };

}

#endif
