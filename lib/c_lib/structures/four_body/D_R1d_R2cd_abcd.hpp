#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__D_R1D_R2cd_abcd_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_BODY__D_R1D_R2cd_abcd_HPP

#include <cmath> // sqrt
#include <math.h> // isnan

#include <meson_deca/lib/c_lib/fct.hpp> // Breit-Wigner, Blatt-Weisskopf, etc.
#include <meson_deca/lib/c_lib/complex.hpp> // Complex numbers

#include <meson_deca/lib/c_lib/structures/four_body/base.hpp> // base class

#include <assert.h>

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
    P_R1d_R2cd_abcd(particle _P, particle _a, particle _b, 
		    particle _c, particle _d,
		    int _l_1, int _l_2, int _l_3,  
		    particle _R_1, particle _R_2,
		    double _W_R_1, double _W_R_2) : 
      resonance_base_4(_P,_a, _b, _c, _d), 
      l_1(_l_1), l_2(_l_2), l_3(_l_3),
      R_1(_R_1), R_2(_R_2), W_R_1(_W_R_1), W_R_2(_W_R_2) {};


    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABCD (not symmetrized)
    template <typename T0, typename T1, typename T2, typename T3, typename T4>
    std::vector<typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type >
    // m2_12 is the invariant square mass of particles a and b.
    // Analogously, m2_34 is i.sq.m. of c and d, m2_23 - of b and c, etc.
    value(int debug, const T0& m2_12, const T1& m2_14, const T2& m2_23,
	  const T3& m2_34, const T4& m2_13) {

      typedef typename boost::math::tools::promote_args<T0,T2,T4>::type T_123;
      typedef typename boost::math::tools::promote_args<T0,T1,T2,T3,T4>::type T_res;

      std::vector<T_res> A(2, 0.0);

      // Invariant square mass of particles a,b,c together
      T_123 m2_123;
      m2_123 = m2_12 + m2_13 + m2_23 - a.m2 - b.m2 - c.m2;
      
      // Check whether we are in the physically relevant phase space region
      if ( ! fct::valid_5d(m2_12, m2_14, m2_23, m2_34, m2_13,
          this->P, this->a, this->b, this->c, this->d) == true)
        return A; // 0

      T0 m_12;
      m_12 = sqrt(m2_12);
      T_123 m_123;
      m_123 = sqrt(m2_123);

      // Form factor P -> R_1 d
      T_123 F_P = fct::blatt_weisskopf(this->l_1, this->P.r2, this->P.m2, 
          m_123, this->d.m) /
          fct::blatt_weisskopf(this->l_1, this->P.r2, this->P.m2,
              this->R_1.m, this->d.m);

      // Form factor R_1 -> R_2 c
      T_123 F_R_1 = fct::blatt_weisskopf(this->l_2, this->R_1.r2, m2_123,
          this->R_2.m, this->c.m) /
          // POSSIBLY m_12 instead of R_2.m above and below
          fct::blatt_weisskopf(this->l_2, this->R_1.r2, this->R_1.m2,
              this->R_2.m, this->c.m);

      // Form factor R_2 -> a b
      T0 F_R_2 = fct::blatt_weisskopf(this->l_3, this->R_2.r2, m2_12,
          this->a.m, this->b.m) /
          fct::blatt_weisskopf(this->l_3, this->R_2.r2, this->R_2.m2,
              this->a.m, this->b.m);

      // Dynamical (Breit-Wigner) form factor of the first resonance
      T_123 width_R_1 = fct::breit_wigner::relativistic_width(this->R_1.m, W_R_1, 
							      this->l_2, this->R_1.r, 
							      m2_123, this->R_2.m2, 
							      c.m2);
      std::vector<T_123> T_R_1 = fct::breit_wigner::value(this->R_1.m, 
							  m2_123,width_R_1);

      // Dynamical (Breit-Wigner) form factor of the 2nd resonance
      T_123 width_R_2 = fct::breit_wigner::relativistic_width(this->R_2.m, W_R_2, 
							      this->l_3, this->R_2.r, 
							      m2_12, a.m2, b.m2);
      std::vector<T_123> T_R_2 = fct::breit_wigner::value(this->R_2.m, 
							  m2_12,width_R_2);

      // Zemach tensors
      // Calculate transformed variables z2, cos2_theta for 
      // the decay D-> R_1 d -> R_2 c d
      // in the rest frame of R_1
      T_123 p2_c = fct::breakup_momentum::p2(m2_123, m2_12, this->c.m);

      T_123 E_c = sqrt(this->c.m2 + p2_c);
      
      T_123 p2_d_rest_frame_of_P = fct::breakup_momentum::p2(this->P.m2, 
							     m2_123, this->d.m);
      T_123 E_d_rest_frame_of_P = sqrt(this->d.m2 + p2_d_rest_frame_of_P);
      // R_1 and d have the same abs. momenta |p2| in rest frame of P
      T_123 E_R_1_rest_frame_of_P = sqrt(m2_123 + p2_d_rest_frame_of_P);

      // Compute p2_d in the rest frame of R_1 by
      // performing a Lorents boost in the direction -v_R_1
      T_123 v_R_1 = sqrt(p2_d_rest_frame_of_P) / E_R_1_rest_frame_of_P;
      T_123 gamma = 1.0 / sqrt(1.0 - v_R_1 * v_R_1);

      T_123 p_d = gamma * ( sqrt(p2_d_rest_frame_of_P) + E_d_rest_frame_of_P * v_R_1 );
      T_123 p2_d = p_d*p_d;
      T_123 E_d = sqrt(this->d.m2 + p2_d);

      T_res p_c_dot_p_d = (-0.5) * (m2_34 - c.m2 - d.m2 - 2.0 * E_c * E_d);
      T_res cos2_theta = p_c_dot_p_d * p_c_dot_p_d / p2_c / p2_d;

      T_123 s = m2_123 + d.m2 + 2.0 * m_123 * E_d;
      T_123 z2 = p2_d / s;

      /*std::cout << "m2_34 \t\t" << "p2_c \t\t" << "E_c \t\t" << "p2_d_rest_frame_of_P \t"
      	<< "E_R_1_rest_frame_of_P  \t" << "v_R_1 \t\t"
      	<< "gamma \t" << "p2_d \t\t" << "E_d \t\t" << "p_c_dot_p_d \t"
      	<< "cos2_theta \t" << "s \t" << "z2 \n";

      std::cout << m2_34 << "\t" << p2_c << " \t" << E_c <<  " \t" << p2_d_rest_frame_of_P <<  " \t\t"
      		<< E_R_1_rest_frame_of_P <<  " \t\t" << v_R_1 <<  " \t"
      		<< gamma <<  " \t" << p2_d << " \t" <<  E_d << " \t" <<  p_c_dot_p_d
      		<<  " \t" << cos2_theta  <<  " \t" << s <<  " \t" << z2 << "\n";
*/
      //assert(cos2_theta >= 0. && cos2_theta <= 1.);
      //assert(z2 >= 0. && z2 <= 1.);

      T_res Z_1(0);
      if (p2_c >= 0 && p2_d_rest_frame_of_P >= 0)
        Z_1 = fct::zemach(this->P.J, this->R_1.J, l_1, z2, cos2_theta);

      //assert(Z_1 >= 0. && Z_1 <= 2.);

      if (std::isnan(Z_1)==true) {
        std::cout << "p2_d_rfo_P : " << p2_d_rest_frame_of_P << "\n";
        std::cout << "E_c : " << E_c << "\n";
        std::cout << "E_R_1_rfo_P : " << E_R_1_rest_frame_of_P << "\n";
        std::cout << "v_R_1 : " << v_R_1 << "\n";
        std::cout << "gamma : " << gamma << "\n";
        std::cout << "p2_d : " << p2_d << "\n";
        std::cout << "E_d : " << E_d << "\n";
        std::cout << "s : " << s << "\n";
        std::cout << "Z_1 :" << Z_1 << "\n";
      }

      // Todo: debug cos2_theta and z2

      // Calculate transformed variables z2, cos2_theta for 
      // the decay R_1 -> R_2 c -> a b c
      // in the rest frame of R_2

      // p2_c above is calculated in the rest frame of R_1.
      // We want it in the rest frame of R_2, so we perform a 
      // Lorentz boost again.
      T_123 E_R_2 = sqrt(m2_12 + p2_c);
      T_123 v_R_2 = p2_c / E_R_2;
      T_123 gamma_2 = 1.0 / sqrt(1.0 - v_R_2 * v_R_2);
      T_123 p2_c_rest_frame_of_R_2 = gamma_2 * (p2_c + E_R_2 * v_R_2);
      T_123 E_c_rest_frame_of_R_2 = sqrt(this->c.m2 + p2_c_rest_frame_of_R_2);
      
      T_123 p2_b = fct::breakup_momentum::p2(m2_12, this->a.m, this->b.m);
      T_123 E_b = sqrt(this->b.m2 + p2_b);
      
      T_res p_b_dot_p_c = (-0.5) * (m2_23 - b.m2 - c.m2 
      			    - 2.0 * E_b * E_c_rest_frame_of_R_2);
      // The variables below are re-used.
      T_123 cos2_theta_2 = p_b_dot_p_c * p_b_dot_p_c / p2_b / p2_c_rest_frame_of_R_2;

      T_123 s_2 = m2_12 + c.m2 + 2.0 * m_12 * E_c_rest_frame_of_R_2;
      T_123 z2_2 = p2_c / s_2;

      // Todo: debug cos2_theta_2 and z2_2

      T_res Z_2(0);
      if (p2_b >= 0 && p2_c_rest_frame_of_R_2 >= 0)
        Z_2 = fct::zemach(this->R_1.J, this->R_2.J, l_2, z2_2, cos2_theta_2);

      // Combine the factors to the decay amplitude
      A = complex::scalar::mult(F_P * F_R_1 * Z_1 * F_R_2  * Z_2,
				complex::scalar::mult(T_R_1,T_R_2));

      A = complex::scalar::complex(1.0, 0.0);
      switch (debug) {
      case 1 : return complex::scalar::complex(F_P, 0.0);
      case 2 : //return complex::scalar::complex(F_R_1, 0.0);
        return complex::scalar::complex(fct::valid_5d(m2_12, m2_14, m2_23, m2_34, m2_13,
                this->P,
                this->a, this->b, this->c, this->d), 0.);
      case 3 : return complex::scalar::complex(F_R_2, 0.0);
      case 4 : return complex::scalar::complex(Z_1, 0.0);
      case 5 : return complex::scalar::complex(Z_2, 0.0);
      case 6 : return T_R_1;
      case 7 : return T_R_2;
      default : return A;
      }
      //} end of if(valid)

      return A;
    }

  };

}

#endif
