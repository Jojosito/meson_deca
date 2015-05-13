#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__FLATTE_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__FLATTE_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/fct.hpp>
#include <meson_deca/lib/c_lib/complex.hpp>
#include <meson_deca/lib/c_lib/structures/three_body/base.hpp>

namespace resonances {

  // Flatte Resonance for the 3-body decay
  struct flatte : resonances::resonance_base_3
  {
    // Flatte has same properties as a particle + 2 widths
    const particle R;
    const double G_pp;
    const double G_kk;

    flatte(particle _P, particle _a, particle _b, particle _c,
	   particle _R, double _G_pp, double _G_kk) :
      resonance_base_3(_P, _a, _b, _c), R(_R), G_pp(_G_pp), G_kk(_G_kk) {};

    // Returns the amplitude of the decay P->abc via Flatte resonance.
    template <typename T>
    std::vector<T>
    value(const T& m2_ab, const T& m2_bc) 
    {
      if (fct::valid(m2_ab, m2_bc, 
		     this->P, this->a, this->b, this->c) == true) {

	T m_ab = sqrt(m2_ab);

        // Form factor P -> Rc
        T F_P = fct::blatt_weisskopf(this->R.J, this->P.r2, 
				     this->P.m2, m_ab, this->c.m) /
	  fct::blatt_weisskopf(this->R.J, this->P.r2,
			       this->P.m2, this->R.m, this->c.m);

        // Form factor R -> ab
        T F_R = fct::blatt_weisskopf(this->R.J, this->R.r2, 
				     m2_ab, this->a.m, this->b.m)/
	  fct::blatt_weisskopf(this->R.J, this->R.r2,
			       this->R.m2, this->a.m, this->b.m);

	std::vector<T> T_R = fct::flatte::value(this->R.m, m2_ab,
						this->G_pp, this->G_kk);
        T Z = fct::zemach(this->R.J, m2_ab, m2_bc, 
			  this->P.m, this->a, this->b, this->c);

	std::vector<T> res(2);
	res = complex::scalar::mult(F_P * F_R * Z, T_R);
        return res;
      }
      else {
	std::vector<T> res(2, 0.0);
        return res;
      }
    }


    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABC (symmetrized, i.e. A==C)
    template <typename T>
    inline
    std::vector<T>
    value_sym(const T& m2_ab, const T& m2_bc) {

      return complex::scalar::add(this->value(m2_ab, m2_bc),
				  this->value(m2_bc, m2_ab));
    }

  };
}

#endif
