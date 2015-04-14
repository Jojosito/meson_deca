#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__STRUCT_RESONANCES_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__STRUCT_RESONANCES_HPP

#include <cmath> // sqrt

#include <meson_deca/lib/c_lib/fct.hpp>
#include <meson_deca/lib/c_lib/complex.hpp>

namespace resonances {

  struct reson_base // Base struct for model-dependent resonances
  {
    const int J; // Spin
    const double M; // Mass
    const double M2; // Squared mass
    const double R; // Radius
    const double R2; // Squared radius;

    reson_base(int _J, double _M, double _R) :
      J(_J), M(_M), M2(_M*_M), R(_R), R2(_R*_R) {};
  };


  struct breit_wigner : resonances::reson_base // Breit-Wigner resonance
  {
    const double W; // Width

    breit_wigner(int _J, double _M, double _W, double _R) :
      reson_base(_J, _M, _R), W(_W) {};


    // Evaluates the resonance at the given point in the Dalitz plot
    // for the decay P -> ABC (not symmetrized)
    template <typename T>
    inline
    std::vector<T>
    value(const T& m2_ab, const T& m2_bc,
          const particle &p, const particle &a, 
          const particle &b, const particle &c) 
    {

      if (fct::valid(m2_ab, m2_bc, p, a, b, c) == true) {

	T m_ab = sqrt(m2_ab);

        // Form factor P -> Rc
        T F_P = fct::blatt_weisskopf(this->J, p.r2, p.m2, m_ab,   c.m) /
                fct::blatt_weisskopf(this->J, p.r2, p.m2, this->M, c.m);

        // Form factor R -> ab
        T F_R = fct::blatt_weisskopf(this->J, this->R2, m2_ab, a.m, b.m)/
                fct::blatt_weisskopf(this->J, this->R2, this->M2,a.m, b.m);

        T width = fct::breit_wigner::relativistic_width(this->M, this->W,
	  this->J, this->R, m2_ab, a.m, b.m);

	std::vector<T> T_R = fct::breit_wigner::value(this->M, m2_ab, width);
        T Z = fct::zemach(this->J, m2_ab, m2_bc, p.m, a, b, c);

	std::vector<T> res(2);
        res[0] = F_P * F_R * T_R[0] * Z;
        res[1] = F_P * F_R * T_R[1] * Z;  

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
    value_sym(const T& m2_ab, const T& m2_bc,
              const particle &p, const particle &a, 
	      const particle &b, const particle &c) {

      return complex::scalar::add(this->value(m2_ab, m2_bc, p, a, b, c),
				  this->value(m2_bc, m2_ab, p, a, b, c));
    }

  };


  struct flatte : resonances::reson_base // Flatte Resonance
  {
    const double G_pp;
    const double G_kk;

    flatte(int _J, double _M, double _G_pp, double _G_kk, double _R) :
      reson_base(_J, _M, _R), G_pp(_G_pp), G_kk(_G_kk) {};

    template <typename T>
    std::vector<T>
    value(const T& m2_ab, const T& m2_bc,
          const particle &p, const particle &a, 
          const particle &b, const particle &c) 
    {

      if (fct::valid(m2_ab, m2_bc, p, a, b, c) == true) {

	T m_ab = sqrt(m2_ab);

        // Form factor P -> Rc
        T F_P = fct::blatt_weisskopf(this->J, p.r2, p.m2, m_ab,   c.m) /
	  fct::blatt_weisskopf(this->J, p.r2, p.m2, this->M, c.m);

        // Form factor R -> ab
        T F_R = fct::blatt_weisskopf(this->J, this->R2, m2_ab, a.m, b.m)/
	  fct::blatt_weisskopf(this->J, this->R2, this->M2, a.m, b.m);

	std::vector<T> T_R = fct::flatte::value(this->M, m2_ab,
						this->G_pp, this->G_kk);
        T Z = fct::zemach(this->J, m2_ab, m2_bc, p.m, a, b, c);

	std::vector<T> res(2);
	res[0] = F_P * F_R * T_R[0] * Z;
	res[1] = F_P * F_R * T_R[1] * Z;
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
    value_sym(const T& m2_ab, const T& m2_bc,
              const particle &p, const particle &a, 
	      const particle &b, const particle &c) {

      return complex::scalar::add(this->value(m2_ab, m2_bc, p, a, b, c),
				  this->value(m2_bc, m2_ab, p, a, b, c));
    }

  };
}

#endif
