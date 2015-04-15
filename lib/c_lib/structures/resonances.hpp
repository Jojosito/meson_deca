#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__RESONANCES_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__RESONANCES_HPP

#include <meson_deca/lib/c_lib/structures/struct_resonances.hpp>
// classes resonances::breit_wigner, resonances::flatte

namespace resonances {

  const double math_pi = 3.1415926535;
  const double deg =  math_pi / 180.;

  //resonance.non_resonant d_3pi_flat(1.36, 150.1*deg);

  resonances::flat flat_res(0.0);

  resonances::breit_wigner toy0_1000(0, 1.,
                                   0.1,
                                   1.0);

  resonances::breit_wigner toy0_1200(0, 1.2,
                                   0.1,
                                   1.0);

  resonances::flatte toy0_flatte(0, 1.,
				 0.329, 2*0.329,
				 1.1);

  resonances::flatte f0_980(0, 0.980, // Spin, mass
                          0.329, 2*0.329, // Width
                          //1.4, 12.*deg,  // Amplitude
                          1.);  // Radius

  resonances::breit_wigner f0_600(0, 0.800, // Spin, mass
                                0.800,  // Width
                                //3.7, -3.*deg, // Amplitude
                                1.); // Radius

  resonances::breit_wigner f0_1370(0, 1.350,
                                 0.350,
                                 //1.3, -21.*deg,
                                 1.0);

  resonances::breit_wigner f0_1500(0, 1.507,
                                 0.109,
                                 //1.1, -44.*deg,
                                 1.0);

  resonances::breit_wigner rho_770(1, 0.770,
                                 0.1491,
				 //1.0, 0.*deg,
				  1.0);

  resonances::breit_wigner f2_1270(2, 1.2754,
                                 0.1852,
			         //2.1, -123.*deg,
			         1.0);

}

#endif
