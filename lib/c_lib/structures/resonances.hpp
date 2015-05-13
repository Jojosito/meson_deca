#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__RESONANCES_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__RESONANCES_HPP

// 3 body decay resonances
#include <meson_deca/lib/c_lib/structures/three_body/bw.hpp>
#include <meson_deca/lib/c_lib/structures/three_body/flat.hpp>
#include <meson_deca/lib/c_lib/structures/three_body/flatte.hpp>

// 4 body decay resonances
#include <meson_deca/lib/c_lib/structures/four_body/D_R1d_R2cd_abcd.hpp>
#include <meson_deca/lib/c_lib/structures/three_body/flat.hpp>

namespace resonances {

  // 3-body decay, Spin 0
  resonances::flat_3 flat_D3pi(particles::d, particles::pi, particles::pi, 
			       particles::pi);

  resonances::breit_wigner toy0_1000(particles::d, particles::pi,
				     particles::pi, particles::pi,
				     particles::toy0_1000, 0.1);

  resonances::breit_wigner toy0_1200(particles::d, particles::pi,
				     particles::pi, particles::pi,
				     particles::toy0_1200, 0.1);

  resonances::flatte toy0_flatte(particles::d, particles::pi,
				 particles::pi, particles::pi,
				 particles::toy0_1000, 0.329, 2*0.329);

  resonances::flatte f0_980(particles::d, particles::pi,
			    particles::pi, particles::pi,
			    particles::f0_980, 0.329, 2*0.329);

  resonances::breit_wigner f0_600(particles::d, particles::pi,
				  particles::pi, particles::pi,
				  particles::f0_600, 0.800);

  resonances::breit_wigner f0_1370(particles::d, particles::pi,
				   particles::pi, particles::pi,
				   particles::f0_1370, 0.350);

  resonances::breit_wigner f0_1500(particles::d, particles::pi,
				   particles::pi, particles::pi,
				   particles::f0_1370, 0.109);


  // 3-body decay, Spin 1
  resonances::breit_wigner rho_770(particles::d, particles::pi,
				   particles::pi, particles::pi,
				   particles::rho_770, 0.1491);

  // 3-body decay, Spin 2
  resonances::breit_wigner f2_1270(particles::d, particles::pi,
				   particles::pi, particles::pi,
				   particles::f2_1270, 0.1852);


  // 4-body resonances

  // Flat
  resonances::flat_4  flat_4_res(particles::D0, particles::pi,
				 particles::pi, particles::pi,
				 particles::pi);

  // Consecutive decays
  // To be adjusted: width of a1
  resonances::P_R1d_R2cd_abcd D_a_rho_S_wave(particles::D0, 
					     particles::pi, particles::pi, 
					     particles::pi, particles::pi,
					     1, 0, 1,
					     particles::a1, 
					     particles::rho_770, 
					     0.1, 0.1491); 

}

#endif
