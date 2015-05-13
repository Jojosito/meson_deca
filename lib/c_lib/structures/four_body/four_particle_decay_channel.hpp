#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_PARTICLE_DECAY_CHANNEL_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__FOUR_PARTICLE_DECAY_CHANNEL_HPP

#include <meson_deca/lib/c_lib/structures/struct_resonances.hpp>
#include <meson_deca/lib/c_lib/structures/struct_four_particle_decay_channel.hpp>
#include <meson_deca/lib/c_lib/structures/particles.hpp>

namespace resonances {

  resonances::P__R1d__R2cd__abcd__bw_decay 
  D_a_rho_S_wave(1, 0, 1, particles::D0, 
		 particles::a1, particles::rho0, 
		 particles::pi, particles::pi, 
		 particles::pi, particles::pi, 0.1, 0.1);

  resonances::flat_5d
  flat_5d_res(0.0);
}


#endif
