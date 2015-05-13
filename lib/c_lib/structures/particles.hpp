#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__PARTICLES_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__PARTICLES_HPP

#include <meson_deca/lib/c_lib/structures/struct_particles.hpp>

namespace particles {

  // "Particles"
  particle pi(0.13957, 5., 0); // Charged pion
  particle k(0.49368, 5., 0);  // Charged kaon
  particle d(1.86960, 5., 0); // Charged D meson
  particle D0(1.86960, 5., 0); // Charged D meson

  // "Resonances"
  /* The description of resonances as particles needs certain 
   * clarification. Namely, if we describe resonances as particles,
   * why do we need the radius, and what about the width of the
   * resonance? Well, we need the radius, because in the four
   * body decay we may use resonance as a parent particle.
   * We do not define the widths here, because width is more
   * model-dependent (in a certain manner); therefore, it seems
   * more convenient to define width of a resonance while
   * instantiating a given decay channel (see 'resonances.hpp').
   */
  // Spin 0
  particle toy0_1000(1000, 5., 0); // Toy resonance
  particle toy0_1200(1200, 5., 0); // Toy resonance

  particle f0_600(0.800, 5., 0); // a.k.a sigma
  particle f0_980(0.980, 5., 0);
  particle f0_1370(1.350, 5., 0);
  particle f0_1500(1.507, 5., 0);

  // Spin 1
  particle a1(1.260, 5., 1); // Needs precision
  particle rho_770(0.770, 5., 1); // Needs precision

  // Spin 2
  particle f2_1270(1.2754, 5., 2);

}

#endif
