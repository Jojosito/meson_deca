#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__BASE_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__THREE_BODY__BASE_HPP

#include <meson_deca/lib/c_lib/structures/struct_particles.hpp>

namespace resonances {

  // Base struct for model-dependent 3-body-decay resonances
  struct resonance_base_3
  {

    const particle P; // Parent particle
    const particle a; // Final state particles
    const particle b; 
    const particle c;

    resonance_base_3(particle _P, particle _a, particle _b, particle _c) :
      P(_P), a(_a), b(_b), c(_c) {};
  };

}
#endif
