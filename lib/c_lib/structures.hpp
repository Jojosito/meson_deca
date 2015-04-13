#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES_HPP

/*
 * structures.hpp
 *
 * Some elements of PWA could be easier stored in structures,
 * not in functions. Basically, your PWA works as follows:
 * you have a bunch of parameters and variables that are
 * passed from STAN, but you also have some constants known
 * from the previous analyses or PDG. Depending on your model,
 * it is convenient to 'glue' all the form-factors together
 * (declaring some new structure, task #1) and insert all the
 * constants (initializing this new structure, task #2).
 * 
 * The files 'struct_*' do the task #1: they contain the information
 * about the structures - for example, about the resonances,
 * which are glued from Blatt-Weisskopf, Breit-Wigner, or other 
 * functions. 
 * 
 * The other files declare these structures. For example, 
 * 'particles.hpp' initializes pion and kaon particles with
 * members that describe their mass, radius, etc.
 */
#include <meson_deca/lib/c_lib/structures/struct_particles.hpp>
#include <meson_deca/lib/c_lib/structures/particles.hpp>
#include <meson_deca/lib/c_lib/structures/struct_resonances.hpp>
#include <meson_deca/lib/c_lib/structures/resonances.hpp>

#endif