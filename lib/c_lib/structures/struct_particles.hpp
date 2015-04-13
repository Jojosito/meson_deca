#ifndef MESON_DECA__LIB__C_LIB__STRUCTURES__STRUCT_PARTICLES_HPP
#define MESON_DECA__LIB__C_LIB__STRUCTURES__STRUCT_PARTICLES_HPP

struct particle 
{
  const double m; // Mass in GeV
  const double m2; // Squared mass
  const double r; // Radius in GeV**-1
  const double r2; // Squared radius

  particle(double _m, double _r) : 
    m(_m), m2(_m*_m), r(_r), r2(_r*_r) {};
};

#endif
