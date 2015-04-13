"""
fct
===

Provides
  1. Blatt-Weisskopf form factors blatt_weisskopf()
  2. Breit-Wigner dynamical function breit_wigner()
  3. Zemach tensors zemach()
  4. Flatte form factor flatte()


Contains

   blatt_weisskopf()
   breakup_momentum_complex_p()
   breakup_momentum_p2()
   breakup_momentum_r2()
   breit_wigner()
   flatte()
   relativistic_breit_wigner_width()
   resonance_bw()
   uniform_dalitz_points()
   valid()
   zemach()


How to use the documentation
----------------------------
Documentation is provided in form of docstrings. The docstrings examples
assume that fct has been imported without a prefix, i.e. start python from
the folder containing 'fct.py' and call

  >>> from fct import *

For additional information on a function func, type 'help(func)'.
"""

import numpy as np
import sys

sys.path.insert(1, "../phys") 
import particles

def blatt_weisskopf(J_R, r2_P, m2_ab, m_a, m_b = -1):
    """
    Blatt Weisskopf form factors.

    Implemented as in: arxiv:1406.6311v2, p.151, eq. (13.2.8).

    J_R   : resonance spin
    r2_P  : parent particle squared radius
    m2_ab : Dalitz plot variable (squared mass)
    m_a   : 1st daughter particle mass of m2_ab
    m_b   : 2nd daughter particle mass of m2_ab
    """
    if (J_R == 0 or J_R > 2):
        return 1;

    p2 = breakup_momentum_p2(m2_ab, m_a, m_b) * r2;
    if (J_R == 1):
        return np.sqrt(1.0/(1.0 + p2));
    elif (J_R == 2):
        return np.sqrt(1.0/(9.0 + 3.0 * p2 + p2 ** 2));


def breakup_momentum_complex_p(m2_R, m_a, m_b):
    """Complex breakup momentum."""
    p2 = breakup_momentum_p2(m2_R, m_a, m_b)
    if (p2>=0):
        return np.sqrt(p2)
    else:
        return 1j*np.sqrt(-p2)


def breakup_momentum_p2(m2_R, m_a, m_b):
    """
    Squared breakup momentum (R -> ab).

    'R' is the parent particle decaying into 'a' and 'b'. The value m_b < 0 
    is interpreted as the case of identical particles 'a' and 'b'.
    """
    if (m_b < 0):
        return m2_R / 4.0 - m_a ** 2;
    return (m2_R - (m_a + m_b) ** 2) * (m2_R - (m_a - m_b) ** 2) / m2_R / 4.0;


def breakup_momentum_r2(m2_N, m2_d, m_a, m_b = -1):
    """Ratio of squared breakup momenta (N -> ab) / (D -> ab)."""
    return breakup_momentum_p2(m2_N, m_a, m_b) / \
               breakup_momentum_p2(m2_d, m_a, m_b);



def breit_wigner(M_R, m2_ab, width_m2_ab):
    """
    Breit Wigner dynamical function.

    M_R     : resonance mass
    m2_ab   : Dalitz plot variable
    """
    return 1.0/(M_R ** 2 - m2_ab - 1.0j * M_R * width_m2_ab);


# TO DO: Eliminate the dependence on specific particles.
def flatte(M, m2, gpp, gkk):
    """
    Flatte resonance form factor.

    Implemented as in: arxiv:1406.6311v2, p. 150, eq. (13.2.5)

    M   : Resonance mass
    m2  : Dalitz plot variable
    gpp : phase-space factor * coupling constant^2 of the channel -> pi + pi 
    gkk : - // -                                   of the cahnnel -> K + K
    """
    pp = gpp * gpp * breakup_momentum_complex_p(m2, pwamath.m_pi, pwamath.m_pi)
    kk = gkk * gkk * breakup_momentum_complex_p(m2, pwamath.m_k, pwamath.m_k)
    return ((M*M - m2) - (2. / np.sqrt(m2)) * 1j * (pp + kk))


def relativistic_breit_wigner_width(M_R, W_R, J_R, r_R, m2_ab, m_a, m_b=-1):
    """
    Relativistic Breit Wigner resonance width.

    Implemented as in: arxiv:1406.6311v2, p.150, eq. (13.2.4).

    M_R   : resonance mass
    W_R   : resonance width
    J_R   : resonance spin
    r_R   : resonance radius
    m2_ab : Dalitz plot variable
    m_a   : 1st daughter mass
    m_b   : 2nd daughter mass
    """
    width = W_R *                                                        \
            M_R / np.sqrt(m2_ab)   *                                     \
            np.power(breakup_momentum_r2(m2_ab, M_R ** 2, m_a, m_b), 
                     J_R + 0.5) *                                          \
            blatt_weisskopf(J_R, r_R ** 2, m2_ab, m_a, m_b) ** 2 /          \
            blatt_weisskopf(J_R, r_R ** 2, M_R ** 2, m_a, m_b) ** 2
    return width


def resonance_bw(m2_ab, m2_bc, A_R = 1.0, 
                               J_R = 0.0, 
                               M_R = 1.0, 
                               W_R = 0.1, 
                               r_R = 1.0,
                               m2_p = particles.m2_d,
                               m2_a = particles.m2_pi,
                               m2_b = particles.m2_pi,
                               m2_c = particles.m2_pi,
                               r2_p = particles.r2_d):
    """
    Resonance amplitude; model-dependent with Breit-Wigner form factor.

    m2_ab : 1st Dalitz coordinate
    m2_bc : 2nd Dalitz coordinate
    A_R   : resonance amplitude (complex)
    J_R   : resonance spin
    M_R   : resonance mass
    W_R   : resonance width
    r_R   : resonance radius
    m2_p  : parent particle square mass
    m2_a  : 1st daughter square mass
    m2_b  : 2nd daughter square mass
    m2_c  : 3rd daughter square mass (bachelor particle)
    r2_p  : parent particle square radius

    Remember that the amplitude is not symmetrized here; for equal daughter
    particles, further symmetrization wrt. m2_ab, m2_bc is required.
    """

    if valid(m2_ab, m2_bc, m2_p, m2_a, m2_b, m2_c) == True:

        m_ab = np.sqrt(m2_ab)
        m_a = np.sqrt(m2_a)
        m_c = np.sqrt(m2_c)
        m_p = np.sqrt(m2_p) 
        if m2_b == (-1):
            m_b = 1
        else:
            m_b = np.sqrt(m2_b)

        # Form factor P -> Rc
        F_D = blatt_weisskopf(J_R, r2_p, m2_p, m_ab, m_c) /     \
              blatt_weisskopf(J_R, r2_p, m2_p, M_R,  m_c)

        # Form factor R -> ab
        F_R = blatt_weisskopf(J_R, r_R ** 2, m2_ab,    m_a, m_b) /     \
              blatt_weisskopf(J_R, r_R ** 2, M_R ** 2, m_a, m_b)

        width = relativistic_breit_wigner_width(M_R, W_R, J_R, r_R, 
                                                m2_ab, m_a, m_b)

        T = breit_wigner(M_R, m2_ab, width)
        Z = zemach(J_R, m2_ab, m2_bc, m_p, m_a, m_b, m_c)

        return A_R *  F_D * F_R * T * Z

    else:
        return 0


def uniform_dalitz_points(D, m2_p = particles.m2_d,
                             m2_a = particles.m2_pi,
                             m2_b = -1,
                             m2_c = -1):
    """
    Returns D points uniformely distributed over the Dalitz plot,
    assuming the decay D -> 3pi.
    """
    d = 0 # Counter for points
    m2_ab_list = np.zeros(D, dtype=float) # Point coordinates
    m2_bc_list = np.zeros(D, dtype=float)

    while d < D:
        m2_ab = np.random.uniform(0,3)
        m2_bc = np.random.uniform(0,3)
        if valid(m2_ab, m2_bc, m2_p, m2_a, m2_b, m2_c) == True:
            m2_ab_list[d] = m2_ab
            m2_bc_list[d] = m2_bc
            d += 1

    return [m2_ab_list, m2_bc_list]


def valid(m2_ab, m2_bc, m2_p, m2_a, m2_b=-1, m2_c=-1):
    """
    This function determines whether the point (m2_ab, m2_bc) is
    within the bound of the Dalitz plot.
    """
    if m2_b < 0:
        m2_b = m2_a
    if m2_c < 0:
        m2_c = m2_a

    #if (m2_ab < 0) or (m2_bc < 0):
    #    return False
    if (m2_ab < m2_a + m2_b + 2 * np.sqrt(m2_a * m2_b)) or \
       (m2_ab > m2_p + m2_c - 2 * np.sqrt(m2_p * m2_c)):
        return False

    E_b = (m2_ab - m2_a + m2_b)/2/np.sqrt(m2_ab)
    E_c = (m2_p - m2_ab - m2_c)/2/np.sqrt(m2_ab)
    P_b = np.sqrt(E_b ** 2 - m2_b)
    P_c = np.sqrt(E_c ** 2 - m2_c)

    if np.abs(m2_bc - m2_b - m2_c - 2 * E_b * E_c) <= 2 * P_b * P_c:
        return True

    return False


def zemach(J_R, m2_ab, m2_bc, m2_R, m2_a, m2_b = -1, m2_c = -1):
    """
    Zemach angular dependence function.

    Implemented as in: arxiv:1406.6311v2, p.151, eq. (13.2.6).

    By default, m2_b and m2_c are set to -1, which is interpreted as
    particles 'b' and 'c' being identical with the particle 'a'.

    J_R   : resonance spin
    m2_ab : 1st Dalitz plot coordinate
    m2_bc : 2nd Dalitz plot coordinate
    m2_R  : resonance squared mass
    m2_a  : 1st daughter squared mass
    m2_b  : 2nd daughter squared mass
    m2_c  : 3rd daughter squared mass
    """
    if (J_R > 2):
        return 0
    elif (J_R == 0):
        return 1

    # By default, particles 'b' and 'c' are identical to 'a'
    if (m2_b < 0):
        m2_b = m2_a;
    if (m2_c < 0):
        m2_c = m2_a;

    if (J_R == 1):
        return m2_R + m2_a + m2_b + m2_c - m2_ab - 2*m2_bc - \
               (m2_R - m2_c) * (m2_a - m2_b)/m2_ab;

    elif (J_R == 2):
        return (zemach(1, m2_ab, m2_bc, m2_R, m2_a, m2_b, m2_c) ** 2) -     \
               (m2_ab - 2 * m2_R - 2 * m2_c + (m2_R - m2_c) ** 2 / m2_ab) * \
               (m2_ab - 2 * m2_a - 2 * m2_b + (m2_a - m2_b) ** 2 / m2_ab)/3.0;
