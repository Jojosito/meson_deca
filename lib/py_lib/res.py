"""
Contains amplitudes of various resonances.
"""

import numpy as np
import sys

sys.path.insert(1, "../math")
import fct
reload(fct)
import particles

def res_playground(m2_ab, m2_bc, A_m, A_ph, J=0, M=1.0, W=0.1):
    """
    Playground Breit-Wigner resonance.

    m2_ab, m2_bc : Dalitz plot coordinates
    A_m, A_ph    : Resonance amplitude modulus and phase
    J, M, W      : Resonance spin, mass, and width, respectively
    """
    return fct.resonance_bw(m2_ab, m2_bc, 
                            A_m * np.exp(1j * A_ph / 180.0 * np.pi),
                            J, M, W)


def H2_playground(m2_ab, m2_bc, A_m_diff, A_ph_diff):
    """
    Toy probability of 2 Breit-Wigner resonances with small mass separation.

    The mass separation is preset to 0.2.

    A_m_diff is the difference in amplitude modulus of both resonances.
    A_ph_diff is the difference in phase of both resonances.
    """
    # Remember to symmetrize the amplitudes
    J = 0
    M1 = 1.0
    deltaM = 0.2
    M2 = M1 + deltaM
    W = 0.1
    return np.square(np.abs(res_playground(m2_ab, m2_bc, 1.0, 0,
                                           J, M1, W) +            \
           res_playground(m2_ab, m2_bc, 1.0 + A_m_diff, A_ph_diff,
                          J, M2, W) +                             \
           res_playground(m2_bc, m2_ab, 1.0, 0,
                          J, M1, W) +                             \
           res_playground(m2_bc, m2_ab, 1.0 + A_m_diff, A_ph_diff,
                          J, M2, W)))


#f_0_500.J = 0
#f_0_500.M = 1.0#0.475
#f_0_500.W = 0.1#0.550
#f_0_500.a = 1.0 #1.3 * np.exp(1j * (-21)/180 * np.pi)
#f_0_500.r = 1
#f_0_500.r2 = f_0_500.r ** 2;



#f_0_1370.J = 0
#f_0_1370.M = 1.2#1.350
#f_0_1370.W = 0.1#0.050
#f_0_1370.a = 1.0 * np.exp(1j * np.pi / 180 * 350)#1.3 * np.exp(1j * (-21)/180 * np.pi)
#f_0_1370.r = 1
#f_0_1370.r2 = f_0_1370.r ** 2;

# Resonance: f_0 (980)
#def f_0_980_sym(m2_ab, m2_bc):
#    return abs(f_0_980(m2_ab, m2_bc) + f_0_980(m2_bc, m2_ab)) ** 2


def f_0_980(m2_ab, m2_bc):

    if fct.valid(m2_ab, m2_bc, particles.m2_d, particles.m2_pi):
        
        m_ab = np.sqrt(m2_ab)
        F_D = fct.blatt_weisskopf(f_0_980.J, particles.r2_d, particles.m2_d, 
                                  m_ab, particles.m_pi)/ \
              fct.blatt_weisskopf(f_0_980.J, particles.r2_d, particles.m2_d,
                                 f_0_980.M, particles.m_pi)

        F_R = fct.blatt_weisskopf(f_0_980.J, f_0_980.r2, m2_ab, 
                                  particles.m_pi)/ \
              fct.blatt_weisskopf(f_0_980.J, f_0_980.r2, f_0_980.M ** 2,
                                  particles.m_pi)

        T = fct.flatte(f_0_980.M, m2_ab ,f_0_980.gpp, f_0_980.gkk)

        Zemach = fct.zemach(f_0_980.J, m2_ab, m2_bc, particles.m2_d, 
                            particles.m2_pi)

        return f_0_980.a * F_D * F_R * T * Zemach
    return 0

f_0_980.J = 0
f_0_980.M = 0.980
f_0_980.gpp = 0.329
f_0_980.gkk = 2. * 0.329
f_0_980.a = 1.4 * np.exp(1j * (12) / 180 * np.pi)
f_0_980.r = 1.
f_0_980.r2 = f_0_980.r ** 2;



# Resonance: f_2 (1270)
#def f_2_1270_sym(m2_ab, m2_bc):
#    return abs(f_2_1270(m2_ab, m2_bc) + f_2_1270(m2_bc, m2_ab)) ** 2


#def f_2_1270(m2_ab, m2_bc):

#    if Valid(m2_ab, m2_bc, particles.m2_d, particles.m2_pi):
        
#        m_ab = np.sqrt(m2_ab)
#        F_D = F(f_2_1270.J, particles.r2_d, particles.m2_d, m_ab, particles.m_pi)/ \
#             F(f_2_1270.J, particles.r2_d, particles.m2_d, f_2_1270.M, particles.m_pi)

#       F_R = F(f_2_1270.J, f_2_1270.r2, m2_ab, particles.m_pi)/ \
#             F(f_2_1270.J, f_2_1270.r2, f_2_1270.M ** 2, particles.m_pi)
#
#        width = rBW_width(f_2_1270.M, f_2_1270.W, f_2_1270.J, f_2_1270.r, m2_ab)
#        T = BW(f_2_1270.M, m2_ab, width)

#        Zemach = W(f_2_1270.J, m2_ab, m2_bc, particles.m2_d, particles.m2_pi)

#        return  F_R#F(f_2_1270.J, particles.r2_d, particles.m2_d, m_ab, particles.m_pi)#F_R#f_2_1270.a * Zemach# * F_D * F_R = *math.r2_d, particles.m2_d, m_ab, particles.m_pi) T * Zemach 
#    return 0


#f_2_1270.J = 2
#f_2_1270.M = 1.2754
#f_2_1270.W = 0.1852
#f_2_1270.a = 2.1 * np.exp(1j * (-123)/180 * np.pi)
#f_2_1270.r = 1
#f_2_1270.r2 = f_0_1370.r ** 2;


 
# Resonance: f_0 (1500)
#def f_0_1500_sym(m2_ab, m2_bc):
#    return abs(f_0_1500(m2_ab, m2_bc) + f_0_1500(m2_bc, m2_ab)) ** 2


#def f_0_1500_s_ca(m2_ab, m2_bc):
#    res = f_0_1500(m2_ab, m2_bc) + f_0_1500(m2_bc, m2_ab)
#    return res#np.asarray([np.abs(res), np.angle(res)])


#def f_0_1500(m2_ab, m2_bc):

#    if Valid(m2_ab, m2_bc, particles.m2_d, particles.m2_pi):
        
#        m_ab = np.sqrt(m2_ab)
#        F_D = F(f_0_1500.J, particles.r2_d, particles.m2_d, m_ab, particles.m_pi)/ \
#             F(f_0_1500.J, particles.r2_d, particles.m2_d, f_0_1500.M, particles.m_pi)
#
#        F_R = F(f_0_1500.J, f_0_1500.r2, m2_ab, particles.m_pi)/ \
#              F(f_0_1500.J, f_0_1500.r2, f_0_1500.M ** 2, particles.m_pi)
#
#        width = rBW_width(f_0_1500.M, f_0_1500.W, f_0_1500.J, f_0_1500.r, m2_ab)
#        T = BW(f_0_1500.M, m2_ab, width)

#        Zemach = W(f_0_1500.J, m2_ab, m2_bc, particles.m2_d, particles.m2_pi)

#        return f_0_1500.a * F_D * F_R * T * Zemach
#    return 0


#f_0_1500.J = 0
#f_0_1500.M = 1.2#1.507
#f_0_1500.W = 0.2#0.109
#f_0_1500.a = 1.0 #1.1 * np.exp(1j * (-44)/180 * np.pi)
#f_0_1500.r = 1
#f_0_1500.r2 = f_0_1500.r ** 2
