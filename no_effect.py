#!/usr/bin/env python3

"""Calculation of the percentage effect on gamma counts for a given
temperature.

Calculates gamma radiation anisotropy for 54Mn and 60Co nuclear orientation
thermometers.

Args:
    T: Temperature in Kelvin
    
Returns:
    The drop in counts in percent from the warm value
     
Raises:
    Some errors/exceptions?
    
Notes:
    Thermometer material selection is contained in module 'material.py'.
    Usually 60Co in hcp Co at Oxford Instruments at the time of writing.

    The detector angle is set in module 'detector_angle.py'. Usually set
    to zero at Oxford Instruments.

    The source-detector distance and detector face diameter are set in
    module 'solid_angle.py'.

    Some physical constants are set up in module 'constants.py'

    The numbers in square brackets e.g. [2] refer to references contained
    in file 'references.md'.

Current status at 15/3/17
Still setting up and getting under git control.
changing modules to all lower case.
getting imports correct.
sorting docstrings
don't forget to git fetch!
changing this module just to calc percentrage effect and move all other temp calcs to another module.
probably therefore not if name = main on this one?!

"""

__author__ = "Simon C. J. Kingsley"
__email__ = "simon.kingsley@oxinst.com"
__email__ = "scjk@btinternet.com"

from math import exp as exp
from math import factorial as fact
from math import sqrt as sqrt

import constants
import material
import detector_angle
import solid_angle

def em(gn, Bhf, J):
    """
    Calculates the Zeeman splitting energy per nuclear-spin
    angular-momentum sub-state. The absolute value is taken because
    it's just a delta. This means positive temperatures are returned
    even with negative hyperfine fields.
    I'm not convinced about the reasons for this... but it's consistent
    with the data in Table I reference [2].
    gn = nuclear magnetic dipole moment (nuclear g-factor),
         dimensionless.
    Bhf = hyperfine field, T.
    J = nuclear-spin angular-momentum, dimensionless.

    >>> em(1, 1, 1)
    5.05084e-27
    >>> em(0.2, -1, 1)
    1.010168e-27
    >>> em(3.799, -21.92, 5)/1.38066e-23
    0.006092796984445121
    >>>

    """
    # Pull in fundamental constants. No shit.
    mun = constants.NUCLEAR_MAGNETON

    # Calculate energy in J.
    result = abs((gn * mun * Bhf) / J)

    return result


def boltzmann(m, J, deltae, T):
    """
    Calculates the relative population of the mth nuclear spin
    sub-state using the Boltzmann distribution function.
    Equation 6a [2].
    m         = nuclear spin sub-state. Runs from -I to I.
    I         = nuclear-spin angular-momentum, dimensionless.
    deltae    = Zeeman energy splitting per nuclear spin sub-state.
    T         = temperature

    >>> boltzmann(1, 5, 1.38e-23, 1)
    0.0015709598596227391
    >>> boltzmann(1, 3, 1.38e-23, 1)
    0.01160724312949314
    >>>

    """
    # Pull in fundamental constant.
    kB = constants.BOLTZMANNS_CONSTANT

    # Calculate top part of Bolztmann population.
    population = exp((-1 * m * deltae) / (kB * T))

    # Calculate bottom part of Boltzmann population - the partition
    # function.
    partition_sum = 0.0
    for n in range(-J, J + 1):
        partition_sum += exp((-1 * n * deltae) / (kB * T))
    partition = partition_sum

    # Calculate the relative population.
    result = population / partition

    return result


def f2(J, T):
    """
    Calculates the second nuclear spin distribution moment.
    See equation 7a [2].
    J = nuclear spin
    T = temperature

    >>> f2(5, 0.02)
    0.11417244178609942
    >>> f2(3, 0.01)
    0.31034212203379946
    >>>

    """
    # Get material moment and hyperfine field
    gn = material.moment
    Bhf = material.field

    # Get Zeeman energy splitting per spin sub-state.
    deltae = em(gn, Bhf, J)

    # Calculate f2. 'sigma' is the sum term in Eq 7a [2].
    # 'boltzmann()' is the 'p(m)' term in [2] Eq 7a.
    sigma = 0.0
    for m in range(-J, J + 1):
        bol = boltzmann(m, J, deltae, T)
        sigma += (m ** 2) * bol

    f2_result = (1 / (J ** 2)) * (sigma - (J * (J + 1) / 3))

    return f2_result


def f4(J, T):
    """
    Calculates the fourth nuclear spin distribution moment.
    See equation 7b [2].
    J = nuclear spin
    T = temperature

    >>> f4(5, 0.02)
    0.001575043210768672
    >>> f4(3, 0.01)
    0.011632331479432023
    >>>

    """
    # Get material moment and hyperfine field
    gn = material.moment
    Bhf = material.field

    # Get Zeeman energy splitting per spin sub-state.
    deltae = em(gn, Bhf, J)

    # Calculate f4. 'sigma1{2}' is the first{second} sum term in
    # [2] Eq 7b.
    # 'boltzmann()' is the 'p(m)' term in [2] Eq 7b.
    sigma1 = 0.0
    sigma2 = 0.0

    for m in range(-J, J + 1):
        bol = boltzmann(m, J, deltae, T)
        sigma1 += (m ** 4) * bol
        sigma2 += (m ** 2) * bol

    a1 = ((6 * J ** 2) + (6 * J) - 5) / 7
    a2 = 3 * J * (J - 1) * (J + 1) * (J + 2) / 35

    f4_result = (1 / (J ** 4)) * (sigma1 - (a1 * sigma2) + a2)

    return f4_result


def B2(J, T):
    """
    Calculates the B2 term.
    J = nuclear spin
    T = temperature
    The following self tests only work for cobalt-60 with the specfic
    nuclear magnetic dipole moment, hyperfine field etc. etc.

    >>> B2(5, 0.00121856)
    1.691513076072361
    >>> B2(5, 0.0030464)
    1.5444638348112194
    >>> B2(5, 0.0060928)
    1.1820069590699673
    >>> B2(5, 0.060928)
    0.04290865998087582
    >>>

    """

    roottop = (2 * J + 1) * 5 * fact(2 * J - 2)
    rootbot = fact(2 * J + 3)

    b2_result = 6 * (J ** 2) * sqrt(roottop / rootbot) * f2(J, T)
    return b2_result


def B4(J, T):
    """
    Calculates the B4 term.
    J = nuclear spin
    T = temperature
    The following self tests only work for cobalt-60 with the specfic
    nuclear magnetic dipole moment, hyperfine field etc. etc.

    >>> B4(5, 0.00121856)
    1.1608400320878072
    >>> B4(5, 0.0030464)
    0.8609733528356602
    >>> B4(5, 0.0060928)
    0.37825253856770413
    >>> B4(5, 0.060928)
    0.0002457048732942802
    >>>

    """

    roottop = (2 * J + 1) * 9 * fact(2 * J - 4)
    rootbot = fact(2 * J + 5)

    b4_result = 70 * (J ** 4) * sqrt(roottop / rootbot) * f4(J, T)
    return b4_result


def Q2(y):
    """
    Calculates second order solid angle correction using arg y as the
    cos(alpha) value.

    >>> Q2(0)
    0.0
    >>> Q2(1)
    1.0
    >>> Q2(2)
    3.0
    >>> Q2(0.990268)
    0.985449355912
    >>>

    """
    q2_result = 0.5 * y * (1 + y)
    return q2_result


def Q4(y):
    """
    Calculates fourth order solid angle correction using arg y as the
    cos(alpha) value.

    >>> Q4(0)
    -0.0
    >>> Q4(1)
    1.0
    >>> Q4(2)
    18.75
    >>> Q4(0.990268)
    0.9520463139363521
    >>>

    """
    q4_result = 0.125 * y * (1 + y) * ((7 * (y ** 2)) - 3)
    return q4_result


def P2(x):
    """
    Calculates second order Legendre polynomial using arg x as the
    cos(theta) value.

    >>> P2(0)
    -0.5
    >>> P2(1)
    1.0
    >>> P2(2)
    5.5
    >>>

    """
    p2_result = (3 * (x ** 2) - 1) / 2
    return p2_result


def P4(x):
    """
    Calculates fourth order Legendre polynomial using x as the
    cos(theta) value.

    >>> P4(0)
    0.375
    >>> P4(1)
    1.0
    >>> P4(2)
    55.375
    >>>
    """
    p4_result = ((35 * (x ** 4)) - (30 * (x ** 2)) + 3) / 8
    return p4_result


def effect(T):
    """
    TODO(SCJK) Temperature in what units????
    """
    # Load data from the 'material' module.
    J = material.spin
    e = material.epsilon
    UF2 = material.U2F2
    UF4 = material.U4F4

    # Load cosine of detector angle and solid angle from modules.
    ct = detector_angle.cos_theta
    ca = solid_angle.cos_alpha

    W = 1 + \
        (1 - 3 * e) * Q2(ca) * B2(J, T) * UF2 * P2(ct) + \
        (1 - 10 * e) * Q4(ca) * B4(J, T) * UF4 * P4(ct)

    effect_result = -100 * (1 - W)

    return effect_result


# START
#T = 0.00364

running = True

while running is True:

    inputvalid = False

    while inputvalid is False:

        print("==========================================")
        inputdata = input("Effect (%)?")

        if inputdata is not "":
            try:
                effdata = float(inputdata)
                if effdata < -98.8 or effdata > -1.0:
                    print("Input out of valid range")
                    inputvalid = False
                else:
                    inputvalid = True
            except ValueError as e:
                print("Input probably not a number...")
        elif inputdata == "":
            raise RuntimeError("No data entered: program terminated.")

    NMAX = 100

    hot = 1.0
    cold = 0.0
    TOL = 0.0000001
    n = 0

    mid = (cold + hot) / 2.0

    while n < NMAX:

        print("h=", format(hot, "0.6f"), \
              "c=", format(cold, "0.6f"), \
              "m=", format(mid, "0.6f"), \
              "eff=", format(effect(mid), "0.3f"), \
              "n=", format(n, "d"))

        if (hot - cold) / 2 < TOL:
            break
        a = effdata
        b = effect(mid)
        if a < b:
            hot = mid
        else:
            cold = mid
        mid = (cold + hot) / 2
        n += 1

    print("T = ", format(mid * 1000, "0.2f"), " mK", sep='')
