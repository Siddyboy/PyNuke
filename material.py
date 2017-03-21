'''
Created on 7 Jan 2014
@author: kingsleys

This module contains data for specific nuclear orientation materials.

In this order:

* Active nuclei species for gamma ray production.
* Host matrix (non-radioactive).
* Nuclear magnetic dipole moment in units of nuclear magnetons [1].
* Hyperfine field in Tesla [2].
* Nuclear spin angular momentum, dimensionless [3].
* Nuclear spin parity [3].
* U2F2 coefficients, dimensionless. [2]
* U4F4 coefficients, dimensionless. [2]
* epsilon. Fraction deviation from full magnetic saturation. [2]

'''

mat_cobalt = {
    'active': "cobalt-60",
    'host': "cobalt (hcp) single crystal",
    'moment': 3.799,
# Alternative, older value [2]. SCJK 7 Jan 2014
    'moment_old': 3.754,  
    'field': -21.92,
    'spin': 5,
    'parity': "+",
    'U2F2': -0.42056,
    'U4F4': -0.24280,
    'epsilon': 0.0
     }
    
mat_manganese = {
    'active': "manganese-54",
    'host': "nickel",
    'moment': 3.2819,
# Alternative older value [2]. SCJK 7 Jan 2014
    'moment_old': 3.302,
    'field': -32.56,
    'spin': 3,
    'parity': "+",
    'U2F2': -0.49486,
    'U4F4': -0.44669,
# Not sure this should be zero for 54Mn. SCJK 21 Mar 2017
    'epsilon': 0.0,  
    }

#    active = "cobalt-60"
#    host = "cobalt (hcp) single crystal"
#    moment = 3.799
#    # Alternative older value [2]. SCJK 7 Jan 2014
#    # moment = 3.754
#    field = -21.92
#    spin = 5
#    parity = "+"
#    U2F2 = -0.42056
#    U4F4 = -0.24280
#    epsilon = 0.0
