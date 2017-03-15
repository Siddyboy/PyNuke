'''
Created on 2 Jan 2014
@author: simon

This module sets up some physical constants for reference.

Uses SI units only.

'''

import math

BOLTZMANNS_CONSTANT = 1.38066e-23  # Units of J/K [3]

NUCLEAR_MAGNETON = 5.05084e-27  # Units of J/T [3]

PERMEABILITY_OF_FREE_SPACE = 4 * math.pi * 1.0e-7  # Units of T.m/A [4]

#======================================================================
# Ideas for improvement:
#    a) Make entries tuples with name, value, units?
#    b) Then printing them in debugging mode could pick up units and
#       names from here as well as values.
# SCJK 10 Jan 2014
#======================================================================
