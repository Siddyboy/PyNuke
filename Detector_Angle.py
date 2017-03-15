'''
Created on 7 Jan 2014
@author: simon

This module defines the source-detector angle and does some trivial
calculations to pass values to the main code.

Specifically the angle is that between source magnetisation axis and
the detector's central axis.

For 60CohcpCo the magnetisation axis is parallel to the crystal's
c-axis and the c-axis has been made parallel with the long axis of the
source.

Not sure what happens with 54Mn. Presumably the externally applied
magnetic field axis defines the source axis.

'''

import math

# Source to detector angle in degrees.
theta_degrees = 0.0

# Convert to radians.
theta_radians = theta_degrees * math.pi / 180

# Take cosine.
cos_theta = math.cos(theta_radians)
