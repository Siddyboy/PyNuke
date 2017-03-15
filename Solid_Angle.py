'''
Created on 9 Jan 2014
@author: kingsleys

This module defines the source-detector distance and the detector face
diameter. It then calculates the cosine of the maximum angle between
the source-detector axis and the source to outer diameter of detector
face. This cosine is the number used in the solid angle correction
function.

'''

import math

# Source to detector distance in mm.
detector_distance = 400

# Detector face diameter in mm.
detector_diameter = 75

# tan(alpha)
tan_alpha = detector_diameter / (2 * detector_distance)

# Calculate maximum angle in radians.
alpha_radians = math.atan(detector_diameter / (2 * detector_distance))

# Convert to degrees for a check.
alpha_degrees = alpha_radians * 180 / math.pi

# cos(alpha)
cos_alpha = math.cos(alpha_radians)
