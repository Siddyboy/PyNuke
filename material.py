"""Contains materials data for specific nuclear orientation thermometry
materials.

The dictionary keys are:
'active': The active nuclei species for gamma ray production.
'host': The host matrix (non-radioactive).
'moment': The nuclear magnetic dipole moment in units of nuclear
     magnetons [1].
'moment_old': An alternative, older value for 'moment'. See ref [2].
'field': The hyperfine field in Tesla [2].
'spin': The nuclear spin angular momentum, dimensionless [3].
'parity': The nuclear spin parity [3].
'U2F2': Coefficients, dimensionless. [2]
'U4F4': Coefficients, dimensionless. [2]
'epsilon': Fraction deviation from full magnetic saturation. [2]
"""

mat_cobalt = {
    'active': "cobalt-60",
    'host': "cobalt (hcp) single crystal",
    'moment': 3.799,
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
    'moment_old': 3.302,
    'field': -32.56,
    'spin': 3,
    'parity': "+",
    'U2F2': -0.49486,
    'U4F4': -0.44669,
    # TODO(SCJK): Not sure epsilon should be zero for 54Mn. Probably needs more
    # research.
    'epsilon': 0.0,
}


# TODO(SCJK): I'm not sure that the following code is strictly necessary. Is
# it redundant?
def main():
    pass


if __name__ == "__main__":
    main()
