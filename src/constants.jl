module Constants

export CKMS, CCMS, H, KB, TCMB

# Speed of light
const CKMS = 299792.458      # km/s
const CCMS = 29979245800.0   # cm/s

# Other physical constants in CGS; from `astropy.constants`
const H    = 6.62607015e-27  # erg s, Planck's constant
const KB   = 1.380649e-16    # erg/K, Boltzmann's constant
const TCMB = 2.72548         # K, T(CMB); Fixsen (2009) ApJ 707 916F

end  # module
