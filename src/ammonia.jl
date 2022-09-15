module Ammonia

import SpectralModels: predict, predict!

export AmmoniaSpectrum,
       predict!,
       predict

using StaticArrays

using ..Constants: CCMS, CKMS, H, KB
using ..Core
using ..Hyperfine


# Ammonia rotation constants [Splat ID: 01709]
# Coudert & Roueff (2006), A&A 449 855-859
const BROT = 298192.92e6
const CROT = 186695.86e6
# Poynter & Kakar (1975), ApJS, 29, 87
#const BROT = 298117.06e6
#const CROT = 186726.36e6

# The following values are taken from:
#   `pyspeckit.pyspeckit.spectrum.models.ammonia`

# Levels to calculate the partition function over
const NPART = 50
const NPARA = 34  # out of 51
const NORTH = 17  # out of 51
# Number of rotational levels: (1,1) thru (9,9)
const N_LEVELS = 9
# Number of model parameters
const N_PARAMS = 6
# J quantum numbers for para states
#const JORTH = [i for i=0:NPART if i % 3 == 0]
const JORTH = SVector{17}([i for i=0:NPART if i % 3 == 0])
# J quantum numbers for ortho states
#const JPARA = [i for i=0:NPART if i % 3 != 0]
const JPARA = SVector{34}([i for i=0:NPART if i % 3 != 0])
# Number of hyperfine transitions per level
const NHF = [
        18,  # (1,1)
        21,  # (2,2)
        26,  # (3,3)
        7,   # (4,4)
        7,   # (5,5)
        7,   # (6,6)
        7,   # (7,7)
        7,   # (8,8)
        1,   # (9,9)
]

# Ammonia inversion transition rest frequencies
# In Hz, note that for (1,1) Erik's custom freq is used,
# see pyspeckit issue 91.
const NU = [
        23.6944955e9,    # (1,1)
        23.722633335e9,  # (2,2)
        23.8701296e9,    # (3,3)
        24.1394169e9,    # (4,4)
        24.53299e9,      # (5,5)
        25.05603e9,      # (6,6)
        25.71518e9,      # (7,7)
        26.51898e9,      # (8,8)
        27.477943e9,     # (9,9)
]

# Transition Einstein A values
# Einstein A values are calculated using the expression:
#   A = 64π^4 / (3h * c^3) * ν^3 * μ0^2 * (j / (j + 1))
# The values computed below use the dipole moment found in
# Coudert & Roueff (2006), A&A 449 855-859
#   mu0 = 1.471 D     (note 1 Debye is 1e-18 statC cm)
# in addition to up-to-date values for `h` and `c`.  These are
# consistent with the JPL values on Splatalogue to approximately
# 4 significant digits.
const EA = [
        1.67524303e-07,  # (1,1)
        2.24162441e-07,  # (2,2)
        2.56915917e-07,  # (3,3)
        2.83423417e-07,  # (4,4)
        3.09910019e-07,  # (5,5)
        3.39590403e-07,  # (6,6)
        3.74750461e-07,  # (7,7)
        4.17525824e-07,  # (8,8)
        4.70284410e-07,  # (9,9)
]
# Values from pyspeckit, originally computed with:
#   mu0 = 1.476 D  (Poynter & Kakar 1975; pg. 9)
#EA = [
#        1.712e-7,        # (1,1)
#        2.291e-7,        # (2,2)
#        2.625e-7,        # (3,3)
#        3.167e-7,        # (4,4)
#        3.099109e-07,    # (5,5)
#        3.395797e-07,    # (6,6)
#        3.747560e-07,    # (7,7)
#        4.175308e-07,    # (8,8)
#        2.602045e-07,    # (9,9)
#]

const TRANSITIONS = [
        begin
            para = i % 3 != 0
            Transition(i, para, NU[i], EA[i])
        end
        for i in 1:9
]

# Velocity offsets of the hyperfine transitions in km/s
const VOFF = [
   [  # (1,1)
        19.851300,  19.315900,  7.8866900,  7.4696700,  7.3513200,
        0.4604090,  0.3220420, -0.0751680, -0.2130030,  0.3110340,
        0.1922660, -0.1323820, -0.2509230, -7.2334900, -7.3728000,
        -7.815260, -19.411700, -19.550000,
], [  # (2,2)
        26.526300,  26.011100,  25.950500,  16.391700,  16.379300,
        15.864200,  0.5625030,  0.5284080,  0.5237450,  0.0132820,
       -0.0037910, -0.0132820, -0.5018310, -0.5313400, -0.5890800,
       -15.854700, -16.369800, -16.382200, -25.950500, -26.011100,
       -26.526300,
], [  # (3,3)
        29.195098,  29.044147,  28.941877,  28.911408,  21.234827,
        21.214619,  21.136387,  21.087456,  1.0051220,  0.8060820,
        0.7780620,  0.6285690,  0.0167540, -0.0055890, -0.0134010,
       -0.6397340, -0.7445540, -1.0319240, -21.125222, -21.203441,
       -21.223649, -21.076291, -28.908067, -28.938523, -29.040794,
       -29.191744,
], [  # (4,4)
        0.0, -30.49783692, 30.49783692, 0.0,  24.25907811,
       -24.25907811, 0.0,
], [  # (5,5)
        31.4053287863, 26.0285409785, 0.0, 0.0, 0.0, -25.9063412556,
       -31.2831290633,
], [  # (6,6)
        31.5872901302, 27.0406347326, 0.0, 0.0, 0.0, -26.9209859064,
       -31.4676413039,
], [  # (7,7)
        31.3605314845, 27.3967468359, 0.0, 0.0, 0.0, -27.5133287373,
       -31.477113386,
], [  # (8,8)
        30.9752235915, 27.4707274918, 0.0, 0.0, 0.0, -27.5837757531,
       -30.9752235915,
], [  # (9,9)
        0.0,
]]

# Optical depth weights of the hyperfine transitions, unitless.
# Weights are taken from pyspeckit and normalized.
# Note that the magnetic hyperfines are not included past (3,3).
# Note the precision on (1,1)-(3,3) is excessive, as the weights are only known
# to about 5-6 digits, but present for the purposes of testing the numerical
# accuracy against pyspeckit.
const TAU_WTS = [
   [  # (1,1)
        3.7036944444583331e-02, 7.4073888889166661e-02,
        4.6296430555354165e-02, 8.3333374999937510e-02,
        9.2594861107708343e-03, 1.8518472222291665e-02,
        9.2594861107708343e-03, 9.2594861107708343e-03,
        4.6296430555354165e-02, 1.6666475000287499e-02,
        1.4999977500033751e-01, 2.3333315000027499e-01,
        1.6666475000287499e-02, 4.6296430555354165e-02,
        9.2594861107708343e-03, 8.3333374999937510e-02,
        3.7036944444583331e-02, 7.4073888889166661e-02,
], [  # (2,2)
        3.3333014814319341e-03, 2.9999713332887409e-02,
        1.6666507407159671e-02, 2.9629434979121079e-02,
        2.0741161893659245e-02, 1.4811134150653125e-03,
        1.6666507407159671e-02, 9.2593477367631464e-03,
        8.4654390943867397e-03, 2.1296340535048242e-01,
        3.9788439670906156e-01, 1.1666714444518766e-01,
        9.2593477367631464e-03, 8.4654390943867397e-03,
        1.6666507407159671e-02, 1.4811134150653125e-03,
        2.0741161893659245e-02, 2.9629434979121079e-02,
        1.6666507407159671e-02, 2.9999713332887409e-02,
        3.3333014814319341e-03,
], [  # (3,3)
        1.0733009496302131e-02, 7.3598529604831297e-03,
        3.0055577436436044e-03, 4.8085422957419802e-03,
        5.8220646798827188e-03, 7.7475821627062281e-03,
        4.3472933350838039e-03, 1.0143100958382566e-02,
        1.6829022799877465e-02, 9.0910682245853580e-03,
        9.4700450746138028e-03, 8.2989803509693240e-03,
        2.5670824033959128e-01, 4.0182836637346286e-01,
        1.5524222134698701e-01, 8.2989803509693240e-03,
        9.4700450746138028e-03, 1.6829022799877465e-02,
        4.3472933350838039e-03, 7.7475821627062281e-03,
        5.8220646798827188e-03, 1.0143100958382566e-02,
        4.8085422957419802e-03, 3.0055577436436044e-03,
        7.3598529604831297e-03, 1.0733009496302131e-02,
], [  # (4,4)
        0.2431, 0.0162, 0.0162, 0.3008, 0.0163, 0.0163, 0.3911,
], [  # (5,5)
        0.0109080940831, 0.0109433143618, 0.311493418617, 0.261847767275,
        0.382955997218,  0.0109433143618, 0.0109080940831,
], [  # (6,6)
        0.0078350431801, 0.00784948916416, 0.317644539734, 0.274246689798,
        0.376739705779, 0.00784948916416, 0.0078350431801,
], [  # (7,7)
        0.00589524944656, 0.00590204051181, 0.371879455317, 0.321515700951,
        0.283010263815, 0.00590204051181, 0.00589524944656,
], [  # (8,8)
        0.00459516014524, 0.00459939439378, 0.324116135075, 0.289534720829,
        0.367960035019, 0.00459939439378, 0.00459516014524,
], [  # (9,9)
        1.0,
]]

const SATELLITES = [Satellites(VOFF[i], TAU_WTS[i]) for i in 1:9]

const PAR_NAMES = ["voff", "trot", "tex", "ntot", "sigm", "orth"]
const PAR_NAMES_SHORT = ["v", "Tk", "Tx", "N", "s", "o"]

const TEX_LABELS = [
        raw"$v_\mathrm{lsr}$",
        raw"$T_\mathrm{rot}$",
        raw"$T_\mathrm{ex}$",
        raw"$\log(N_\mathrm{p})$",
        raw"$\sigma_\mathrm{v}$",
        raw"$f_\mathrm{o}$",
]

const TEX_LABELS_WITH_UNITS = [
        raw"$v_\mathrm{lsr} \ [\mathrm{km\, s^{-2}}]$",
        raw"$T_\mathrm{rot} \ [\mathrm{K}]$",
        raw"$T_\mathrm{ex} \ [\mathrm{K}]$",
        raw"$\log(N) \ [\log(\mathrm{cm^{-2}})]$",
        raw"$\sigma_\mathrm{v} \ [\mathrm{km\, s^{-1}}]$",
        raw"$f_\mathrm{o}$",
]


"""
    AmmoniaSpectrum(xarr, data, noise; trans_id=1)

# Arguments
- `xarr`: Frequency axis array. Channels must be in ascending order. **units**: Hz.
- `data`: Brightness temperature intensity values. **units**: K.
- `noise`: The brightness temperature baseline RMS noise level. **units**: K.
- `trans_id`: NH3 meta-stable transition ID, i.e.:

    1 -> (1,1)
    2 -> (2,2)
    ...
    9 -> (9,9)
"""
struct AmmoniaSpectrum{F <: Real, V <: AbstractArray}
    hfs::HfSpectrum{F, V}
    trans_id::Int
    function AmmoniaSpectrum(hfs::HfSpectrum{F, V}, trans_id) where
            {F <: Real, V <: AbstractArray}
        @assert trans_id in 1:9
        new{F, V}(hfs, trans_id)
    end
end
function AmmoniaSpectrum(xarr, data, noise::Real; trans_id::Integer=1)
    @assert noise > zero(noise)
    hfs = HfSpectrum(
            SpectrumData(xarr, data),
            TRANSITIONS[trans_id],
            SATELLITES[trans_id],
    )
    AmmoniaSpectrum(hfs, trans_id)
end


"""
    function swift_convert(tkin)

Convert a gas kinetic temperature in Kelvin to ammonia rotation temperature
using the "cold ammonia" approximation derived in Swift et al. (2005) by
equation A6.
    https://ui.adsabs.harvard.edu/abs/2005ApJ...620..823S/abstract
"""
@inline function swift_convert(tkin)
    tkin / (1.0 + (tkin / 41.18) * log(1.0 + 0.6 * exp(-15.7 / tkin)))
end


@inline function partition_level(j, trot)
    return (
            (2 * j + 1)
            * exp(-H * (BROT * j * (j + 1)
            + (CROT - BROT) * j^2) / (KB * trot))
    )
end


@inline function partition_func(para, trot)
    # NOTE could likely interpolate for improved performance
    Qtot = zero(trot)
    if para
        Qtot = sum(partition_level.(JPARA, trot))
    else
        Qtot = 2 * sum(partition_level.(JORTH, trot))
    end
    Qtot
end


function predict!(work::HfWorkspace, spec::AmmoniaSpectrum,
        trot, tex, ntot, sigm, voff, orth; cold=false, lte=false)
    s = spec.hfs
    t = s.transition
    work.pred .= zero(eltype(work.pred))
    trot = cold ? swift_convert(trot) : trot
    tex  = lte  ? trot : tex
    # Calculate the partition function and the level populations
    zlev = partition_level(t.id, trot)
    qtot = partition_func(t.para, trot)
    # Calculate the main line optical depth
    species_frac = t.para ? one(orth) - orth : orth
    pop_rotstate = exp10(ntot) * species_frac * zlev / qtot
    nu = t.rest_freq
    expterm = (
               (one(tex) - exp(-H * nu / (KB * tex))) /
               (one(tex) + exp(-H * nu / (KB * tex)))
    )
    fracterm = CCMS^2 * t.ea / (8π * nu^2)
    widthterm = CKMS / (sigm * nu * sqrt(2π))
    tau_main = pop_rotstate * fracterm * expterm * widthterm
    params = HfModelParams(voff, tex, tau_main, sigm)
    predict!(work, s, params)
end

function predict(spec::AmmoniaSpectrum, args...; kwargs...)
    work = HfWorkspace(spec.hfs)
    predict!(work, spec, args...; kwargs...)
end


end  # module
