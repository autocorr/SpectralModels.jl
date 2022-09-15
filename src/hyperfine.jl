module Hyperfine
# FIXME
# - target performance of ~11us

import SpectralModels: predict, predict!
import SpectralModels.Core: getnu

export Satellites,
       HfSpectrum,
       HfWorkspace,
       HfModelParams,
       predict!,
       predict,
       getnu

using LoopVectorization

using ..Core
using ..Constants: CKMS, H, KB, TCMB


# Values to use for interpolating `1 / (exp(x) - 1)`
const T0_SIZE = 1000
const T0_LO = H * 23.0e9 / KB  # Frequency in Hz
const T0_HI = H * 28.0e9 / KB
const T0_XMIN = T0_LO / 8.0    # Upper Tex = 8.0 K
const T0_XMAX = T0_HI / 2.7    # Lower Tex = 2.7 K
const T0_X = [i * (T0_XMAX - T0_XMIN) / T0_SIZE for i=1:T0_SIZE]
const T0_Y = @. 1.0 / expm1(T0_X)
const T0_INV_DX = 1.0 / (T0_X[2] - T0_X[1])


"""
# Fields
- `n`:       Number of hyperfine transitions.
- `voff`:    Velocity offsets of each hyperfine transition.
- `tau_wts`: Optical depth weigth for each hyperfine transition.
"""
struct Satellites{V <: AbstractArray}
    n::Int
    voff::V
    tau_wts::V
    function Satellites(voff::V, tau_wts::V) where {V <: AbstractArray}
        @assert size(voff) == size(tau_wts)
        n = length(voff)
        new{V}(n, voff, tau_wts)
    end
end


struct HfSpectrum{F <: Real, V <: AbstractArray}
    spectrum::SpectrumData{F, V}
    transition::Transition{F}
    satellites::Satellites{V}
end


"""
# Fields
- `pred`:    Predicted model brightness temperature per channel.
- `tarr`:    Array of summed optical depth values per channel.
- `tbg_arr`: Background radiation temperature values per channel.
"""
struct HfWorkspace{V <: AbstractArray}
    pred::V
    tarr::V
    tbg_arr::V
    function HfWorkspace(pred::V, tarr::V, tbg_arr::V) where
            {V <: AbstractArray}
        @assert size(pred) == size(tarr) == size(tbg_arr)
        new{V}(pred, tarr, tbg_arr)
    end
end

function HfWorkspace(hfs::HfSpectrum)
    pred = zero(hfs.spectrum.data)
    tarr = zero(hfs.spectrum.data)
    tbg_arr = zero(hfs.spectrum.data)
    for i in eachindex(tbg_arr)
        nu = getnu(hfs, i)
        T0 = H * nu / KB
        tbg_arr[i] = one(eltype(tbg_arr)) / expm1(T0 / TCMB)
    end
    HfWorkspace(pred, tarr, tbg_arr)
end


"""
# Fields
- `vcen`:     Line velocity offset
- `tex`:      Line excitation temperature
- `tau_main`: Main line optical depth
- `sigm`:     Line velocity dispersion
"""
struct HfModelParams{F <: Real}
    vcen::F
    tex::F
    tau_main::F
    sigm::F
    function HfModelParams(
            vcen::F, tex::F, tau_main::F, sigm::F) where {F <: Real}
        @assert sigm > 0
        new{F}(vcen, tex, tau_main, sigm)
    end
end


@inline function getnu(hfs::HfSpectrum, i)
    getnu(hfs.spectrum, i)
end


@doc raw"""
Predict the spectral profile for a line given the excitation temperature
and main-line optical depth. The function results in a mutation of the
``Spectrum.pred`` array. The slabs are assumed to be optically thin with
respect to each other.

Brightness temperature is predicted according to the expression:

```math
T_\mathrm{b} =
    \frac{h \nu}{k} \left( \frac{1}{\exp(h\nu/kT_\mathrm{ex}} - 1)
    - \frac{1}{\exp(h\nu/kT_\mathrm{bg}) - 1} \right)
    \times \left(1 - \exp(-\tau) \right)
```
"""
function predict!(
        work::HfWorkspace,
        hfs::HfSpectrum,
        params::HfModelParams,
        approx::Union{Exact,Approximate}=Approximate()
)
    calc_tau_profile!(approx, work, hfs, params)
    calc_brightness_temp_profile!(work, hfs, params.tex)
    work
end

function predict(hfs::HfSpectrum, params::HfModelParams; kwargs...)
    work = HfWorkspace(hfs)
    predict!(work, hfs, args...; kwargs...)
    work.pred
end
function predict(hfs::HfSpectrum, args...; kwargs...)
    params = HfModelParams(args...)
    predict(hfs, params; kwargs...)
end


"""
Calculate the velocity/frequency related constants for the hyperfine
transitions.
"""
function calc_tau_profile!(approx, work::HfWorkspace, hfs::HfSpectrum,
        params::HfModelParams)
    tarr = work.tarr
    tarr .= zero(eltype(tarr))
    for (voff, tau_wt) in zip(hfs.satellites.voff, hfs.satellites.tau_wts)
        # Hyperfine (HF) dependent values
        hf_freq   = (one(voff) - voff / CKMS) * hfs.transition.rest_freq
        hf_width  = params.sigm / CKMS * hf_freq
        hf_offset = params.vcen / CKMS * hf_freq
        hf_nucen  = hf_freq - hf_offset
        hf_tau    = params.tau_main * tau_wt
        hf_idenom = 0.5 / hf_width^2
        add_component_tau!(approx, tarr, hfs.spectrum, hf_nucen, hf_idenom, hf_tau)
    end
    tarr
end

@inline function add_component_tau!(::Exact,
        tarr, s::SpectrumData, hf_nucen, hf_idenom, hf_tau)
    @turbo for i in eachindex(tarr)
        nu = getnu(s, i)
        tarr[i] += hf_tau * exp(-(nu - hf_nucen)^2 * hf_idenom)
    end
    tarr
end

@doc raw"""
For each HF line, sum the optical depth in each channel. The Gaussians are
approximated by only computing them within the range of ``\exp(-12.5)`` (5σ,
3.7e-6) away from the HF line center. Thus, inverting the expression:

```math
\exp(-ν\^2 × ξ) = exp(-12.5)
```

where ξ is ``1/2σ^2``, yields a cut-off half-width of ``\sqrt{12.5 / ξ}``.
"""
@inline function add_component_tau!(::Approximate,
        tarr, s::SpectrumData, hf_nucen, hf_idenom, hf_tau)
    nu_cutoff = sqrt(12.5 / hf_idenom)
    nu_lo = hf_nucen - s.freq_min - nu_cutoff
    nu_hi = hf_nucen - s.freq_min + nu_cutoff
    # Get the lower and upper indices then check bounds
    nu_lo_ix = ceil(Int, nu_lo/s.δchan)
    nu_hi_ix = ceil(Int, nu_hi/s.δchan)
    # If the profile lies outside of the range, then no-op
    if nu_hi_ix < 1 || nu_lo_ix > s.nchan
        return tarr
    end
    nu_lo_ix = max(      1, nu_lo_ix)
    nu_hi_ix = min(s.nchan, nu_hi_ix)
    # Calculate the Gaussian tau profile over the interval
    @turbo for i in nu_lo_ix:nu_hi_ix
        nu = getnu(s, i)
        tarr[i] += hf_tau * exp(-(nu - hf_nucen)^2 * hf_idenom)
    end
    tarr
end


function calc_brightness_temp_profile!(work::HfWorkspace, s::HfSpectrum, tex)
    pred = work.pred
    tarr = work.tarr
    tbg_arr = work.tbg_arr
    @turbo for i in eachindex(pred)
        T0 = H * getnu(s, i) / KB
        pred[i] += (
            (T0 / (exp(T0 / tex) - one(T0)) - T0 * tbg_arr[i])
            * (one(T0) - exp(-tarr[i]))
        )
    end
    work
end


@doc raw"""
Use a linear interpolation to approximate the function:

```math
f(x) = \frac{1}{\exp(x) - 1}
```

over the domain in ``x`` from (0.138, 0.498). This corresponds to frequency
values between 23-28 GHz and excitation temperatures between 2.7-8.0 K.
Outside of this interval, the exact solution is used. Using a linear
interpolation with N=1000 results in a relative numerical precision of ~1.8e-6
(|f'-f|/f) and a ~1.3x speed increase.
"""
@inline function iemtex_interp(x)
    if T0_XMIN < x < T0_XMAX
        i_lo = ceil(Int, (x - T0_XMIN) * T0_INV_DX)
        i_hi = i_lo + one(i_lo)
        x_lo = T0_X[i_lo]
        y_lo = T0_Y[i_lo]
        y_hi = T0_Y[i_hi]
        slope = (y_hi - y_lo) * T0_INV_DX
        return slope * (x - x_lo) + y_lo
    else
        return one(x) / (exp(x) - one(x))
    end
end


end  # module
