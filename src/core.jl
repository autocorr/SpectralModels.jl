module Core

import SpectralModels: getnu

export SpectrumData,
       Transition,
       Exact,
       Approximate,
       getnu


struct Exact end
struct Approximate end


struct SpectrumData{F <: Real, V <: AbstractArray}
    nchan::Int
    δchan::F
    freq_min::F
    freq_max::F
    data::V
    function SpectrumData(data::V, δchan::F, freq_min::F) where
            {F <: Real, V <: AbstractArray}
        @assert ndims(data) == 1
        nchan = length(data)
        @assert nchan > 1
        freq_max = freq_min + δchan * (nchan - 1)
        @assert freq_max > freq_min
        new{F, V}(nchan, δchan, freq_min, freq_max, data)
    end
end
function SpectrumData(xarr::AbstractArray, data::AbstractArray)
    @assert size(xarr) == size(data)
    δchan = xarr[2] - xarr[1]
    @assert δchan > 0
    freq_min = xarr[1]
    SpectrumData(data, δchan, freq_min)
end


struct Transition{F <: Real}
    id::Int                       # Transition ID number (1 indexed)
    para::Bool                    # Is a para-transition?
    rest_freq::F                  # Rest frequency in Hz
    ea::F                         # Einstein A coefficient
end


@inline function getnu(s::SpectrumData, i)
    s.freq_min + s.δchan * (i - one(i))
end


end  # module
