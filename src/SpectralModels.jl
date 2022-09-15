module SpectralModels

function predict! end
function predict end
function getnu end

include("constants.jl")
include("core.jl")
include("hyperfine.jl")
include("ammonia.jl")

using .Core
using .Hyperfine
using .Ammonia

export Spectrum,
       AmmoniaSpectrum,
       getnu,
       Transition,
       Satellites,
       HfSpectrum,
       HfWorkspace,
       predict!,
       predict

end
