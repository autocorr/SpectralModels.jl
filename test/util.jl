module Util

export get_psk_ammonia,
       get_ammonia

using PyCall
using SpectralModels

@pyinclude(joinpath(@__DIR__, "util.py"))


function get_psk_ammonia()
    return py"PskAmmoniaSpec()"
end


function get_ammonia()
    psk = Util.get_psk_ammonia()
    xarr = psk.xa11
    data = zero(xarr)
    spec = AmmoniaSpectrum(xarr, data, 0.1; trans_id=1)
end

end  # module
