using Profile
using InteractiveUtils

using BenchmarkTools
using Cthulhu

using SpectralModels
using SpectralModels.Ammonia
using SpectralModels.Hyperfine

include("util.jl")


const PARAMS = [15.0, 6.0, 15.0, 0.6, 0.0, 0.0]


function bench_psk_call()
    psk = Util.get_psk_ammonia()
    @benchmark begin
        spec = $psk.get_spec11()
    end
end


function ammonia_call()
    spec = Util.get_ammonia()
    work = Hyperfine.HfWorkspace(spec.hfs)
    predict!(work, spec, PARAMS...)
end


function profile_ammonia_call()
    spec = Util.get_ammonia()
    work = Hyperfine.HfWorkspace(spec.hfs)
    Profile.clear()
    @profile (for i=1:1_000_000; predict!(work, spec, PARAMS...); end)
    Profile.print()
    work
end


function bench_ammonia_call()
    spec = Util.get_ammonia()
    work = Hyperfine.HfWorkspace(spec.hfs)
    @benchmark begin
        predict!($work, $spec, 15.0, 6.0, 15.0, 0.6, 0.0, 0.0)
    end
end


function check_types()
    spec = Util.get_ammonia()
    work = Hyperfine.HfWorkspace(spec.hfs)
    @code_warntype getnu(spec.hfs.spectrum, 10)
    @code_warntype getnu(spec.hfs, 10)
    @descend_code_warntype predict!(work, spec, PARAMS...)
    params = Hyperfine.HfModelParams(0.0, 6.0, 7.6615, 0.6)
    approx = Hyperfine.Approximate()
    @code_warntype Hyperfine.calc_tau_profile!(approx, work, spec.hfs, params)
    @code_warntype Hyperfine.calc_brightness_temp_profile!(work, spec.hfs, params.tex)
    return nothing
end

