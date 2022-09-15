using Test

using SpectralModels
using SpectralModels.Ammonia

include("util.jl")


@testset "SpectralModels.jl" begin
    @testset "PySpecKit validation" begin
        psk = Util.get_psk_ammonia()
        @test psk !== nothing
        @test psk.vchan ≈ 0.158
        @test psk.xa11[2] - psk.xa11[1] > 0
        spec11, spec22 = psk.get_spectra()
        @test maximum(spec11) ≈ 3.161_447_43
        @test maximum(spec22) ≈ 1.852_223_11
    end
    @testset "Ammonia validation" begin
        @test Ammonia.swift_convert(15) ≈ 14.023_487_576
        @test Ammonia.partition_level(1, 15.0) ≈ 0.635_853_477
        @test Ammonia.partition_func(false, 15.0) ≈ 2.003_701_059
        @test Ammonia.partition_func( true, 15.0) ≈ 0.703_882_578
        psk = Util.get_psk_ammonia()
        ammonia = Util.get_ammonia()
        spec11 = predict(ammonia, 15.0, 6.0, 15.0, 0.6, 0.0, 0.0)
        println(maximum(spec11.pred))
        @test spec11 !== nothing
        # NOTE These numbers are different (~0.01 K) than pyspeckit because the
        # latter uses a different CMB background temperature (2.7315 K) and
        # this can't be changed in the `ammonia` module because it is bound as
        # a default parameter in the function argument and is thus unaffected by
        # re-assignment.
        @test maximum(spec11.pred) ≈ 3.167_224_500
        @test maximum(spec11.tarr) ≈ 3.635_084_016
    end
end
