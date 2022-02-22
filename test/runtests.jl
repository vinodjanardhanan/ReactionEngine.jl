using ReactionEngine
using Test

@testset "ReactionEngine.jl" begin
    # Write your tests here.
    
    @testset "Validate Inspect" begin
        retcode = inspect("inspect/inspect.xml","lib/")
        @test retcode == Symbol("Success")
    end

    @testset "Validate Thermoprobe" begin
        retcode = thermoprobe("thermo_probe/thermo_probe.xml","lib/")
        @test retcode == Symbol("Success")
    end

    @testset "Validate Transport properties" begin
        retcode = transport_properties("transport/transport.xml","lib/")
        @test retcode == Symbol("Success")
    end


    @testset "Validate CSTR" begin
        retcode = cstr("cstr/cstr.xml","lib/")
        @test retcode == Symbol("Success")
    end

    @testset "Validate Plug" begin
        retcode = plug("plug/plug.xml","lib/")
        @test retcode == Symbol("Success")
    end

    @testset "Validate Batch" begin
        retcode = batch("batch/batch.xml","lib/")
        @test retcode == Symbol("Success")
    end

    @testset "Validate Equilibrium" begin
        retcode = equilibrate("equil/ch4/equil.xml","lib/")
        @test retcode == Symbol("Success")
    end

end
