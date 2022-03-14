# import Pkg; Pkg.activate(@__DIR__)

using Random
using LinearAlgebra
using Test
using JSON
using ArgParse

include(joinpath(@__DIR__, "..", "src", "gen_controllable.jl"))
include(joinpath(@__DIR__, "..", "src", "lqrdata.jl"))

function writetojson(data::LQRData{Nx,Nu}; indent=2, reg=1e-8, outputdir=joinpath(@__DIR__)) where {Nx,Nu}
    N = length(data.Q)
    Np = N*Nx + (N-1)*Nu
    Nd = N*Nx
    basename = "$Nx-$Nu-$(N)_"
    for form in ("banded", "kkt")
        A,b = build_Ab(data, form, reg=reg)
        json_data = Dict(
            "nx"=>Nx, 
            "nu"=>Nu, 
            "nprimals"=>Np,
            "nduals"=>Nd,
            "colptr"=>A.colptr, 
            "rowval"=>A.rowval,
            "nzval"=>A.nzval,
            "b"=>b
        )
        filename = basename * form * ".json"
        open(joinpath(outputdir, filename), "w") do file
            JSON.print(file, json_data, indent)
        end
    end
end

function runtest()
    data = rand(LQRData{10,5}, 2^8)
    A1,b1 = build_Ab(data, "banded", reg=1e-8)
    A2,b2 = build_Ab(data, "kkt", reg=1e-8)

    @test norm(A1) ≈ norm(A2)
    @test norm(b1) ≈ norm(b2)
    p = permvec(data)
    A0 = permute(A2, p, p)
    @test A0 == A1
    @test b1 == b2[p]
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--outputdir", "-o"
            help = "ouput directory for json files"
            arg_type = String
            default = joinpath(@__DIR__, "..", "data")
        "--runtests", "-t"
            help = "run test to check the matrix building"
            arg_type = Bool
            default = false
    end
    parsed_args = parse_args(ARGS, s)

    if parsed_args["runtests"]
        runtest()
    end

    sizes = [(2,1), (10,5), (14,4), (14,7)]
    lengths = [32, 64, 128, 256, 512, 1024]
    for (n,m) in sizes
        for N in lengths
            data = rand(LQRData{n,m}, N)
            writetojson(data, outputdir=parsed_args["outputdir"])
        end
    end
end

main()