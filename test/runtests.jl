using MINDFulCompanion
using MINDFul
using Test, TestSetExtensions
using Graphs
import AttributeGraphs as AG
using JLD2, UUIDs
using Unitful, UnitfulData

import MINDFul: ReturnCodes

const MINDF = MINDFul
const MINDFC = MINDFulCompanion


import JET
import JET: @test_opt

TESTDIR = @__DIR__

# if you don't want JET tests do `push!(ARGS, "--nojet")` before `include`ing
RUNJET = !any(==("--nojet"), ARGS)

# get the test module from MINDFul
TM = Base.get_extension(MINDFul, :TestModule)
@test !isnothing(TM)

include("testsuite/uniformrandom.jl")
# include("testsuite/uniformrandom_multidomain.jl")

nothing
