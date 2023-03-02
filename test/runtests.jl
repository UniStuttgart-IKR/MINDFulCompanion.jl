using MINDFulCompanion, MINDFul
using Test, TestSetExtensions

using Graphs, GraphIO, NestedGraphs, NestedGraphsIO
using Unitful

const testdir =  dirname(@__FILE__)
const MINDFC = MINDFulCompanion
const MINDF = MINDFul

include("testutils.jl")

@testset ExtendedTestSet "MINDFulCompanion.jl" begin
    @includetests ["jointrmsaheuristic"]
end
