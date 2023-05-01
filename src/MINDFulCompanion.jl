module MINDFulCompanion

using MINDFul, Unitful, UUIDs

using Graphs, WrappedMultiGraphs, AttributeGraphs, NestedGraphs
import MetaGraphsNext as MGN
import MetaGraphs as MG

using DocStringExtensions

const MINDF = MINDFul

include("CompilationAlgorithms/compilationalgorithms.jl")
include("utils.jl")
include("defaultbehavior.jl")

end
