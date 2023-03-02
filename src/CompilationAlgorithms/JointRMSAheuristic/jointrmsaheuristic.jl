#=
# Initial methodology inspired from
# V. Gkamas, K. Christodoulopoulos and E. Varvarigos, "A Joint Multi-Layer Planning Algorithm for IP Over Flexible Optical Networks," 
# in Journal of Lightwave Technology, vol. 33, no. 14, pp. 2965-2977, 15 July15, 2015, doi: 10.1109/JLT.2015.2424920.
# 
# modified and extended to work with intents and multi domain networks
=#

export jointrmsagenerilizeddijkstra!, mlnodegraphtomlgraph, getmultilayernodegroups

include("typedef.jl")
include("multilayergraph.jl")
include("intentcompilation.jl")
