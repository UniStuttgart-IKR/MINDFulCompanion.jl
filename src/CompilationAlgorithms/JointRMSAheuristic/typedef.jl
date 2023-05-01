@enum(LinkType, opticallink, virtuallink, virtualtoopticallink, opticaltovirtuallink)
@enum(NodeType, routernodetype, oxcnodetype)

# I need to define nodes that could be
# - IP routers
# - transponders
# - OXC
# 
# and edges that can connect them

"""
$(TYPEDEF)
$(TYPEDFIELDS)

There are 4 types of links:
- optical link: connects 2 optical switches
- virtual link: connects 2 IP routers, a.k.a. lightpath
- inter-layer links:
 - optical to virtual: connects an optical switch of a multilayer node to the IP router
 - virtual to optical: connects an IP router of a multilayer node to the optical switch

A node might have many inter-layer links (each for every available transponder).
"""
mutable struct LinkCostVector{T<:Union{MINDF.TransmissionModuleView, Missing}, I<:Union{UUID, Missing},E<:AbstractEdge}
    const linktype::LinkType
    "The length of the link. It is non-zero only for optical links. It is zero for virtual and inter-layer links."
    const length::typeof(1.0u"km")
    "Cost of the transponder. It is non-zero only for inter-layer links."
    const costopt::Float64
    """
    Cost of the modular router. It is non-zero only for inter-layer links.
    It varies based on the state situation.
    If there are not readily available ports in the router, a new linecard must be added, and the cost of the linecard is considered here.
    If there are no available linecards, a new linecard chassis must be added and then this cost is considered instead.
    If a port is available in a linecard no cost is considered.
    Due to the dynamicity this variable must be modified every time a router configuration changes.
    It also depends on the requested demand rate.
    That's why it is always initialized to 0, and must be modified every time before pathfinding.
    """
    costip::Float64
    "Transmission module"
    const transmodule::T
    "Transmission modes for a particular transponder"
    const transmods::Vector{MINDFul.TransmissionProps}
    "1 if virtual link, 0 otherwise"
    const isvirtuallink::Int
    "if virtual link and already used, then the IBN DAG node UUID"
    const lpidnids::I
    "Spectrum slos. Valid only for optical links. All other (inter-layer and virtual links) have all elements `1`, i.e., `true`"
    const spectrum::BitVector
    "Physical path"
    const phypath::Vector{E}
end
LinkCostVector(::Val{opticallink}, length, spectrum, edg::E) where E<:AbstractEdge = 
    LinkCostVector(opticallink, length, 0.0, 0.0, missing, Vector{MINDFul.TransmissionProps}(), 0, missing, spectrum, Vector{E}([edg]))
LinkCostVector(::Val{opticaltovirtuallink}, costopt, edg::E) where E<:AbstractEdge = 
LinkCostVector(opticaltovirtuallink, 0.0u"km", costopt, 0.0, missing, Vector{MINDFul.TransmissionProps}(), 0, missing, BitVector(fill(true, 320)), Vector{E}([edg]))
LinkCostVector(::Val{virtualtoopticallink}, costopt, transmodule, edg::E) where E<:AbstractEdge = 
LinkCostVector(virtualtoopticallink, 0.0u"km", costopt, 0.0, transmodule, deepcopy(MINDF.gettransmissionmodes(transmodule)), 0, missing, BitVector(fill(true, 320)), Vector{E}([edg]))
LinkCostVector(::Val{virtuallink}, edgs::Vector{E}, spslots, idnuuid=missing) where E<:AbstractEdge = 
LinkCostVector(virtuallink, 0.0u"km", 0.0, 0.0, missing, Vector{MINDFul.TransmissionProps}(), 1, idnuuid,
               BitVector([el in spslots for el in 1:320]), edgs)

"""
$(TYPEDEF)
$(TYPEDFIELDS)

Adding up several `LinkCostVector` together.
"""
mutable struct PathCostVector{E<:AbstractEdge, T<:MINDFul.TransmissionModuleView}
    length::typeof(1.0u"km")
    costopt::Float64
    costip::Float64
    transmods::Vector{MINDFul.TransmissionProps}
    "Number of virtual links contained"
    virtuallinks::Int
    "Available spectrum on the path, considering contiguity"
    spectrum::BitVector
    "Last transmission module used"
    transmodule::Union{Missing,T}
    "The chosen rates, based on the `transmods` for all included virtual links"
    chosentransmodls::Vector{T}
    "Path composed by a vector of identifiable consecutive node pairs."
    path::Vector{E}
    "Physical path. Different only for virtual links"
    phypath::Vector{Edge{Int}}
end

getlength(pcv::PathCostVector) = pcv.length
getcostopt(pcv::PathCostVector) = pcv.costopt
getcostip(pcv::PathCostVector) = pcv.costip
getcost(pcv::PathCostVector) = pcv.costopt + pcv.costip
getvirtuallinks(pcv::PathCostVector) = pcv.virtuallinks
maximumtransrate(pcv::PathCostVector) = maximum(getrate.(pcv.chosentransmodls); init=0)

startend(pcv::PathCostVector) = (src(pcv.path[1]), dst(pcv.path[end]))

function PathCostVector(lpr::MINDF.LightpathRequirements, edg::Edge, edgspslots)
    trmodview = MINDF.TransmissionModuleView("unknown", 
                     MINDF.TransmissionModuleDummy([MINDF.TransmissionProps(lpr.optreach, lpr.rate, length(lpr.spslots))], 1, 0.0))
    spslots = BitVector([b && i in lpr.spslots for (i,b) in enumerate(edgspslots)])
    PathCostVector(lpr.dist, 0.0, 0.0, MINDF.gettransmissionmodes(trmodview),
                   0, spslots, trmodview, Vector{MINDF.TransmissionModuleView}(), [SingleMultiEdge(edg)], [edg])
end

PathCostVector(lcv::LinkCostVector, e::AbstractEdge) = PathCostVector(lcv.length, lcv.costopt, lcv.costip, lcv.transmods, lcv.isvirtuallink, lcv.spectrum, lcv.transmodule, Vector{MINDF.TransmissionModuleView{MINDF.TransmissionModuleDummy}}(), [e], lcv.phypath)

function PathCostVector(mlg::NestedGraph, es::AbstractEdge...)
    linkcost1 = getedgeattr(mlg, Tuple(es[1])...)
    pathcost1 = PathCostVector(linkcost1, es[1])
    if length(es) > 1
        add(pathcost1, mlg, es[2:end]...)
    else
        pathcost1
    end
end

function add(pcv::PathCostVector, mlg::NestedGraph, es::AbstractEdge...)
    pcvadd = deepcopy(pcv)
    for e in es
        add!(pcvadd, mlg, e)
    end
    pcvadd
end

function add(pcv::PathCostVector, mlg::NestedGraph, e::AbstractEdge, rate::Float64; minrate=false)
    pcvadd = deepcopy(pcv)
    add!(pcvadd, mlg, e, rate; minrate)
end

function add!(pcv::PathCostVector, mlg::NestedGraph, e::AbstractEdge, rate::Float64; minrate=false)
    lcv = getedgeattr(mlg, Tuple(e)...)
    if lcv.linktype == opticallink || lcv.linktype == virtuallink
        addoptvirt!(pcv, lcv, rate)
    elseif lcv.linktype == virtualtoopticallink
        addvirt2opt!(pcv, lcv)
    elseif lcv.linktype == opticaltovirtuallink
        addopt2virt!(pcv, lcv; minrate)
    end
    push!(pcv.path, e)
    push!(pcv.phypath, lcv.phypath...)
    return pcv
end

function addoptvirt!(pcv::PathCostVector, lcv::LinkCostVector, rate::Float64)
    pcv.length += lcv.length
    pcv.costopt += lcv.costopt
    pcv.costip += lcv.costip
    pcv.virtuallinks += lcv.isvirtuallink
    pcv.spectrum .&= lcv.spectrum
    filter!(pcv.transmods) do tp
        getoptreach(tp) >= pcv.length && MINDF.longestconsecutiveblock(==(true), pcv.spectrum) >= getfreqslots(tp) && getrate(tp) >= rate
    end
end
function addopt2virt!(pcv::PathCostVector, lcv::LinkCostVector; minrate=false)
    pcv.length += lcv.length
    pcv.costopt += lcv.costopt
    pcv.costip += lcv.costip
    pcv.virtuallinks += lcv.isvirtuallink
    if length(pcv.transmods) > 0
        if minrate
            rateind = findmin(tm -> tm.rate, pcv.transmods)[2]
        else
            rateind = findmax(tm -> tm.rate, pcv.transmods)[2]
        end
        chosentransmod = pcv.transmods[rateind]
        selectedmode = findfirst(==(chosentransmod), MINDF.gettransmissionmodes(pcv.transmodule))
        MINDF.setselection!(pcv.transmodule, selectedmode)
        push!(pcv.chosentransmodls, pcv.transmodule)
    end
    pcv.spectrum .&= lcv.spectrum
    pcv.transmods = empty(pcv.transmods)
    pcv.transmodule = missing
end
function addvirt2opt!(pcv::PathCostVector, lcv::LinkCostVector)
    pcv.length = 0.0u"km"
    pcv.costopt += lcv.costopt
    pcv.costip += lcv.costip
    pcv.transmods = deepcopy(lcv.transmods)
    pcv.transmodule =  lcv.transmodule
    pcv.virtuallinks += lcv.isvirtuallink
    pcv.spectrum .= true 
end

function finalizepathcostvec!(pcv::PathCostVector; minrate=false)
    if length(pcv.transmods) > 0
        if minrate
            rateind = findmin(tm -> tm.rate, pcv.transmods)[2]
        else
            rateind = findmax(tm -> tm.rate, pcv.transmods)[2]
        end
        chosentransmod = pcv.transmods[rateind]
        selectedmode = findfirst(==(chosentransmod), MINDF.gettransmissionmodes(pcv.transmodule))
        MINDF.setselection!(pcv.transmodule, selectedmode)
        push!(pcv.chosentransmodls, pcv.transmodule)
    end
end

"""
$(TYPEDSIGNATURES)

Checks if `pcv1` is dominating `pcv2`.
As a requirement they need to have the same end-nodes (source and destination).
"""
function isdominating(pcv1::PathCostVector, pcv2::PathCostVector; gothrough=Int[])
    startend(pcv1) == startend(pcv2) || error("`pcv1` and `pcv2` must have the same end nodes")

    if !isempty(gothrough) 
        gothroughflags1 = gothrough .∈ [pathify(pcv1.phypath)]
        gothroughflags2 = gothrough .∈ [pathify(pcv2.phypath)]
        all(gothroughflags1) && !all(gothroughflags2) && return true
        !all(gothroughflags1 .== gothroughflags2) && return false
    end
    # go here only if both paths gothrough the needed nodes or `gothrough` are empty

    (pcv1.length < pcv2.length && 
     pcv1.costopt + pcv1.costip <= pcv2.costopt + pcv2.costip &&
    pcv1.virtuallinks <= pcv2.virtuallinks && 
    maximumtransrate(pcv1) >= maximumtransrate(pcv2) &&
    all(pcv1.spectrum .>= pcv2.spectrum)) ||
    (pcv1.length <= pcv2.length && 
     pcv1.costopt + pcv1.costip < pcv2.costopt + pcv2.costip &&
    pcv1.virtuallinks <= pcv2.virtuallinks && 
    maximumtransrate(pcv1) >= maximumtransrate(pcv2) &&
    all(pcv1.spectrum .>= pcv2.spectrum)) ||
    (pcv1.length <= pcv2.length && 
     pcv1.costopt + pcv1.costip <= pcv2.costopt + pcv2.costip &&
    pcv1.virtuallinks < pcv2.virtuallinks && 
    maximumtransrate(pcv1) >= maximumtransrate(pcv2) &&
    all(pcv1.spectrum .>= pcv2.spectrum)) ||
    (pcv1.length <= pcv2.length && 
     pcv1.costopt + pcv1.costip <= pcv2.costopt + pcv2.costip &&
    pcv1.virtuallinks <= pcv2.virtuallinks && 
    maximumtransrate(pcv1) > maximumtransrate(pcv2) &&
    all(pcv1.spectrum .>= pcv2.spectrum)) ||
    (pcv1.length <= pcv2.length && 
     pcv1.costopt + pcv1.costip <= pcv2.costopt + pcv2.costip &&
    pcv1.virtuallinks <= pcv2.virtuallinks && 
    maximumtransrate(pcv1) >= maximumtransrate(pcv2) &&
    all(pcv1.spectrum .>= pcv2.spectrum) &&
    any(pcv1.spectrum .> pcv2.spectrum))
end

