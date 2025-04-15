
"""
$(TYPEDEF)
$(TYPEDFIELDS)

Do everything uniformly randomly:
First, pick randomly among the candidate paths.
Then, pick randomly a spectrum slot, transmission module, transmission mode, and eventualy border node.
"""
struct UniformRandomCompilation <: MINDF.IntentCompilationAlgorithm
    "The starting seed"
    seed::Int
    "The number of candidate paths to choose upon"
    candidatepathsnum::Int
    "The random generator"
    rng::MersenneTwister
end

"""
$(TYPEDSIGNATURES)
"""
function getseed(uniformrandomcompilation::UniformRandomCompilation)
    return uniformrandomcompilation.seed
end

"""
$(TYPEDSIGNATURES)
"""
function getcandidatepathsnum(uniformrandomcompilation::UniformRandomCompilation)
    return uniformrandomcompilation.candidatepathsnum
end

"""
$(TYPEDSIGNATURES)
"""
function getrng(uniformrandomcompilation::UniformRandomCompilation)
    return uniformrandomcompilation.rng
end

"""
$(TYPEDSIGNATURES)
"""
function UniformRandomCompilation(seed::Int, candidatepathsnum::Int)
    return UniformRandomCompilation(seed, candidatepathsnum, MersenneTwister(seed))
end

"The keyword for [UniformRandomPolicy](@ref)"
const UniformRandomAlg = :uniformrandompolicy

"""
$(TYPEDSIGNATURES)

Give back the algorithm mapped to the symbol
"""
function MINDF.getcompilationalgorithmtype(s::Val{UniformRandomAlg})
    return UniformRandomCompilation
end

"""
$(TYPEDSIGNATURES)

Give back the symbol mapped to the algorithm
"""
function MINDF.getcompilationalgorithmkeyword(::Type{UniformRandomCompilation})
    return UniformRandomAlg
end

"""
$(TYPEDSIGNATURES)
"""
function MINDF.getdefaultcompilationalgorithmargs(s::Val{UniformRandomAlg})
    return (0, 5)
end

"""
$(TYPEDSIGNATURES)
"""
function MINDF.compileintent!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomcomp::UniformRandomCompilation)
    intradomaincompalgorithm = MINDF.intradomaincompilationtemplate(
        prioritizepaths = prioritizepaths_random,
        prioritizerouterport = MINDF.prioritizerouterports_first,
        prioritizetransmdlandmode = prioritizetransmdlmode_random,
        choosespectrum = choosespectrum_randomfit,
        chooseoxcadddropport = MINDF.chooseoxcadddropport_first,
    )
    MINDF.compileintenttemplate!(ibnf, idagnode, uniformrandomcomp;
        intradomainalgfun = intradomaincompalgorithm,
        externaldomainalgkeyword = MINDF.getcompilationalgorithmkeyword(uniformrandomcomp),
        prioritizesplitnodes = prioritizesplitnodes_random,
        prioritizesplitbordernodes = prioritizesplitbordernodes_random 
        )
end

"""
$(TYPEDSIGNATURES)

Return the a random [`GlobalNode`](@ref) contained in a random path.
The [`GlobalNode`](@ref) is used to break up the [`ConnectiityIntent`](@ref) into two.
Not several candidates are returned but only a single choice.
"""
function prioritizesplitnodes_random(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomcomp::UniformRandomCompilation)
    ibnag = getibnag(ibnf)
    opticalinitiateconstraint = getfirst(x -> x isa OpticalInitiateConstraint, getconstraints(getintent(idagnode)))
    @assert !isnothing(opticalinitiateconstraint)
    opticalreach = getopticalreach(opticalinitiateconstraint)
    sourceglobalnode = getsourcenode(getintent(idagnode))
    sourcelocalnode = getlocalnode(ibnag, sourceglobalnode)
    destinationglobalnode = getdestinationnode(getintent(idagnode))
    destlocalnode = getlocalnode(ibnag, destinationglobalnode)
    yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, getweights(ibnag), getcandidatepathsnum(intentcompilationalgorithm))
    # customize per yenstate priority order
    yenidxs = randperm(getrng(uniformrandomcomp), length(yenstate.paths))
    for (dist, path) in zip(yenstate.dists[yenidxs], yenstate.paths[yenidxs])
        # the accumulated distance from 2nd up to vorletzten node in path
        diststopathnodes = accumulate(+, getindex.([getweights(ibnag)], path[1:end-2], path[2:end-1]))
        nodeidxs = ramdperm(getrng(uniformrandomcomp), length(diststopathnodes))
        for nodeinpathidx in nodeidxs
            if opticalreach > diststopathnodes[nodeinpathidx]
                # +1 because we start measuring from the second node
                return [getglobalnode(ibnag, path[nodeinpathidx+1])]
            end
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
Return a single choice of a random border node [`GlobalNode`](@ref) and not several candidates.
"""
function prioritizesplitbordernodes_random(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomcomp::UniformRandomCompilation)
    ibnag = getibnag(ibnf)
    sourceglobalnode = getsourcenode(getintent(idagnode))
    sourcelocalnode = getlocalnode(ibnag, sourceglobalnode)
    destinationglobalnode = getdestinationnode(getintent(idagnode))
    borderlocals = getbordernodesaslocal(ibnf);
    # pick closest border node

    borderlocalsofdestdomain = filter(localnode -> getibnfid(getglobalnode(ibnag, localnode)) == getibnfid(destinationglobalnode), borderlocals)
    if !isempty(borderlocalsofdestdomain)
        # known domain
        randlocalnode = rand(getrng(uniformrandomcomp), borderlocalsofdestdomain)
        return [getglobalnode(ibnag, randlocalnode)]
    else
        # if unknown domain give it shortest distance border node
        borderlocalsofsrcdomain = filter(localnode -> getibnfid(getglobalnode(ibnag, localnode)) == getibnfid(sourceglobalnode), borderlocals)
        randlocalnode = rand(getrng(uniformrandomcomp), borderlocalsofsrcdomain)
        return [getglobalnode(ibnag, randlocalnode)]
    end
end

"""
$(TYPEDSIGNATURES)
"""
function prioritizepaths_random(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomcomp::UniformRandomCompilation)
    ibnag = getibnag(ibnf)
    distweights = getweights(ibnag)
    sourcelocalnode = getlocalnode(ibnag, getsourcenode(getintent(idagnode)))
    destlocalnode = getlocalnode(ibnag, getdestinationnode(getintent(idagnode)))
    if sourcelocalnode == destlocalnode
        yenstate = Graphs.YenState([u"0.0km"], [[destlocalnode]])
    else
        yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, distweights, getcandidatepathsnum(uniformrandomcomp))
    end
    return shuffle(getrng(uniformrandomcomp), yenstate.paths)
end


"""
$(TYPEDSIGNATURES)

Return the random available indices.
If non is find return `nothing`.
"""
function prioritizetransmdlmode_random(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomcomp::UniformRandomCompilation, node::LocalNode, path::Union{Nothing, Vector{LocalNode}}, transmdlcompat::Union{Nothing, TransmissionModuleCompatibility}=nothing)
    nodeview = getnodeview(getibnag(ibnf), node)
    demandrate = getrate(getintent(idagnode))
    availtransmdlidxs = getavailabletransmissionmoduleviewindex(nodeview)
    transmissionmoduleviewpool = gettransmissionmoduleviewpool(nodeview)
    returnpriorities = Tuple{Int,Int}[]
    transmdlperm = randperm(getrng(uniformrandomcomp), length(transmissionmoduleviewpool))
    filter!(i -> i ∈ availtransmdlidxs, transmdlperm)
    for transmdlidx in transmdlperm
        transmissionmodule = transmissionmoduleviewpool[transmdlidx]
        transmodes = gettransmissionmodes(transmissionmodule)
        transmodeidxs = randperm(getrng(uniformrandomcomp), length(transmodes))
        for transmodeidx in transmodeidxs
            transmode = transmodes[transmodeidx]
            if !isnothing(path) && isnothing(transmdlcompat)
                if getopticalreach(transmode) >= MINDF.getpathdistance(getibnag(ibnf), path) && getrate(transmode) >= demandrate
                    push!(returnpriorities, (transmdlidx, transmodeidx))
                end
            elseif isnothing(path) && !isnothing(transmdlcompat)
                if istransmissionmoduleandmodecompatible(transmissionmodule, transmodeidx, transmdlcompat)
                    push!(returnpriorities, (transmdlidx, transmodeidx))
                end
            end
        end
    end
    return returnpriorities
end

"""
$(TYPEDSIGNATURES)
"""
function choosespectrum_randomfit(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, path::Vector{LocalNode}, demandslotsneeded::Int)
    pathspectrumavailability = getpathspectrumavailabilities(ibnf, path)
    return randomfit(pathspectrumavailability, demandslotsneeded)
end


"""
$(TYPEDSIGNATURES)

Finds a random contiguous slot range of length `lengthrequire` that satisfies the `boolvec`.
Return the starting index of the range or `nothing` if none available
"""
function randomfit(boolvec::AbstractVector{Bool}, lengthrequire::Int)
    for i in shuffle(eachindex(boolvec))
        slotssatisfies = boolvec[i]
        for j in i:(i+lengthrequire-1)
            slotssatisfies &= j <= length(boolvec)
            slotssatisfies || break
            slotssatisfies &= boolvec[j]
            slotssatisfies || break
        end
        slotssatisfies && return i
    end
    return nothing
end

