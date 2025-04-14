
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
function MINDF.compileintent!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation)
    defaultintradomaincompalgorithm = intradomaincompilationtemplate()
    compileintenttemplate!(ibnf, idagnode, uniformrandomalg;
        intradomainalgfun = defaultintradomaincompalgorithm,
        externaldomainalgkeyword = MINDF.getcompilationalgorithmkeyword(uniformrandomalg),
        prioritizesplitpathsfun = uniformrandomprioritizesplitpathsfun, 
        prioritizesplitnodesfun = uniformrandomprioritizesplitnodesfun,
        prioritizesplitbordernodesfun = getrandomsplitintentbordernode 
        )
end

# this is a combination of random known border node and shortest distance unknown
"""
$(TYPEDSIGNATURES)
"""
function getrandomsplitintentbordernode(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm)
    sourceglobalnode = getsourcenode(getintent(idagnode))
    destinationglobalnode = getdestinationnode(getintent(idagnode))
    # randomly pick a border node
    # TODO-tomorrow
    dglobalbordernode = getfirst(shuffle(getrng(intentcompilationalgorithm), getbordernodesasglobal(ibnf))) do globalbordernode
        getibnfid(globalbordernode) == getibnfid(destinationglobalnode)
    end
    # if unknown domain give it shortest distance border node
    if isnothing(dglobalbordernode)
        sourcelocalnode = getlocalnode(ibnf, sourceglobalnode)
        borderlocals = getbordernodesaslocal(ibnf);
        hopdists = Graphs.dijkstra_shortest_paths(getibnag(ibnf), sourcelocalnode).dists
        borderlocalminidx = argmin(hopdists[borderlocals])
        return getglobalnode(ibnf, borderlocals[borderlocalminidx])
    else
        return dglobalbordernode
    return end
end

"""
$(TYPEDSIGNATURES)
"""
function uniformrandomprioritizesplitpathsfun(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation, paths::Vector{Vector{LocalNode}})
    return randperm(getrng(uniformrandomalg), length(paths))
end

"""
$(TYPEDSIGNATURES)
"""
function uniformrandomprioritizesplitnodesfun(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation, path::Vector{LocalNode})
    return randperm(getrng(uniformrandomalg), length(path)-2)
end

"""
$(TYPEDSIGNATURES)

Interfaces required:
 - `getcandidatepathsnum -> Int`
"""
function uniformrandom!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation)
    # needed variables
    ibnag = getibnag(ibnf)
    idag = getidag(ibnf)
    idagnodeid = getidagnodeid(idagnode)
    intent = getintent(idagnode)
    sourceglobalnode = getsourcenode(intent)
    sourcelocalnode = getlocalnode(ibnag, sourceglobalnode)
    sourcenodeview = getnodeview(ibnag, sourcelocalnode)
    destinationglobalnode = getdestinationnode(intent)
    destlocalnode = getlocalnode(ibnag, destinationglobalnode)
    destnodeview = getnodeview(ibnag, destlocalnode)
    demandrate = getrate(intent)
    constraints = getconstraints(intent)

    returncode::Symbol = ReturnCodes.FAIL
    # start algorthim
    ## work around Graphs.jl bug
    if sourcelocalnode == destlocalnode
        yenstate = Graphs.YenState([u"0.0km"], [[destlocalnode]])
    else
        yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, getweights(ibnag), getcandidatepathsnum(uniformrandomalg))
    end

    lowlevelintentstoadd = LowLevelIntent[]
    ## define a TransmissionModuleCompatibility for the destination node
    transmissionmodulecompat = nothing
    opticalinitiateconstraint = getfirst(x -> x isa OpticalInitiateConstraint, constraints)
    if !isnothing(opticalinitiateconstraint)
        # find router port 
        yenidxs = randperm(length(yenstate.dists))
        for (dist, path) in zip(yenstate.dists[yenidxs], yenstate.paths[yenidxs])
            # find transmission module and mode
            spectrumslotsrange = getspectrumslotsrange(opticalinitiateconstraint)
            if length(path) > 1
                if getopticalreach(opticalinitiateconstraint) < dist
                    returncode = ReturnCodes.FAIL_OPTICALREACH_OPTINIT
                    continue
                end
                pathspectrumavailability = getpathspectrumavailabilities(ibnf, path)
                if !all(pathspectrumavailability[spectrumslotsrange])
                    returncode = ReturnCodes.FAIL_SPECTRUM
                    continue
                end
            end

            transmissionmodulecompat = gettransmissionmodulecompat(opticalinitiateconstraint)
            sourceadddropport = nothing
            opticalinitincomingnode = something(getlocalnode(ibnag, getglobalnode_input(opticalinitiateconstraint)))

            oxcadddropbypassspectrumllis = generatelightpathoxcadddropbypassspectrumlli(path, spectrumslotsrange; sourceadddropport, opticalinitincomingnode, destadddropport = nothing)
            foreach(oxcadddropbypassspectrumllis) do lli
                push!(lowlevelintentstoadd, lli)
            end
            
            # successful source-path configuration
            opticalterminateconstraint = getfirst(x -> x isa OpticalTerminateConstraint, constraints)
            if !isnothing(opticalterminateconstraint)
                # no need to do something more. add intents and return true
                foreach(lowlevelintentstoadd) do lli
                    addidagnode!(idag, lli; parentid = idagnodeid, intentissuer = MachineGenerated())
                end
                return ReturnCodes.SUCCESS
            else
                opticalincomingnode = length(path) == 1 ? opticalinitincomingnode : path[end-1]
                return uniformrandomintradomain_destination!(ibnf, idagnode, lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode, spectrumslotsrange, uniformrandomalg)
            end
        end
    else
        sourcerouterindex = getuniformrandomavailablerouterportindex(getrouterview(sourcenodeview), getrng(uniformrandomalg))
        if isnothing(sourcerouterindex)
            return ReturnCodes.FAIL_SRCROUTERPORT
        end
        sourcerouterportlli = RouterPortLLI(sourcelocalnode, sourcerouterindex)
        push!(lowlevelintentstoadd, sourcerouterportlli)

        for (dist, path) in zip(yenstate.dists, yenstate.paths)
            # find transmission module and mode
            sourceavailtransmdlidxs = getavailabletransmissionmoduleviewindex(sourcenodeview)
            sourcetransmissionmoduleviewpool = gettransmissionmoduleviewpool(sourcenodeview)
            for sourcetransmdlidx in sourceavailtransmdlidxs
                sourcetransmissionmodule = sourcetransmissionmoduleviewpool[sourcetransmdlidx]
                sourcetransmissiomodeidx = getuniformrandomtransmissionmode(sourcetransmissionmodule, demandrate, dist, getrng(uniformrandomalg))

                if isnothing(sourcetransmissiomodeidx)
                    returncode = ReturnCodes.FAIL_SRCTRANSMDL
                    continue
                end
                sourcetransmissionmode = gettransmissionmode(sourcetransmissionmodule, sourcetransmissiomodeidx)
                demandslotsneeded = getspectrumslotsneeded(sourcetransmissionmode)
                transmissionmoderate = getrate(sourcetransmissionmode)
                transmissionmodulename = getname(sourcetransmissionmodule)

                transmissionmodulecompat = TransmissionModuleCompatibility(transmissionmoderate, demandslotsneeded, transmissionmodulename)

                # find oxc configuration
                pathspectrumavailability = getpathspectrumavailabilities(ibnf, path)
                startingslot = randomfit(pathspectrumavailability, demandslotsneeded)
                if isnothing(startingslot)
                    returncode = ReturnCodes.Fail_SPECTRUM
                    continue
                end

                # are there oxc ports in the source ?
                sourceadddropport = getuniformrandomavailableoxcadddropport(sourcenodeview, getrng(uniformrandomalg))
                if isnothing(sourceadddropport)
                    returncode = ReturnCodes.FAIL_SRCOXCADDDROPPORT
                    continue
                end

                sourcetransmissionmodulelli = TransmissionModuleLLI(sourcelocalnode, sourcetransmdlidx, sourcetransmissiomodeidx, sourcerouterindex, sourceadddropport)
                push!(lowlevelintentstoadd, sourcetransmissionmodulelli)

                opticalinitincomingnode = nothing
                spectrumslotsrange = startingslot:(startingslot + demandslotsneeded - 1)
                oxcadddropbypassspectrumllis = generatelightpathoxcadddropbypassspectrumlli(path, spectrumslotsrange; sourceadddropport, opticalinitincomingnode, destadddropport = nothing)

                foreach(oxcadddropbypassspectrumllis) do lli
                    push!(lowlevelintentstoadd, lli)
                end
    
                # successful source-path configuration
                opticalterminateconstraint = getfirst(x -> x isa OpticalTerminateConstraint, constraints)
                if !isnothing(opticalterminateconstraint)
                    # no need to do something more. add intents and return true
                    foreach(lowlevelintentstoadd) do lli
                        addidagnode!(idag, lli; parentid = idagnodeid, intentissuer = MachineGenerated())
                    end
                    return ReturnCodes.SUCCESS
                else
                    # need to allocate a router port, a transmission module and mode, and an OXC configuration
                    opticalincomingnode = path[end-1]
                    return uniformrandomintradomain_destination!(ibnf, idagnode, lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode, spectrumslotsrange, uniformrandomalg)
                end
            end
        end
    end
    return returncode
end

"""
$(TYPEDSIGNATURES)
    Takes care of the final node (destination) for the case of no `OpticalTerminateConstraint`
"""
function uniformrandomintradomain_destination!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode::Int, spectrumslotsrange::UnitRange{Int}, uniformrandomalg::UniformRandomCompilation)
    ibnag = getibnag(ibnf)
    idag = getidag(ibnf)
    idagnodeid = getidagnodeid(idagnode)
    intent = getintent(idagnode)
    destinationglobalnode = getdestinationnode(intent)
    destlocalnode = getlocalnode(destinationglobalnode)
    destnodeview = getnodeview(ibnag, destlocalnode)

    # need to allocate a router port and a transmission module and mode
    destrouterindex = getuniformrandomavailablerouterportindex(getrouterview(destnodeview), getrng(uniformrandomalg))
    !isnothing(destrouterindex) || return ReturnCodes.FAIL_DSTROUTERPORT
    destrouterportlli = RouterPortLLI(destlocalnode, destrouterindex)
    push!(lowlevelintentstoadd, destrouterportlli)

    destavailtransmdlidxs = getavailabletransmissionmoduleviewindex(destnodeview)
    desttransmissionmoduleviewpool = gettransmissionmoduleviewpool(destnodeview)
    destavailtransmdlmodeidx = getuniformrandomcompatibletransmoduleidxandmodeidx(desttransmissionmoduleviewpool, destavailtransmdlidxs, transmissionmodulecompat, getrng(uniformrandomalg))
    !isnothing(destavailtransmdlmodeidx) || return ReturnCodes.FAIL_DSTTRANSMDL
    destavailtransmdlidx, desttransmodeidx = destavailtransmdlmodeidx[1], destavailtransmdlmodeidx[2] 

    # allocate OXC configuration
    destadddropport = getuniformrandomavailableoxcadddropport(destnodeview, getrng(uniformrandomalg))
    !isnothing(destadddropport) || return ReturnCodes.FAIL_DSTOXCADDDROPPORT
    oxclli = OXCAddDropBypassSpectrumLLI(destlocalnode, opticalincomingnode, destadddropport, 0, spectrumslotsrange)
    push!(lowlevelintentstoadd, oxclli)

    desttransmissionmodulelli = TransmissionModuleLLI(destlocalnode, destavailtransmdlidx, desttransmodeidx, destrouterindex, destadddropport)
    push!(lowlevelintentstoadd, desttransmissionmodulelli)

    foreach(lowlevelintentstoadd) do lli
        addidagnode!(idag, lli; parentid = idagnodeid, intentissuer = MachineGenerated())
    end
    return ReturnCodes.SUCCESS
end


"""
$(TYPEDSIGNATURES)
"""
function getuniformrandomavailablerouterportindex(nodeview::NodeView, rng::AbstractRNG)
    return getuniformrandomavailablerouterportindex(getrouterview(nodeview), rng)
end

"""
$(TYPEDSIGNATURES)

Return the uniformly random available router port index and `nothing` if non available.
"""
function prioritizerandomrouterports(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, node::LocalNode)
    routerview = getrouterview(getnodeview(getibnag(ibnf), node))
    reservedrouterports = getrouterportindex.(values(getreservations(routerview)))
    return filter(x -> x ∉ reservedrouterports, shuffle(getrng(intentcompilationalgorithm), 1:getportnumber(routerview)))
end

"""
$(TYPEDSIGNATURES)
"""
function getuniformrandomavailableoxcadddropport(nodeview::NodeView, rng::AbstractRNG)
    return getuniformrandomavailableoxcadddropport(getoxcview(nodeview), rng)
end

"""
$(TYPEDSIGNATURES)

Return the uniformly random available oxc add/drop port and `nothing` if none found
"""
function getuniformrandomavailableoxcadddropport(oxcview::OXCView, rng::AbstractRNG)
    reservedoxcadddropports = getadddropport.(values(getreservations(oxcview)))
    for adddropport in shuffle(1:getadddropportnumber(oxcview))
        adddropport ∉ reservedoxcadddropports && return adddropport
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return a integer tuple `(Int, Int)` denoting a uniformly randomly available and compatible transmission module index and its transmission mode index.
If non found return `nothing`.
"""
function getuniformrandomcompatibletransmoduleidxandmodeidx(transmissionmoduleviewpool::Vector{<:TransmissionModuleView}, availtransmdlidxs::Vector{Int}, transmissionmodulecompat::TransmissionModuleCompatibility, rng::AbstractRNG)
    for availtransmdlidx in shuffle(rng, availtransmdlidxs)
        transmissionmoduleview = transmissionmoduleviewpool[availtransmdlidx] 
        transmissionmodes = gettransmissionmodes(transmissionmoduleview)
        for transmodeidx in shuffle(rng, eachindex(transmissionmodes))
            if istransmissionmoduleandmodecompatible(transmissionmoduleview, transmodeidx, transmissionmodulecompat)
                return (availtransmdlidx, transmodeidx)
            end
        end
    end
    return nothing
end


"""
$(TYPEDSIGNATURES)

Return the index of a transmision mode that is uniform randomly selected with GBPS rate that can get deployed for the given demand rate and distance.
If non is found return `nothing`.
"""
function getuniformrandomtransmissionmode(transmissionmoduleview::TransmissionModuleView, demandrate::GBPSf, demanddistance::KMf, rng::AbstractRNG)
    transmodes = gettransmissionmodes(transmissionmoduleview)
    sps = randperm(rng, length(transmodes))
    for sp in sps
        transmode = transmodes[sp]
        getopticalreach(transmode) >= demanddistance && getrate(transmode) >= demandrate && return sp
    end
    return nothing
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

