
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
    candidatepaths::Int
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
function getcandidatepaths(uniformrandomcompilation::UniformRandomCompilation)
    return uniformrandomcompilation.candidatepaths
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
function UniformRandomCompilation(seed::Int, candidatepaths::Int)
    return UniformRandomCompilation(seed, candidatepaths, MersenneTwister(seed))
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

A template compilation function that can be extended

Pass in the intent compilation algorithm `intentcompilationalgorithm`

Give in the following hook functions:
- `intradomainalgfun` is used as compilation algorithm for the intents handled internally. 
It should return a `Symbol` as a return code. 
Common return codes are found in `MINDFul.ReturnCodes`
```
intradomainalgfun(
    ibnf::IBNFramework, 
    idagnode::IntentDAGNode{<:ConnectivityIntent},
    intentcompilationalgorithm::IntentCompilationAlgorithm
) -> Symbol
```

- `prioritizesplitpathsfun` is called when optical reach is not enough to have a lightpath end-to-end to serve the intent.
When this happens several paths are considered which they can be broken in two.
Before we settle on a node as a split point, we need therefore to choose a path.
This function should return a vector of indices with decreasing priority of which path  of `paths` should be chosen.
```
prioritizesplitpathsfun(
    ibnf::IBNFramework,
    idagnode::IntentDAGNode,
    intentcompilationalgorithm::IntentCompilationAlgorithm,
    paths::Vector{Vector{LocalNode}}) -> Vector{Int}
) -> Vector{Int}
```

- `prioritizesplitnodesfun` is called when optical reach is not enough to have a lightpath end-to-end to serve the intent and a path to split was already selected.
The node selected will break the intent into two pieces with the node standing in between.
This function should return a vector of indices with decreasing priority of which node of `path` should be chosen.
```
prioritizesplitnodesfun(
    ibnf::IBNFramework,
    idagnode::IntentDAGNode,
    intentcompilationalgorithm::IntentCompilationAlgorithm,
    path::Vector{LocalNode}) -> Vector{Int}
) -> Vector{Int}
```

- `externaldomainalgkeyword` is called to select the border node to work as the source node for the delegated intent in a neighboring domain.
The function should return the node in a global representation.
```
externaldomainalgkeyword(
    ibnf::IBNFramework,
    idagnode::IntentDAGNode{<:ConnectivityIntent},
    intentcompilationalgorithm::IntentCompilationAlgorithm)
) -> GlobalNode
```
"""
function compileintenttemplate!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm; intradomainalgfun::F1, externaldomainalgkeyword::Symbol, prioritizesplitpathsfun::F2, prioritizesplitnodesfun::F3, prioritizesplitbordernodesfun::F4) where{F1<:Function, F2<:Function, F3<:Function, F4<:Function}
    sourceglobalnode = getsourcenode(getintent(idagnode))
    destinationglobalnode = getdestinationnode(getintent(idagnode))

    returncode::Symbol = ReturnCodes.FAIL

    if getibnfid(ibnf) == getibnfid(sourceglobalnode) == getibnfid(destinationglobalnode)
        # intra-domain
        returncode = intradomainalgfun(ibnf, idagnode, intentcompilationalgorithm)
        if returncode === ReturnCodes.FAIL_OPTICALREACH_OPTINIT || returncode === ReturnCodes.FAIL_OPTICALREACH
            # uncompile
            @assert MINDF.uncompileintent!(ibnf, getidagnodeid(idagnode)) 
            # find shortest distance neighbor j

            # get a node in between the shortest paths
            splitglobalnode = getsplitintentnode(ibnf, idagnode, intentcompilationalgorithm, prioritizesplitpathsfun, prioritizesplitnodesfun)
            returncode = splitandcompileintradomainconnecivityintent!(ibnf, idagnode, intentcompilationalgorithm, intradomainalgfun, splitglobalnode)
        end
        updateidagnodestates!(ibnf, idagnode)
    elseif getibnfid(ibnf) == getibnfid(sourceglobalnode) && getibnfid(ibnf) !== getibnfid(destinationglobalnode)
        # source intra-domain , destination cross-domain
        # border-node
        if isbordernode(ibnf, destinationglobalnode)
            #TODO-tomorrow
            returncode = splitandcompilecrossdomainconnectivityintent(ibnf, idagnode, intentcompilationalgorithm, intradomainalgfun, externaldomainalgkeyword,  destinationglobalnode)
        else
            # select border node
            destinationglobalbordernode = prioritizesplitbordernodesfun(ibnf, idagnode, intentcompilationalgorithm)

            returncode = splitandcompilecrossdomainconnectivityintent(ibnf, idagnode, intentcompilationalgorithm, intradomainalgfun, externaldomainalgkeyword,  destinationglobalbordernode)
        end
    end
    return returncode
end

"""
$(TYPEDSIGNATURES)
"""
function MINDF.compileintent!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation)
    compileintenttemplate!(ibnf, idagnode, uniformrandomalg;
        intradomainalgfun = uniformrandom!,
        externaldomainalgkeyword = MINDF.getcompilationalgorithmkeyword(uniformrandomalg),
        prioritizesplitpathsfun = uniformrandomprioritizesplitpathsfun , 
        prioritizesplitnodesfun = uniformrandomprioritizesplitnodesfun ,
        prioritizesplitbordernodesfun = getsplitintentbordernode 
        )
end

"""
$(TYPEDSIGNATURES)

Splits connectivity intent on `splitglobalnode`
"""
function splitandcompileintradomainconnecivityintent!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm,intradomainalgfun::F, splitglobalnode::GlobalNode) where {F}
    intent = getintent(idagnode)
    idag = getidag(ibnf)
    firsthalfintent = ConnectivityIntent(sourceglobalnode, splitglobalnode, getrate(intent), getconstraints(intent))
    firsthalfidagnode = addidagnode!(idag, firsthalfintent; parentid = getidagnodeid(idagnode), intentissuer = MachineGenerated())
    returncode = intradomainalgfun(ibnf, firsthalfidagnode, intentcompilationalgorithm)
    updateidagnodestates!(ibnf, firsthalfidagnode)
    issuccess(returncode) || return returncode

    secondhalfintent = ConnectivityIntent(splitglobalnode, destinationglobalnode, getrate(intent), filter(x -> !(x isa OpticalInitiateConstraint), getconstraints(intent)))
    secondhalfidagnode = addidagnode!(idag, secondhalfintent; parentid = getidagnodeid(idagnode), intentissuer = MachineGenerated())
    returncode = intradomainalgfun(ibnf, secondhalfidagnode, intentcompilationalgorithm)
    updateidagnodestates!(ibnf, secondhalfidagnode)
    return returncode
end

# this is a combination of random known border node and shortest distance unknown
function getsplitintentbordernode(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm)
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
function splitandcompilecrossdomainconnectivityintent(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, intradomainalgfun::F1, externaldomainalgkeyword::Symbol, mediatorbordernode::GlobalNode) where {F1}
    idag = getidag(ibnf)
    intent = getintent(idagnode)
    returncode::Symbol = ReturnCodes.FAIL

    internalintent = ConnectivityIntent(getsourcenode(intent), mediatorbordernode, getrate(intent), vcat(getconstraints(intent), OpticalTerminateConstraint()))

    internalidagnode = addidagnode!(idag, internalintent; parentid = getidagnodeid(idagnode), intentissuer = MachineGenerated())
    returncode = intradomainalgfun(ibnf, internalidagnode, intentcompilationalgorithm)
    updateidagnodestates!(ibnf, internalidagnode)

    issuccess(returncode) || return returncode
    
    # need first to compile that to get the optical choice
    opticalinitiateconstraint = getopticalinitiateconstraint(ibnf, getidagnodeid(internalidagnode))
    externalintent = ConnectivityIntent(mediatorbordernode, getdestinationnode(intent), getrate(intent), vcat(getconstraints(intent), opticalinitiateconstraint))
    externalidagnode = addidagnode!(idag, externalintent; parentid = getidagnodeid(idagnode), intentissuer = MachineGenerated())
    remoteibnfid = getibnfid(getdestinationnode(intent))
    internalremoteidagnode = remoteintent!(ibnf, externalidagnode, remoteibnfid)
    # getintent brings in the internal RemoteIntent
    externalremoteidagnodeid = getidagnodeid(getintent(internalremoteidagnode))

    # compile internalremoteidagnode
    remoteibnfhandler = getibnfhandler(ibnf, remoteibnfid)
    # compilationaglorithmkeyword = MINDF.getcompilationalgorithmkeyword(intentcompilationalgorithm)
    returncode = requestcompileintent_init!(ibnf, remoteibnfhandler, externalremoteidagnodeid, externaldomainalgkeyword, MINDF.getdefaultcompilationalgorithmargs(Val(externaldomainalgkeyword)))

    # check state of current internalremoteidagnode
    return returncode
end

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
        yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, getweights(ibnag), getcandidatepaths(uniformrandomalg))
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
function getuniformrandomavailablerouterportindex(routerview::RouterView, rng::AbstractRNG)
    reservedrouterports = getrouterportindex.(values(getreservations(routerview)))
    for routerportindex in shuffle(rng, 1:getportnumber(routerview))
        routerportindex ∉ reservedrouterports && return routerportindex
    end
    return nothing
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

"""
$(TYPEDSIGNATURES)

return the indices of the yenstate
prioritizesplitpathsfun(::IBNFramework, ::IntentDAGNode, ::IntentCompilationAlgorithm, ::Vector{Vector{LocalNode}}) -> Vector{Int}

return the indices of the nodes
prioritizesplitnodesfun(::IBNFramework, ::IntentDAGNode, ::IntentCompilationAlgorithm, ::Vector{LocalNode}) -> Vector{Int}

Return the [`GlobalNode`](@ref) to break up the [`ConnectiityIntent`](@ref) into
"""
function getsplitintentnode(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, prioritizesplitpathsfun::F1, prioritizesplitnodesfun::F2) where {F1<:Function, F2<:Function}
    ibnag = getibnag(ibnf)
    opticalinitiateconstraint = getfirst(x -> x isa OpticalInitiateConstraint, getconstraints(getintent(idagnode)))
    @assert !isnothing(opticalinitiateconstraint)
    opticalreach = getopticalreach(opticalinitiateconstraint)
    sourceglobalnode = getsourcenode(getintent(idagnode))
    sourcelocalnode = getlocalnode(ibnag, sourceglobalnode)
    destinationglobalnode = getdestinationnode(getintent(idagnode))
    destlocalnode = getlocalnode(ibnag, destinationglobalnode)
    yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, getweights(ibnag), getcandidatepaths(uniformrandomalg))
    # customize per yenstate priority order
    yenidxs = prioritizesplitpathsfun(ibnf, idagnode, intentcompilationalgorithm, yenstate.paths)
    # yenidxs = randperm(length(yenstate.dists))
    for (dist, path) in zip(yenstate.dists[yenidxs], yenstate.paths[yenidxs])
        # the accumulated distance from 2nd up to vorletzten node in path
        # diststopathnodes = accumulate(+, getindex.([getweights(ibnag)], path[1:end-2], path[2:end-1]))
        diststopathnodesidxs = prioritizesplitnodesfun(ibnf, idagnode, intentcompilationalgorithm, path)
        # diststopathnodesidxs = randperm(length(diststopathnodes))
        for nodeinpathidx in diststopathnodesidxs
            if opticalreach > diststopathnodes[nodeinpathidx]
                # +1 because we start measuring from the second node
                return getglobalnode(ibnag, path[nodeinpathidx+1])
            end
        end
    end
    return nothing
end

function uniformrandomprioritizesplitpathsfun(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation, paths::Vector{Vector{LocalNode}})
    return randperm(getrng(uniformrandomalg), length(paths))
end

function uniformrandomprioritizesplitnodesfun(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation, path::Vector{LocalNode})
    return randperm(getrng(uniformrandomalg), length(path))
end
