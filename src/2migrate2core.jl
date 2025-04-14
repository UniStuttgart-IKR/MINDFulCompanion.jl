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
    paths::Vector{Vector{LocalNode}})
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
    path::Vector{LocalNode})
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
    yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, getweights(ibnag), getcandidatepathsnum(intentcompilationalgorithm))
    # customize per yenstate priority order
    yenidxs = prioritizesplitpathsfun(ibnf, idagnode, intentcompilationalgorithm, yenstate.paths)
    # yenidxs = randperm(length(yenstate.dists))
    for (dist, path) in zip(yenstate.dists[yenidxs], yenstate.paths[yenidxs])
        # the accumulated distance from 2nd up to vorletzten node in path
        diststopathnodes = accumulate(+, getindex.([getweights(ibnag)], path[1:end-2], path[2:end-1]))
        @show diststopathnodes
        diststopathnodesidxs = prioritizesplitnodesfun(ibnf, idagnode, intentcompilationalgorithm, path)
        @show diststopathnodesidxs
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

"""
$(TYPEDSIGNATURES)

Splits connectivity intent on `splitglobalnode`
"""
function splitandcompileintradomainconnecivityintent!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm,intradomainalgfun::F, splitglobalnode::GlobalNode) where {F}
    sourceglobalnode = getsourcenode(getintent(idagnode))
    destinationglobalnode = getdestinationnode(getintent(idagnode))
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


"""
$(TYPEDSIGNATURES)

Interfaces required:
 - `getcandidatepathsnum -> Int`
"""
function getclosestsplitintentbordernode(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm)
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
    end
end

"""
$(TYPEDSIGNATURES)

Return a intra domain compilation algorithm as customized with the signature
```
intradomainalgfun(
    ibnf::IBNFramework, 
    idagnode::IntentDAGNode{<:ConnectivityIntent},
    intentcompilationalgorithm::IntentCompilationAlgorithm
) -> Symbol
```

Intra domain compilation algorithm template

The major selection process is made on the source.
The destination chooses the first that are compatible.

Interfaces needed:
```
getcandidatepathsnum(
    intentcompilationalgorithm::IntentCompilationAlgorithm)
 -> Int
```
```
prioritizepaths(
    ibnf::IBNFramework,
    idagnode::IntentDAGNode{<:ConnectivityIntent},
    intentcompilationalgorithm::IntentCompilationAlgorithm,
    paths::Vector{Int}
) -> Vector{Int}
```
```
prioritizerouterport(
    ibnf::IBNFramework,
    idagnode::IntentDAGNode{<:ConnectivityIntent},
    intentcompilationalgorithm::IntentCompilationAlgorithm,
    node::LocalNode
) -> Vector{Int}
```
```
prioritizetransmdlandmode(
    
) -> Vector{Int}
```
```
prioritizetransmode(
    
) -> Vector{Int}
```
```
choosespectrum(
    
) -> Vector{Int}
```
```
chooseoxcadddropport(
    
) -> Vector{Int}
```
```
chooserouterport(
    
) -> Int
```
```
choosetransmissionmdlandmode(
    
) -> Int
```
```
chooseoxcadddropport(
    
) -> Int
```
"""
function intradomaincompilationtemplate(;
    prioritizepaths = prioritizepaths_shortest,
    prioritizerouterport = prioritizerouterports_first,
    prioritizetransmdlandmode = prioritizetransmdlmode_cheaplowrate,
    choosespectrum = choosespectrum_firstfit,
    chooseoxcadddropport = chooseoxcadddropport_first,
    )
    return function(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm)
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
        # if sourcelocalnode == destlocalnode
        #     yenstate = Graphs.YenState([u"0.0km"], [[destlocalnode]])
        # else
        #     yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, getweights(ibnag), getcandidatepathsnum(intentcompilationalgorithm))
        # end
        # yenidxs = prioritizepaths(ibnfs, idagnode, intentcompilationalgorithm, paths)
        candidatepaths = prioritizepaths(ibnf, idagnode, intentcompilationalgorithm)

        lowlevelintentstoadd = LowLevelIntent[]
        ## define a TransmissionModuleCompatibility for the destination node
        transmissionmodulecompat = nothing
        opticalinitiateconstraint = getfirst(x -> x isa OpticalInitiateConstraint, constraints)
        if !isnothing(opticalinitiateconstraint)
            # template: prioritizepaths
            # yenidxs = randperm(length(yenstate.dists))
            for path in candidatepaths
                # find transmission module and mode
                spectrumslotsrange = getspectrumslotsrange(opticalinitiateconstraint)
                if length(path) > 1
                    if getopticalreach(opticalinitiateconstraint) < getpathdistance(ibnag, path)
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
                    return intradomaincompilationtemplate_destination!(ibnf, idagnode, intentcompilationalgorithm,lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode, spectrumslotsrange, prioritizerouterport, prioritizetransmdlandmode, chooseoxcadddropport)
                end
            end
        else
            # template prioritizerouterport
            returncode = ReturnCodes.FAIL_SRCROUTERPORT
            sourcerouteridxs = prioritizerouterport(ibnf, idagnode, intentcompilationalgorithm, sourcelocalnode)
            for sourcerouteridx in sourcerouteridxs
                sourcerouterportlli = RouterPortLLI(sourcelocalnode, sourcerouteridx)
                # TODO-NOW put pushes at the end and avoid abstract vector
                push!(lowlevelintentstoadd, sourcerouterportlli)

                for path in candidatepaths
                    # find transmission module and mode
                    sourcetransmissionmoduleviewpool = gettransmissionmoduleviewpool(sourcenodeview)
                    # template prioritizetransmdlandmode
                    returncode = ReturnCodes.FAIL_SRCTRANSMDL
                    for (sourcetransmdlidx, sourcetransmissiomodeidx) in prioritizetransmdlandmode(ibnf, idagnode, intentcompilationalgorithm, sourcelocalnode, path)
                    # for sourcetransmdlidx in sourceavailtransmdlidxs
                        sourcetransmissionmodule = sourcetransmissionmoduleviewpool[sourcetransmdlidx]
                        sourcetransmissionmode = gettransmissionmode(sourcetransmissionmodule, sourcetransmissiomodeidx)
                        demandslotsneeded = getspectrumslotsneeded(sourcetransmissionmode)
                        transmissionmoderate = getrate(sourcetransmissionmode)
                        transmissionmodulename = getname(sourcetransmissionmodule)

                        transmissionmodulecompat = TransmissionModuleCompatibility(transmissionmoderate, demandslotsneeded, transmissionmodulename)

                        # template choosespectrum
                        startingslot = choosespectrum(ibnf, idagnode, intentcompilationalgorithm, path, demandslotsneeded)
                        if isnothing(startingslot)
                            returncode = ReturnCodes.Fail_SPECTRUM
                            continue
                        end

                        # are there oxc ports in the source ?
                        # template chooseoxcadddropport
                        sourceadddropport = chooseoxcadddropport(ibnf, idagnode, intentcompilationalgorithm, sourcelocalnode)
                        if isnothing(sourceadddropport)
                            returncode = ReturnCodes.FAIL_SRCOXCADDDROPPORT
                            continue
                        end

                        sourcetransmissionmodulelli = TransmissionModuleLLI(sourcelocalnode, sourcetransmdlidx, sourcetransmissiomodeidx, sourcerouteridx, sourceadddropport)
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
                            return intradomaincompilationtemplate_destination!(ibnf, idagnode, intentcompilationalgorithm, lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode, spectrumslotsrange, prioritizerouterport, prioritizetransmdlandmode, chooseoxcadddropport)
                        end
                    end
                end
            end
        end
        return returncode
    end
end

"""
$(TYPEDSIGNATURES)
    Takes care of the final node (destination) for the case of no `OpticalTerminateConstraint`

Functions to pass in
```
chooserouterport(
     
 ) -> Int
```
```
choosetransmissionmdlandmode(
) -> Int
```

```
 chooseoxcadddropport(
     
 )-> Int
```
"""
function intradomaincompilationtemplate_destination!(
    ibnf::IBNFramework, 
    idagnode::IntentDAGNode{<:ConnectivityIntent},
    intentcompilationalgorithm::IntentCompilationAlgorithm,
    lowlevelintentstoadd,
    transmissionmodulecompat,
    opticalincomingnode::Int,
    spectrumslotsrange::UnitRange{Int},
    prioritizerouterport::F1,
    prioritizetransmdlmode::F2,
    chooseoxcadddropport::F3) where {F1<:Function, F2<:Function, F3<:Function}

    ibnag = getibnag(ibnf)
    idag = getidag(ibnf)
    idagnodeid = getidagnodeid(idagnode)
    intent = getintent(idagnode)
    destinationglobalnode = getdestinationnode(intent)
    destlocalnode = getlocalnode(destinationglobalnode)
    destnodeview = getnodeview(ibnag, destlocalnode)

    # need to allocate a router port and a transmission module and mode
    # template chooserouterport
    destrouteridxs = prioritizerouterport(ibnf, idagnode, intentcompilationalgorithm, destlocalnode)
    !isempty(destrouteridxs) || return ReturnCodes.FAIL_DSTROUTERPORT
    destrouteridx = first(destrouteridxs)
    destrouterportlli = RouterPortLLI(destlocalnode, destrouteridx)
    push!(lowlevelintentstoadd, destrouterportlli)

    destavailtransmdlidxs = getavailabletransmissionmoduleviewindex(destnodeview)
    desttransmissionmoduleviewpool = gettransmissionmoduleviewpool(destnodeview)
    # template intradomaincompilationtemplate_destination!
    destavailtransmdlmodeidxs = prioritizetransmdlmode(ibnf, idagnode, intentcompilationalgorithm, destlocalnode, nothing, transmissionmodulecompat)
    !isempty(destavailtransmdlmodeidxs) || return ReturnCodes.FAIL_DSTTRANSMDL
    destavailtransmdlmodeidx = first(destavailtransmdlmodeidxs)
    destavailtransmdlidx, desttransmodeidx = destavailtransmdlmodeidx[1], destavailtransmdlmodeidx[2] 

    # allocate OXC configuration
    # template chooseoxcadddropport
    destadddropport = chooseoxcadddropport(ibnf, idagnode, intentcompilationalgorithm, destlocalnode)
    !isnothing(destadddropport) || return ReturnCodes.FAIL_DSTOXCADDDROPPORT
    oxclli = OXCAddDropBypassSpectrumLLI(destlocalnode, opticalincomingnode, destadddropport, 0, spectrumslotsrange)
    push!(lowlevelintentstoadd, oxclli)

    desttransmissionmodulelli = TransmissionModuleLLI(destlocalnode, destavailtransmdlidx, desttransmodeidx, destrouteridx, destadddropport)
    push!(lowlevelintentstoadd, desttransmissionmodulelli)

    foreach(lowlevelintentstoadd) do lli
        addidagnode!(idag, lli; parentid = idagnodeid, intentissuer = MachineGenerated())
    end
    return ReturnCodes.SUCCESS
end

function prioritizepaths_shortest(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm)
    ibnag = getibnag(ibnf)
    distweights = getweights(ibnag)
    sourcelocalnode = getlocalnode(ibnag, getsourcenode(getintent(idagnode)))
    destlocalnode = getlocalnode(ibnag, getdestinationnode(getintent(idagnode)))
    # yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, distweights, getcandidatepathsnum(intentcompilationalgorithm))

    if sourcelocalnode == destlocalnode
        yenstate = Graphs.YenState([u"0.0km"], [[destlocalnode]])
    else
        yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, distweights, getcandidatepathsnum(intentcompilationalgorithm))
    end

    return yenstate.paths
end

function prioritizerouterports_first(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, node::LocalNode)
    routerview = getrouterview(getnodeview(getibnag(ibnf), node))
    reservedrouterports = getrouterportindex.(values(getreservations(routerview)))
    return filter(x -> x ∉ reservedrouterports, 1:getportnumber(routerview))
end


"""
$(TYPEDSIGNATURES)

Return the index with the lowest GBPS rate that can get deployed for the given demand rate and distance.
If non is find return `nothing`.
"""
function prioritizetransmdlmode_cheaplowrate(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, node::LocalNode, path::Union{Nothing, Vector{LocalNode}}, transmdlcompat::Union{Nothing, TransmissionModuleCompatibility}=nothing)
    nodeview = getnodeview(getibnag(ibnf), node)
    demandrate = getrate(getintent(idagnode))
    availtransmdlidxs = getavailabletransmissionmoduleviewindex(nodeview)
    transmissionmoduleviewpool = gettransmissionmoduleviewpool(nodeview)
    returnpriorities = Tuple{Int,Int}[]
    transmdlperm = sortperm(by = x -> MINDF.getcost(x) , transmissionmoduleviewpool)
    filter!(i -> i ∈ availtransmdlidxs, transmdlperm)
    for transmdlidx in transmdlperm
        transmissionmodule = transmissionmoduleviewpool[transmdlidx]
        transmodes = gettransmissionmodes(transmissionmodule)
        transmodeidxs = sortperm(transmodes; by = getrate)
        for transmodeidx in transmodeidxs
            transmode = transmodes[transmodeidx]
            if !isnothing(path) && isnothing(transmdlcompat)
                if getopticalreach(transmode) >= getpathdistance(getibnag(ibnf), path) && getrate(transmode) >= demandrate
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
function choosespectrum_firstfit(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, path::Vector{LocalNode}, demandslotsneeded::Int)
    pathspectrumavailability = getpathspectrumavailabilities(ibnf, path)
    return randomfit(pathspectrumavailability, demandslotsneeded)
end

"""
$(TYPEDSIGNATURES)

Return the uniformly random available oxc add/drop port and `nothing` if none found
"""
function chooseoxcadddropport_first(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, intentcompilationalgorithm::IntentCompilationAlgorithm, node::LocalNode)
    oxcview = getoxcview(getnodeview(getibnag(ibnf), node))
    reservedoxcadddropports = getadddropport.(values(getreservations(oxcview)))
    for adddropport in shuffle(1:getadddropportnumber(oxcview))
        adddropport ∉ reservedoxcadddropports && return adddropport
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)
"""
function getpathdistance(ibnag::IBNAttributeGraph, path::Vector{Int})
    ws = getweights(ibnag)
    return sum([getindex(ws, nodepair...) for nodepair in zip(path[1:end-1], path[2:end])])
end

