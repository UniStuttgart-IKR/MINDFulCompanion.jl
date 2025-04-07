"""
$(TYPEDEF)
$(TYPEDFIELDS)

Do everything uniformly randomly:
First, pick randomly among the candidate paths.
Then, pick randomly a spectrum slot, transmission module, transmission mode, and eventualy border node.
"""
struct UniformRandomCompilation
    seed::Int
    candidatepaths::Int
    rng::MersenneTwister
end

"The keyword for [UniformRandomPolicy](@ref)"
const UniformRandoAlg = :uniformrandompolicy

"""
$(TYPEDSIGNATURES)

Give back the algorithm mapped to the symbol
"""
function MINDF.getcompilationalgorithmtype(s::Val{UniformRandoAlg})
    return UniformRandomCompilation
end

"""
$(TYPEDSIGNATURES)
"""
function MINDF.getdefaultcompilationalgorithmargs(s::Val{UniformRandoAlg})
    return (0,)
end

"""
$(TYPEDSIGNATURES)
"""
function MINDF.compileintent!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation)
    sourceglobalnode = getsourcenode(getintent(idagnode))
    destinationglobalnode = getdestinationnode(getintent(idagnode))

    if getibnfid(ibnf) == getibnfid(sourceglobalnode) == getibnfid(destinationglobalnode)
        # intra-domain
        compiledflag = uniformrandom!(ibnf, idagnode, uniformrandomalg)
        updateidagnodestates!(ibnf, idagnode)
        return compiledflag
    elseif getibnfid(ibnf) == getibnfid(sourceglobalnode) && getibnfid(ibnf) !== getibnfid(destinationglobalnode)
        # source intra-domain , destination cross-domain
        # border-node
        if isbordernode(ibnf, destinationglobalnode)
            return kspffcrossdomain!(ibnf, idagnode, uniformrandomalg, destinationglobalnode)
        else
            destinationglobalbordernode =  
            let 
                # randomly pick a border node
                dglobalbordernode = getfirst(getbordernodesasglobal(ibnf)) do globalbordernode
                    getibnfid(globalbordernode) == getibnfid(destinationglobalnode)
                end
                # if unknown domain give it randomly
                if isnothing(dglobalbordernode)
                    sourcelocalnode = getlocalnode(ibnf, sourceglobalnode)
                    borderlocals = getbordernodesaslocal(ibnf);
                    hopdists = Graphs.dijkstra_shortest_paths(getibnag(ibnf), sourcelocalnode).dists
                    borderlocalminidx = argmin(hopdists[borderlocals])
                    getglobalnode(ibnf, borderlocals[borderlocalminidx])
                else
                    dglobalbordernode
                end
            end
            return kspffcrossdomain!(ibnf, idagnode, uniformrandomalg, destinationglobalbordernode)
        end
    end
    return false
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

    # start algorthim
    ## work around Graphs.jl bug
    if sourcelocalnode == destlocalnode
        yenstate = Graphs.YenState([u"0.0km"], [[destlocalnode]])
    else
        yenstate = Graphs.yen_k_shortest_paths(ibnag, sourcelocalnode, destlocalnode, getweights(ibnag), uniformrandomalg.candidatepaths)
    end

    lowlevelintentstoadd = LowLevelIntent[]
    ## define a TransmissionModuleCompatibility for the destination node
    transmissionmodulecompat = nothing
    opticalinitiateconstraint = getfirst(x -> x isa OpticalInitiateConstraint, constraints)
    if !isnothing(opticalinitiateconstraint)
        # find router port 
        for (dist, path) in zip(yenstate.dists, yenstate.paths)
            # find transmission module and mode
            spectrumslotsrange = getspectrumslotsrange(opticalinitiateconstraint)
            if length(path) > 1
                getopticalreach(opticalinitiateconstraint) >= dist || continue
                pathspectrumavailability = getpathspectrumavailabilities(ibnf, path)
                all(pathspectrumavailability[spectrumslotsrange]) || continue
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
                return true
            else
                opticalincomingnode = length(path) == 1 ? opticalinitincomingnode : path[end-1]
                return kspffintradomain_destination!(ibnf, idagnode, lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode, spectrumslotsrange)
            end
        end
    else
        sourcerouterindex = getfirstavailablerouterportindex(getrouterview(sourcenodeview))
        !isnothing(sourcerouterindex) || return false
        sourcerouterportlli = RouterPortLLI(sourcelocalnode, sourcerouterindex)
        push!(lowlevelintentstoadd, sourcerouterportlli)

        for (dist, path) in zip(yenstate.dists, yenstate.paths)
            # find transmission module and mode
            sourceavailtransmdlidxs = getavailabletransmissionmoduleviewindex(sourcenodeview)
            sourcetransmissionmoduleviewpool = gettransmissionmoduleviewpool(sourcenodeview)
            for sourcetransmdlidx in sourceavailtransmdlidxs
                sourcetransmissionmodule = sourcetransmissionmoduleviewpool[sourcetransmdlidx]
                sourcetransmissiomodeidx = getlowestratetransmissionmode(sourcetransmissionmodule, demandrate, dist)

                !isnothing(sourcetransmissiomodeidx) || continue
                sourcetransmissionmode = gettransmissionmode(sourcetransmissionmodule, sourcetransmissiomodeidx)
                demandslotsneeded = getspectrumslotsneeded(sourcetransmissionmode)
                transmissionmoderate = getrate(sourcetransmissionmode)
                transmissionmodulename = getname(sourcetransmissionmodule)

                sourcetransmissionmodulelli = TransmissionModuleLLI(sourcelocalnode, sourcetransmdlidx, sourcetransmissiomodeidx)
                push!(lowlevelintentstoadd, sourcetransmissionmodulelli)

                transmissionmodulecompat = TransmissionModuleCompatibility(transmissionmoderate, demandslotsneeded, transmissionmodulename)

                # find oxc configuration
                pathspectrumavailability = getpathspectrumavailabilities(ibnf, path)
                startingslot = firstfit(pathspectrumavailability, demandslotsneeded)
                !isnothing(startingslot) || continue

                # are there oxc ports in the source ?
                sourceadddropport = getfirstavailableoxcadddropport(sourcenodeview)
                !isnothing(sourceadddropport) || continue

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
                    return true
                else
                    # need to allocate a router port, a transmission module and mode, and an OXC configuration
                    opticalincomingnode = path[end-1]
                    return kspffintradomain_destination!(ibnf, idagnode, lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode, spectrumslotsrange)
                end
            end
        end
    end
    return false
end

"""
$(TYPEDSIGNATURES)
    Takes care of the final node (destination) for the case of no `OpticalTerminateConstraint`
"""
function kspffintradomain_destination!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, lowlevelintentstoadd, transmissionmodulecompat, opticalincomingnode::Int, spectrumslotsrange::UnitRange{Int})
    ibnag = getibnag(ibnf)
    idag = getidag(ibnf)
    idagnodeid = getidagnodeid(idagnode)
    intent = getintent(idagnode)
    destinationglobalnode = getdestinationnode(intent)
    destlocalnode = getlocalnode(destinationglobalnode)
    destnodeview = getnodeview(ibnag, destlocalnode)

    # need to allocate a router port and a transmission module and mode
    destrouterindex = getfirstavailablerouterportindex(getrouterview(destnodeview))
    !isnothing(destrouterindex) || return false
    destrouterportlli = RouterPortLLI(destlocalnode, destrouterindex)
    push!(lowlevelintentstoadd, destrouterportlli)

    destavailtransmdlidxs = getavailabletransmissionmoduleviewindex(destnodeview)
    desttransmissionmoduleviewpool = gettransmissionmoduleviewpool(destnodeview)
    destavailtransmdlmodeidx = getfirstcompatibletransmoduleidxandmodeidx(desttransmissionmoduleviewpool, destavailtransmdlidxs, transmissionmodulecompat)
    !isnothing(destavailtransmdlmodeidx) || return false
    destavailtransmdlidx, desttransmodeidx = destavailtransmdlmodeidx[1], destavailtransmdlmodeidx[2] 
    desttransmissionmodulelli = TransmissionModuleLLI(destlocalnode, destavailtransmdlidx, desttransmodeidx)
    push!(lowlevelintentstoadd, desttransmissionmodulelli)

    # allocate OXC configuration
    destadddropport = getfirstavailableoxcadddropport(destnodeview)
    !isnothing(destadddropport) || return false
    oxclli = OXCAddDropBypassSpectrumLLI(destlocalnode, opticalincomingnode, destadddropport, 0, spectrumslotsrange)
    push!(lowlevelintentstoadd, oxclli)

    foreach(lowlevelintentstoadd) do lli
        addidagnode!(idag, lli; parentid = idagnodeid, intentissuer = MachineGenerated())
    end
    return true
end
