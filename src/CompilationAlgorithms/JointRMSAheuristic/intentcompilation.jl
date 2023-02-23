function jointrmsagenerilizeddijkstra!(ibn::IBN, idagnode::IntentDAGNode{R}, ::MINDF.IntraIntent; time) where R<:ConnectivityIntent
    # implement (Not)GoThroughConstraints
    conint = getintent(idagnode)

    mlg = IBNSims.mlnodegraphtomlgraph(ibn, getrate(conint))

    sourcemlg, initiate_mlg = getmlgsrc(ibn, mlg, conint)
    destmlg = getmlgdst(ibn, mlg, conint)
    gothroughmlg, notgothroughmlg = getgothroughsmlg(mlg, conint)

    mfc = IBNSims.computenondominatedpaths(mlg, sourcemlg, destmlg, getrate(conint); 
                                           gothrough=gothroughmlg, notgothrough=notgothroughmlg, initiate_pathcost=initiate_mlg)

#    !isnothing(initiate_mlg) && @show(mfc)
    mfc1 = optimizenondominatedpaths(ibn, mlg, mfc)

    groomedenableddagexpansion!(ibn, idagnode, mlg, mfc1; time)
    return MINDF.getstate(idagnode)
end

function groomedenableddagexpansion!(ibn::IBN, idagnode::IntentDAGNode{R}, mlg, mfc1::PathCostVector; time) where R<:Intent
    dag = getintentdag(ibn)
    lps, uuids, lightpathtypes = getseparatelightpaths(mlg, mfc1)
    i = 0
    for (lp,uuid,lptype) in zip(lps,uuids,lightpathtypes)
        if !ismissing(uuid)
            add_edge!(dag, getid(idagnode), uuid, nothing)
            continue
        end
        i += 1
        lpu = unique(getindex.(getvertexattr.([mlg], lp), 2))
        lpint = MINDF.getcompliantintent(ibn, getintent(idagnode), LightpathIntent, lpu, mfc1.chosentransmodls[i], lptype)
        isnothing(lpint) && error("Could not create a LightpathIntent")
        # needed mostly for fault detection, since ports and transponders are infinite
        MINDF.isavailable(ibn, lpint) || error("intent resources are not available")
        lpintnode = addchild!(dag, idagnode.id, lpint)
        for lli in MINDF.lowlevelintents(lpintnode.intent)
            addchild!(dag, lpintnode.id, lli)
        end

        # spectrum allocation
        if lptype in [MINDF.borderinitiatelightpath, MINDF.border2borderlightpath]
            bicidx = findfirst(c-> c isa MINDF.BorderInitiateConstraint, getconstraints(getintent(lpintnode)))
            lpr = MINDF.getreqs(getconstraints(getintent(lpintnode))[bicidx])
            spslots = lpr.spslots
        else
            fs = [get_prop(ibn.ngr, e, :link) for e in edgeify(lpintnode.intent.path)]
            trmdlslots = getfreqslots(lpintnode.intent.transmodl)
            trmdlrate = getrate(lpintnode.intent.transmodl)
            startingslot = MINDF.firstfit(fs, trmdlslots)
            startingslot === nothing && error("Not enough slots for transmission module chosen")
            spslots = startingslot:startingslot+trmdlslots-1
        end
        speint = MINDF.getcompliantintent(ibn, lpintnode.intent, MINDF.SpectrumIntent, lpintnode.intent.path, getrate(getintent(lpintnode)), spslots)
        speint === nothing && error("Could not create a SpectrumAllocationIntent")

        if speint !== nothing && MINDF.isavailable(ibn, speint)
            speintnode = addchild!(dag, lpintnode.id, speint)
            for lli in MINDF.lowlevelintents(speintnode.intent)
                addchild!(dag, speintnode.id, lli)
            end
            MINDF.try2setstate!(speintnode, ibn, Val(MINDF.compiled); time)
            MINDF.try2setstate!(idagnode, ibn, Val(MINDF.compiled); time)
        end
    end

end

"$(TYPEDSIGNATURES) Handles interdomain connectivity intents"
function jointrmsagenerilizeddijkstra!(myibn::IBN, neibn::IBN, idagnode::IntentDAGNode{T}, iid::MINDF.InterIntent{R} ;
                time)  where {T<:ConnectivityIntent, R<:MINDF.IntentDirection}
    dag = getintentdag(myibn)
    iidforward = R <: MINDF.IntentForward
    conint = getintent(idagnode)

    if iidforward
        if getdst(conint) in globalnode.([myibn], bordernodes(myibn; subnetwork_view=false))
            myintent = ConnectivityIntent(getsrc(conint), getdst(conint), getrate(conint),
                                          vcat(MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        else
            myintent = MINDF.DomainConnectivityIntent(getsrc(conint), getid(neibn), getrate(conint),
                                          vcat(MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        end
    else
        if getsrc(conint) in globalnode.([myibn], bordernodes(myibn; subnetwork_view=false))
            myintent = ConnectivityIntent(getdst(conint), getsrc(conint), getrate(conint),
                          vcat(MINDF.ReverseConstraint(), MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        else
            myintent = MINDF.DomainConnectivityIntent(getdst(conint), getid(neibn), getrate(conint),
                                  vcat(MINDF.ReverseConstraint(), MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        end
    end
    domint = addchild!(dag, getid(idagnode), myintent)
    state = MINDF.compile!(myibn,  domint, jointrmsagenerilizeddijkstra!; time)

    # create an intent for fellow ibn
    if state == MINDF.compiled
        globalviewpath = MINDF.getcompiledintentpath(myibn, getid(domint))
        updatedconstraints = MINDF.adjustNpropagate_constraints!(myibn, idagnode)
        transnode = globalviewpath[end]
        lpr = MINDF.getlastlightpathrequirements(myibn, getid(domint))
        initconstr = MINDF.BorderInitiateConstraint(NestedEdge(globalviewpath[end-1:end]...), lpr)
        if iidforward
            remintent = ConnectivityIntent(transnode, getdst(conint), getrate(conint),
                                           vcat(initconstr, updatedconstraints), getconditions(myintent))
        else
            remintent = ConnectivityIntent(transnode, getsrc(conint), getrate(conint),
                           vcat(MINDF.ReverseConstraint() , initconstr, updatedconstraints), getconditions(myintent))
        end
        success = MINDF.delegateintent!(myibn, neibn, idagnode, remintent, jointrmsagenerilizeddijkstra!; time)
    end
    MINDF.try2setstate!(idagnode, myibn, Val(MINDF.compiled); time)
    return getstate(idagnode)
end

"$(TYPEDSIGNATURES)"
function jointrmsagenerilizeddijkstra!(ibn::IBN, idagnode::IntentDAGNode{R}; time) where {R<:MINDF.DomainConnectivityIntent}
    dag = getintentdag(ibn)
    conint = getintent(idagnode)

    srcdsts = getintrasrcdst(ibn, getintent(idagnode))
    constraints = getconstraints(conint)

    mlg = IBNSims.mlnodegraphtomlgraph(ibn, getrate(conint))

    mfcfinals = [ let
        sourcemlg, initiate_mlg = getmlgsrc(ibn, mlg, source, constraints)
        destmlg = getmlgdst(mlg, dest, constraints)
        gothroughmlg, notgothroughmlg = getgothroughsmlg(mlg, constraints)
        mfc = IBNSims.computenondominatedpaths(mlg, sourcemlg, destmlg, getrate(conint); 
                                               gothrough=gothroughmlg, notgothrough=notgothroughmlg, initiate_pathcost=initiate_mlg)
        mfc1 = optimizenondominatedpaths(ibn, mlg, mfc)
    end for (source,dest) in srcdsts]

    # compare between unsimilar endnodes pathcost vectors
    mfcwinner = first(mfcfinals)

    # find out one single mfc1
    groomedenableddagexpansion!(ibn, idagnode, mlg, mfcwinner; time)
    return getstate(idagnode)
end

"$(TYPEDSIGNATURES) Return a collection of valid sources and destinations combinations in the intranet"
function getintrasrcdst(ibn::IBN, intent::MINDF.DomainConnectivityIntent{Tuple{Int,Int}, Int})
    neibnidx = MINDF.getlocalibnindex(ibn, getdst(intent))
    dests = MINDF.nodesofcontroller(ibn, neibnidx)
    [(MINDF.localnode(ibn, MINDF.getsrc(intent); subnetwork_view=false),d) for d in dests]
end

"$(TYPEDSIGNATURES) Return a collection of valid sources and destinations combinations in the intranet"
function getintrasrcdst(ibn::IBN, intent::MINDF.DomainConnectivityIntent{Int, Tuple{Int,Int}})
    neibnidx = MINDF.getlocalibnindex(ibn, getsrc(intent))
    sours = MINDF.nodesofcontroller(ibn, neibnidx)
    [(s, MINDF.localnode(ibn, MINDF.getdst(intent); subnetwork_view=false)) for s in sours]
end

"""
$(TYPEDSIGNATURES) 

Get lightpaths as `Vector{Vector{Int}}`
Lightpaths start/end on routers, except if they go/come from a different domain.
So only the first or last node might be the exceptions.
"""
function getseparatelightpaths(mlg, cv::PathCostVector)
    isrouter(i) = getvertexattr(mlg, i)[1] == routernodetype
    lps = Vector{Vector{Int}}()
    lpuuids = Vector{Union{Missing,UUID}}()
    lptypes = Vector{MINDF.LightpathType}()

    firstsrc = src(cv.path[1])
    lp = [firstsrc]
    for (i,ed) in enumerate(cv.path)
        push!(lp, dst(ed))
        if isrouter(dst(ed))
            push!(lps, lp)
            if isrouter(lp[1])
                push!(lptypes, MINDF.fulllightpath)
            else
                push!(lptypes, MINDF.borderinitiatelightpath)
            end
            lp = Vector{Int}()
            push!(lp, dst(ed))
            push!(lpuuids, getedgeattr(mlg, ed).lpidnids)
        elseif i == length(cv.path) #final iteration lightpath is a HalfLightpath
            push!(lps, lp)
            if isrouter(lp[1])
                push!(lptypes, MINDF.borderterminatelightpath)
            else
                push!(lptypes, MINDF.border2borderlightpath)
            end
            push!(lpuuids, getedgeattr(mlg, ed).lpidnids)
        end
    end
    return lps, lpuuids, lptypes
end

