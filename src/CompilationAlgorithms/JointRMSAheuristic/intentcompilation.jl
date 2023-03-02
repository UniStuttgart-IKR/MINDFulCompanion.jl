function jointrmsagenerilizeddijkstra!(ibn::IBN, idagnode::IntentDAGNode{R}, ::MINDF.IntraIntent; optimizepaths=legacyoptimize,time) where R<:ConnectivityIntent
    conint = getintent(idagnode)


    mlg = mlnodegraphtomlgraph(ibn, getrate(conint))

    sourcemlg, initiate_mlg = getmlgsrc(ibn, mlg, conint)
    destmlg = getmlgdst(ibn, mlg, conint)
    gothroughmlg, notgothroughmlg = getgothroughsmlg(ibn, mlg, conint)
    bordernodesmlg = getbordernodesmlg(ibn, mlg)

    mfc = computenondominatedpaths(mlg, sourcemlg, destmlg, getrate(conint); 
                   gothrough=gothroughmlg, notgothrough=notgothroughmlg, initiate_pathcost=initiate_mlg, bordernodes=bordernodesmlg)

    mfc1 = optimizepaths(ibn, mlg, mfc)

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
            MINDF.try2setstate!(idagnode, ibn, Val(MINDF.compiled); time)
            continue
        end
        i += 1
        lpu = unique(getindex.(getvertexattr.([mlg], lp), 2))
        lpintnode = MINDF.compile!(ibn, idagnode, LightpathIntent, lpu, mfc1.chosentransmodls[i], lptype)
        getstate(lpintnode) ∈ [MINDF.compiled, MINDF.installed, MINDF.installfailed] && continue

        speintnode = MINDF.compile!(ibn, lpintnode, MINDF.SpectrumIntent, lptype, MINDF.firstfit)
        MINDF.try2setstate!(speintnode, ibn, Val(MINDF.compiled); time)
    end

end

"$(TYPEDSIGNATURES) Handles interdomain connectivity intents"
function jointrmsagenerilizeddijkstra!(myibn::IBN, neibn::IBN, idagnode::IntentDAGNode{T}, iid::MINDF.InterIntent{R} ;
                time)  where {T<:ConnectivityIntent, R<:MINDF.IntentDirection}
    dag = getintentdag(myibn)
    iidforward = R <: MINDF.IntentForward
    conint = getintent(idagnode)

    globalizedbordernodes = filter(gn -> gn[1] == getid(neibn) ,globalnode.([myibn], bordernodes(myibn; subnetwork_view=false)))
    bordergtc = MINDF.getfirst(gtc -> getnode(gtc) in globalizedbordernodes ,filter(c -> c isa GoThroughConstraint, getconstraints(conint)))

    if iidforward
        if getdst(conint) in globalizedbordernodes
            myintent = ConnectivityIntent(getsrc(conint), getdst(conint), getrate(conint),
                                          vcat(MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        elseif !isnothing(bordergtc)
            myintent = ConnectivityIntent(getsrc(conint), getnode(bordergtc), getrate(conint),
                                          vcat(MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        else
            myintent = MINDF.DomainConnectivityIntent(getsrc(conint), getid(neibn), getrate(conint),
                                          vcat(MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        end
    else
        if getsrc(conint) in globalizedbordernodes
            myintent = ConnectivityIntent(getdst(conint), getsrc(conint), getrate(conint),
                          vcat(MINDF.ReverseConstraint(), MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
        elseif !isnothing(bordergtc)
            myintent = ConnectivityIntent(getdst(conint), getnode(bordergtc), getrate(conint),
                                          vcat(MINDF.BorderTerminateConstraint() , getconstraints(conint)), getconditions(conint))
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
function jointrmsagenerilizeddijkstra!(ibn::IBN, idagnode::IntentDAGNode{R}; optimizepaths=legacyoptimize, time) where {R<:MINDF.DomainConnectivityIntent}
    dag = getintentdag(ibn)
    conint = getintent(idagnode)

    srcdsts = MINDF.getintrasrcdst(ibn, getintent(idagnode))
    constraints = getconstraints(conint)

    mlg = mlnodegraphtomlgraph(ibn, getrate(conint))

    mfcfinals = [ let
        sourcemlg, initiate_mlg = getmlgsrc(ibn, mlg, source, constraints)
        destmlg = getmlgdst(mlg, dest, constraints)
        gothroughmlg, notgothroughmlg = getgothroughsmlg(ibn, mlg, constraints)
        bordernodesmlg = getbordernodesmlg(ibn, mlg)
        mfc = computenondominatedpaths(mlg, sourcemlg, destmlg, getrate(conint); 
                       gothrough=gothroughmlg, notgothrough=notgothroughmlg, initiate_pathcost=initiate_mlg, bordernodes=bordernodesmlg)
        mfc1 = optimizepaths(ibn, mlg, mfc)
    end for (source,dest) in srcdsts]

    # compare between unsimilar endnodes pathcost vectors
    mfcwinner = optimizediversepaths(ibn, mlg, mfcfinals)

    # find out one single mfc1
    groomedenableddagexpansion!(ibn, idagnode, mlg, mfcwinner; time)
    return getstate(idagnode)
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

