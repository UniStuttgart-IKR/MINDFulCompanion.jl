function jointrmsagenerilizeddijkstra!(ibn::IBN, idagnode::IntentDAGNode{R}, ::MINDF.IntraIntent; optimizepaths=legacyoptimize, minrate=false, bordergroomingenabled=true, time) where R<:ConnectivityIntent
    conint = getintent(idagnode)

    mlg = mlnodegraphtomlgraph(ibn, getrate(conint); bordergroomingenabled)

    sourcemlg, initiate_mlg = getmlgsrc(ibn, mlg, conint; iuuid=getid(idagnode))
#    sourcemlg, initiate_mlg = getmlgsrc(ibn, mlg, conint)
    destmlg = getmlgdst(ibn, mlg, conint)
    gothroughmlg, notgothroughmlg = getgothroughsmlg(ibn, mlg, conint)
    bordernodesmlg = getbordernodesmlg(ibn, mlg)

    mfc = computenondominatedpaths(mlg, sourcemlg, destmlg, getrate(conint); 
                   gothrough=gothroughmlg, notgothrough=notgothroughmlg, 
                   initiate_pathcost=initiate_mlg, bordernodes=bordernodesmlg, minrate, iuuid=(getid(ibn),getid(idagnode)))

    if length(mfc) == 0
        @warn("Didn't find any non-dominated paths")
        return MINDF.getstate(idagnode)
    end

    mfc1 = optimizepaths(ibn, mlg, mfc; iuuid=(getid(ibn), getid(idagnode)))

    groomedenableddagexpansion!(ibn, idagnode, mlg, mfc1; time)
    return MINDF.getstate(idagnode)
end

function groomedenableddagexpansion!(ibn::IBN, idagnode::IntentDAGNode{R}, mlg, mfc1::PathCostVector; time) where R<:Intent
    dag = getintentdag(ibn)
    lps, uuids, lightpathtypes = getseparatelightpaths(mlg, mfc1; iuuid=getid(idagnode))
    i = 0
    for (lp,uuid,lptype) in zip(lps,uuids,lightpathtypes)
        if !ismissing(uuid)
            add_edge!(dag, getid(idagnode), uuid, nothing)
            MINDF.syncnodefromdescendants!(idagnode, ibn; time)
            continue
        end
        i += 1
        lpu = unique(getindex.(getvertexattr.([mlg], lp), 2))
        if lptype == MINDF.border2borderlightpath
            lpintnode = MINDF.compile!(ibn, idagnode, LightpathIntent, lpu, mfc1.transmodule, lptype)
        else
            lpintnode = MINDF.compile!(ibn, idagnode, LightpathIntent, lpu, mfc1.chosentransmodls[i], lptype)
        end
        isnothing(lpintnode) && return getstate(idagnode)
        getstate(lpintnode) ∈ [MINDF.compiled, MINDF.installed, MINDF.installfailed] && continue

        speintnode = MINDF.compile!(ibn, lpintnode, MINDF.SpectrumIntent, lptype, MINDF.firstfit)
        isnothing(speintnode) && return getstate(idagnode)
        MINDF.try2setstate!(speintnode, ibn, Val(MINDF.compiled); time)
    end

end

"$(TYPEDSIGNATURES) Handles interdomain connectivity intents"
function jointrmsagenerilizeddijkstra!(myibn::IBN, neibn::IBN, idagnode::IntentDAGNode{T}, iid::MINDF.InterIntent{R} ;
                args...)  where {T<:ConnectivityIntent, R<:MINDF.IntentDirection}
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
#    state = MINDF.compile!(myibn,  domint, jointrmsagenerilizeddijkstra!; time=args[:time])
    state = MINDF.compile!(myibn,  domint, jointrmsagenerilizeddijkstra!; args...)

    # create an intent for fellow ibn
    if state == MINDF.compiled || state == MINDF.installed
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
#        success = MINDF.delegateintent!(myibn, neibn, idagnode, remintent, jointrmsagenerilizeddijkstra!; time=args[:time])
        success = MINDF.delegateintent!(myibn, neibn, idagnode, remintent, jointrmsagenerilizeddijkstra!; args...)
    end
    MINDF.try2setstate!(idagnode, myibn, Val(MINDF.compiled); time=args[:time])
    return getstate(idagnode)
end

"$(TYPEDSIGNATURES)"
function jointrmsagenerilizeddijkstra!(ibn::IBN, idagnode::IntentDAGNode{R}; optimizepaths=legacyoptimize, minrate=false, bordergroomingenabled=true, time) where {R<:MINDF.DomainConnectivityIntent}
    dag = getintentdag(ibn)
    conint = getintent(idagnode)

    srcdsts = MINDF.getintrasrcdst(ibn, getintent(idagnode))
    constraints = getconstraints(conint)

    mlg = mlnodegraphtomlgraph(ibn, getrate(conint); bordergroomingenabled)

    mfcfinals = [ let
        sourcemlg, initiate_mlg = getmlgsrc(ibn, mlg, source, constraints)
        destmlg = getmlgdst(mlg, dest, constraints)
        gothroughmlg, notgothroughmlg = getgothroughsmlg(ibn, mlg, constraints)
        bordernodesmlg = getbordernodesmlg(ibn, mlg)
        mfc = computenondominatedpaths(mlg, sourcemlg, destmlg, getrate(conint); 
                       gothrough=gothroughmlg, notgothrough=notgothroughmlg,
                       initiate_pathcost=initiate_mlg, bordernodes=bordernodesmlg, minrate, iuuid=(getid(ibn),getid(idagnode)))
        length(mfc) == 0 ? missing : optimizepaths(ibn, mlg, mfc)
    end for (source,dest) in srcdsts]

    if all(ismissing, mfcfinals)
        @warn("Didn't find any non-dominated paths $((getid(ibn), getid(idagnode)))")
        return MINDF.getstate(idagnode)
    end
    # compare between unsimilar endnodes pathcost vectors
    mfcwinner = optimizediversepaths(ibn, mlg, skipmissing(mfcfinals))

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
function getseparatelightpaths(mlg, cv::PathCostVector; iuuid::UUID=UUID(0x0))
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
            push!(lpuuids, getedgeattr(mlg, Tuple(ed)...).lpidnids)
        elseif i == length(cv.path) #final iteration lightpath is a HalfLightpath
            push!(lps, lp)
            if isrouter(lp[1])
                push!(lptypes, MINDF.borderterminatelightpath)
            else
                push!(lptypes, MINDF.border2borderlightpath)
            end
            push!(lpuuids, getedgeattr(mlg, Tuple(ed)...).lpidnids)
        end
    end
    return lps, lpuuids, lptypes
end

