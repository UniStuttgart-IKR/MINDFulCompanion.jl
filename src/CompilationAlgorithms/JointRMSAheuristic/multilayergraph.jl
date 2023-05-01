"""
$(TYPEDSIGNATURES) 

Convert a nested graph composed by `MLNode`s to a multi layer nested graph.
New multi layer graph loses all info about nested stuctures
Return also a data structure to nodes of the new nested graph to the nodes of the old one.

There are 2 layers in the multi layer graph: PHY, IP
"""
function mlnodegraphtomlgraph(ngr::NestedGraph)
    mlg = NestedGraph([OAttributeGraph(MultiDiGraph(); vertex_type=Tuple{NodeType,Int}, edge_type=LinkCostVector), OAttributeGraph(MultiDiGraph(); vertex_type=Tuple{NodeType,Int},edge_type=LinkCostVector)])
    for v in vertices(ngr)
        if MG.has_prop(ngr, v, :mlnode)
            mlnode = MINDF.getmlnode(ngr, v)
            # oxc
            addvertex!(mlg; subgraphs=1)[end]
            oxcnode = mlg.vmap[end]
            addvertexattr!(mlg, length(mlg.vmap), (oxcnodetype, v))
            # router
            addvertex!(mlg; subgraphs=2)[end]
            routernode = mlg.vmap[end]
            addvertexattr!(mlg, length(mlg.vmap), (routernodetype, v))
            # transmission modules
            for (i,tm) in enumerate(mlnode.transmodulespool)
                # virtual to optical
                ne_down = NestedEdge(routernode, oxcnode)
                addedge!(mlg, ne_down)
                lc_down = LinkCostVector(Val(virtualtoopticallink), MINDF.getcost(tm), deepcopy(tm), edge(mlg,ne_down))
                addedgeattr!(mlg, edge(mlg,ne_down), i, lc_down)
                # optical to virtual
                addedge!(mlg, reverse(ne_down))
                lc_up = LinkCostVector(Val(opticaltovirtuallink), MINDF.getcost(tm), edge(mlg, reverse(ne_down)))
                addedgeattr!(mlg, edge(mlg, reverse(ne_down)), i, lc_up)
            end
        else
            oxcnode = addvertex!(mlg; subgraphs=1)[end]
            addvertexattr!(mlg, length(mlg.vmap), (oxcnodetype, v))
        end
    end
    for e in edges(ngr)
        ne = NestedEdge((1,e.src), (1,e.dst))
        addedge!(mlg, ne)
        fv = MINDF.getlink(ngr, e.src, e.dst)
        lcv = LinkCostVector(Val(opticallink), MINDF.getdistance(fv), MINDF.getspectrumslots(fv), edge(mlg, ne))
        addedgeattr!(mlg, edge(mlg,ne), lcv)
    end
    return mlg
end

"$(TYPEDSIGNATURES) Add virtual links and ip cost"
function mlnodegraphtomlgraph(ibn::IBN, rate::Real; bordergroomingenabled=true)
    mlg = mlnodegraphtomlgraph(ibn.ngr)
    # get ip port cost
    for ed in asmultiedges(edges(mlg))
        lc = getedgeattr(mlg, Tuple(ed)...)
        if lc.linktype == virtualtoopticallink || lc.linktype == opticaltovirtuallink
            ibnnode = getvertexattr(mlg, src(ed))[2]
            router = MINDF.getrouter(ibn, ibnnode)
            lc.costip = MINDF.getportusagecost(router, rate)
        end
    end
    # do virtual links
    lps, lpuuids, rescap, lptypes, spectrums = findlightpathsNusage(ibn)
    for (lp, lpuuid, resc, lptype, spslots) in zip(lps, lpuuids, rescap, lptypes, spectrums)
        !bordergroomingenabled && lptype in [MINDF.borderinitiatelightpath, MINDF.borderterminatelightpath ,MINDF.border2borderlightpath] && continue
        if resc >= rate
            if lptype in [MINDF.borderinitiatelightpath, MINDF.border2borderlightpath]
                mlgvsrc = findfirst( x -> x==(oxcnodetype, lp[1]), vertex_attr(mlg))
            else
                mlgvsrc = findfirst( x -> x==(routernodetype, lp[1]), vertex_attr(mlg))
            end

            if lptype in [MINDF.borderterminatelightpath, MINDF.border2borderlightpath]
                mlgvdst = findfirst( x -> x==(oxcnodetype, lp[end]), vertex_attr(mlg))
            else
                mlgvdst = findfirst( x -> x==(routernodetype, lp[end]), vertex_attr(mlg))
            end
            mlgvbetween = [findfirst( x -> x==(oxcnodetype, lp[i]), vertex_attr(mlg)) for i in 2:length(lp)-1]
            mlgpath = vcat(mlgvsrc, mlgvbetween, mlgvdst)
            multipl = multiplicity(AttributeGraphs.getgraph(mlg.flatgr), mlgvsrc, mlgvdst)
            lc_virtual = LinkCostVector(Val(virtuallink), edgeify(mlgpath), spslots, lpuuid)
            add_edge!(mlg, mlgvsrc, mlgvdst)
            addedgeattr!(mlg, mlgvsrc, mlgvdst, multipl+1, lc_virtual)
        end
    end
    return mlg
end

"$(TYPEDSIGNATURES) Return multilayer nodes in Vector{Vector{Int}} format. Used for plotting"
function getmultilayernodegroups(mlg)
    initializeupuntil(mlns, v) = v > length(mlns) && push!(mlns, fill(Int[], v-length(mlns))...)
    mlnodes = Vector{Vector{Int}}()
    for (i,c) in enumerate(getfield.(vertex_attr(mlg),2))
        initializeupuntil(mlnodes, c)
        push!(mlnodes[c], i)
    end
    return mlnodes
end

"$(TYPEDSIGNATURES) Return `Vector{Vector{Int}}` of available lightpaths and `Vector{Float}` of remaining rate"
function findlightpathsNusage(ibn::IBN)
    lps = Vector{Vector{Int}}()
    lpuuids = Vector{UUID}()
    restratevec = Vector{Float64}()
    lightpathtypes = Vector{MINDF.LightpathType}()
    spectrums = Vector{UnitRange{Int}}()
    for lpi in filter(x -> getintent(x) isa LightpathIntent, getallintentnodes(ibn))
        getstate(lpi) == MINDF.installed || continue
        # also add UUID
        if any(c -> c isa MINDF.BorderInitiateConstraint, getconstraints(getintent(lpi)))
            if any(c -> c isa MINDF.BorderTerminateConstraint, getconstraints(getintent(lpi)))
                push!(lightpathtypes, MINDF.border2borderlightpath)
            else
                push!(lightpathtypes, MINDF.borderinitiatelightpath)
            end
        elseif any(c -> c isa MINDF.BorderTerminateConstraint, getconstraints(getintent(lpi)))
            if any(c -> c isa MINDF.BorderInitiateConstraint, getconstraints(getintent(lpi)))
                push!(lightpathtypes, MINDF.border2borderlightpath)
            else
                push!(lightpathtypes, MINDF.borderterminatelightpath)
            end
        else
            push!(lightpathtypes, MINDF.fulllightpath)
        end
        push!(lps, MINDF.getpath(getintent(lpi)))
        push!(lpuuids, getid(lpi))
        usedrate = getlightpathusedrate(getintentdag(ibn), lpi)
        restrate = getrate(MINDF.gettransmodl(getintent(lpi))) - usedrate
        push!(restratevec, restrate)
        speint = getintent(getintentnodedescendant(ibn, lpi, MINDF.SpectrumIntent))
        push!(spectrums, speint.spectrumalloc)
    end
    return lps, lpuuids, restratevec, lightpathtypes, spectrums
end

function getintentnodedescendant(ibn::IBN, idn::IntentDAGNode, typeget::Type{R}) where R
    wantedintentnodes = filter(ind -> getintent(ind) isa R, MINDF.descendants(getintentdag(ibn), idn))
    first(wantedintentnodes)
end

getmlgnode(mlg, v, nodetype=oxcnodetype) = findfirst(x -> x[1] == nodetype && x[2] == v, vertex_attr(mlg))
getflatnode(mlg, v) = getvertexattr(mlg, v)[2]
getflatedge(mlg, ed) = Edge( getflatnode(mlg, src(ed)), getflatnode(mlg, dst(ed)) )

function getlightpathusedrate(dag, lpidagn)
    pars = MINDF.parents(dag, lpidagn)
    sum([getrate(getintent(par)) for par in pars])
end

function getmlgsrc(ibn, mlg, conint::ConnectivityIntent; iuuid::UUID=UUID(0x0))
    source = MINDF.localnode(ibn, getsrc(conint); subnetwork_view=false)
    constraints = getconstraints(conint)
    getmlgsrc(ibn, mlg, source, constraints; iuuid)
end

function getmlgsrc(ibn, mlg, source, constraints; iuuid::UUID=UUID(0x0))
#function getmlgsrc(ibn, mlg, source, constraints)
    bic = MINDF.getfirst(c -> c isa MINDF.BorderInitiateConstraint, constraints)
    if isnothing(bic)
        (getmlgnode(mlg, source, routernodetype), nothing)
    else
        # todo: grooming initiate border intent
        edg = Edge(MINDF.localnode(ibn, src(bic.edg); subnetwork_view=false), MINDF.localnode(ibn, dst(bic.edg); subnetwork_view=false))
        edgspslots = MINDF.getspectrumslots(MINDF.getlink(ibn, edg))
        mlgedge = Edge(getmlgnode(mlg, src(edg), oxcnodetype), getmlgnode(mlg, dst(edg), oxcnodetype))

        initialedg = searchforavailablelightpath(mlg, mlgedge, MINDF.getreqs(bic); iuuid=(getid(ibn), iuuid))
        if isnothing(initialedg) 
            pcv = PathCostVector(MINDF.getreqs(bic), mlgedge, edgspslots)
        else # otherwise point to the appropriate lightpath
            pcv = PathCostVector(mlg, initialedg)
        end
        (findfirst(x -> x[1] == oxcnodetype && x[2] == source, vertex_attr(mlg)), pcv)
    end
end

"$(TYPEDSIGNATURES) Search if there is a lightpath created for `bic`. If yes, return the half-lightpath edge. Else return `nothing`"
function searchforavailablelightpath(mlg, edg::AbstractEdge, lpr::MINDF.LightpathRequirements; iuuid=UUID(0x0))
    edglcvs = filter(edge_attr(mlg)) do (_,lcv)
#        edg in lcv.phypath && all([isequal(i in lpr.spslots, b) for (i,b) in enumerate(lcv.spectrum)])
        isequal(lcv.linktype, virtuallink) && mlgedginphypath(mlg, edg, lcv.phypath; iuuid) && all([isequal(i in lpr.spslots, b) for (i,b) in enumerate(lcv.spectrum)])
    end

    length(edglcvs) > 1 && @warn("Found more than one appropriate lightpath to map to")
    length(edglcvs) == 0 && return nothing
    return SingleMultiEdge(first(edglcvs)[1])
end

function mlgedginphypath(mlg, edg, phpath; iuuid=UUID(0x0))
    fledg = getflatedge(mlg, edg)
    for phedg in phpath
        flphedg = getflatedge(mlg, phedg)
        fledg == flphedg && return true
    end
    return false
end

function getmlgdst(ibn, mlg, conint::ConnectivityIntent)
    dest = MINDF.localnode(ibn, getdst(conint); subnetwork_view=false)
    constraints = getconstraints(conint)
    res = getmlgdst(mlg, dest, constraints)
    return res
end

function getmlgdst(mlg, dest, constraints)
    bic = MINDF.getfirst(c -> c isa MINDF.BorderTerminateConstraint, constraints)
    if isnothing(bic)
        getmlgnode(mlg, dest, routernodetype)
    else
        getmlgnode(mlg, dest, oxcnodetype)
    end
end

function getbordernodesmlg(ibn, mlg)
    [getmlgnode(mlg, bd, oxcnodetype) for bd in MINDF.bordernodes(ibn; subnetwork_view=false)]
end

getgothroughsmlg(ibn, mlg, conint::ConnectivityIntent) = getgothroughsmlg(ibn, mlg, getconstraints(conint))
function getgothroughsmlg(ibn, mlg, constraints)
    gothroughs = Vector{Int}()
    notgothroughs = Vector{Int}()
     for gtc in filter(c -> c isa GoThroughConstraint, constraints)
         getid(ibn) == getnode(gtc)[1] || continue
         if MINDF.getreqlayer(gtc) == MINDF.signalElectrical
             # must go through router
             push!(gothroughs, getmlgnode(mlg, getnode(gtc)[2], routernodetype))
         elseif MINDF.getreqlayer(gtc) == MINDF.signalUknown
             # must go through OXC (always does)
             push!(gothroughs, getmlgnode(mlg, getnode(gtc)[2], oxcnodetype))
         elseif MINDF.getreqlayer(gtc) == MINDF.signalOXCbypass
             # must NOT go through router
             push!(gothroughs, getmlgnode(mlg, getnode(gtc)[2], oxcnodetype))
             push!(notgothroughs, getmlgnode(mlg, getnode(gtc)[2], routernodetype))
         end
     end
     return (gothroughs, notgothroughs)
end

"""
$(TYPEDSIGNATURES) 

`source` and `dest` must be a router for now.
(In MD case this should not be a restriction)
"""
function computenondominatedpaths(mlg::NestedGraph, source::Int, dest::Int, rate::Float64; gothrough=Int[], notgothrough=Int[], initiate_pathcost=nothing, bordernodes=Int[], minrate=false, iuuid)
    Mf = Vector{PathCostVector}()
    M = Vector{PathCostVector}()

    # initialize M
    if !isnothing(initiate_pathcost)
        source = dst(initiate_pathcost.path[end])
        if source == dest || source ∈ bordernodes
            return [initiate_pathcost]
        end
    end
    for ed in Iterators.filter(e-> src(e)==source , asmultiedges(edges(mlg)))
        dst(ed) in notgothrough && continue
        dst(ed) in bordernodes && dst(ed) != dest && continue
#TODO        validnewedgesfun(mlg, ed, source, pc) || continue
        pc = isnothing(initiate_pathcost) ? let lc = getedgeattr(mlg, Tuple(ed)...);
            PathCostVector(lc, ed)
        end : add(initiate_pathcost, mlg, ed, rate; minrate)
        # do not loop
        dst(ed) ∈ pathify(pc.path)[1:end-1] && continue

        # discard Vpj if odesn follow reqs
        longestblock = MINDF.longestconsecutiveblock(==(1), pc.spectrum)
        getedgeattr(mlg, Tuple(ed)...).linktype in [opticaltovirtuallink, virtuallink] || 
        any(tm -> getrate(tm) >= rate &&  longestblock >= getfreqslots(tm), pc.transmods) || continue

        push!(M, pc)
    end

    
    while length(M) > 0
        p = first(M)
        p in Mf && error("path already in Mf")

        # delete dominated paths in Mf
        inds2delete = Vector{Int}()
        for (i,pj) in enumerate(Mf)
            startend(p) == startend(pj) || continue
            if isdominating(p, pj; gothrough)
                push!(inds2delete, i)
            end
        end
        deleteat!(Mf, inds2delete)
        push!(Mf, p)

        deleteat!(M, 1)
        n = dst(p.path[end])

        # delete dominated paths in M from `p`
        inds2delete = Vector{Int}()
        for (i,pj) in enumerate(M)
            startend(p) == startend(pj) || continue
            if isdominating(p, pj; gothrough)
                push!(inds2delete, i)
            end
        end
        deleteat!(M, inds2delete)


        #expand p
        # for all outgoing edges `en` of `n` that don't make a circle
        for ed in Iterators.filter( e -> validnewedgesfun(mlg, e, n, p), asmultiedges(edges(mlg)))
            dst(ed) in notgothrough && continue
            src(ed) == dest && continue
            dst(ed) in bordernodes && dst(ed) != dest && continue
            pj_new = add(p, mlg, ed, rate; minrate)

            # discard Vpj if odesn follow reqs
            longestblock = MINDF.longestconsecutiveblock(==(1), pj_new.spectrum)
            getedgeattr(mlg, Tuple(ed)...).linktype in [opticaltovirtuallink, virtuallink] || 
            any(tm -> getrate(tm) >= rate &&  longestblock >= getfreqslots(tm), pj_new.transmods) || continue

            # check dominance and add if dominance or go back to another outgoing edge
            breakncontinue = false
            for pj in Mf
                startend(pj) == startend(pj_new) || continue
                if isdominating(pj, pj_new; gothrough)
                    breakncontinue = true
                    break
                end
            end
            breakncontinue && continue
            inds2delete = Vector{Int}()
            for (i,pj) in enumerate(M)
                startend(pj) == startend(pj_new) || continue
                if isdominating(pj, pj_new; gothrough)
                    breakncontinue = true
                    break
                elseif isdominating(pj_new, pj; gothrough)
                    push!(inds2delete, i)
                end
            end
            deleteat!(M, inds2delete)
            breakncontinue && continue
            pj_new in M && error("path already in M")
            if dst(ed) == dest
                finalizepathcostvec!(pj_new; minrate)
            end
            push!(M, pj_new)
        end
    end
    if isnothing(initiate_pathcost)
        filter(m-> startend(m) == (source,dest), Mf)
    else
        filter(m-> startend(m) == (src(initiate_pathcost.phypath[1]),dest), Mf)
    end
end

"""
new edges must
- start from the current path ending
- don't make cycles, except for regeneration
- regeneration can be done once per node
"""
function validnewedgesfun(mlg, e,n,p)
    startfromcurrentpathend = src(e)==n
    # dst shouldn't be in the path.
    # if it is, it should be before 2 nodes (regeneration) and nowhere time before
    noillegalcycles = let
        finalnode = dst(e) 
        currentpath = pathify(p.path)
        onelayerpath = unique(getflatnode.([mlg], currentpath))
#        if getflatnode(mlg, finalnode) ∉ onelayerpath[1:end-1]
        if finalnode ∉ currentpath
            true
        elseif getvertexattr(mlg, finalnode)[1] == oxcnodetype && finalnode == currentpath[end-1] && 
            getvertexattr(mlg, currentpath[end])[1] == routernodetype && finalnode ∉ currentpath[1:end-2]
#            getflatnode(mlg, finalnode) == getflatnode(mlg, currentpath[end])
            true
        else
            false
        end
    end
    startfromcurrentpathend && noillegalcycles
end

"$(TYPEDSIGNATURES) Pick path out of candidate paths. Current implementation picks the one with lower cost"
function legacyoptimize(ibn, mlg, vpcs::Vector{<:PathCostVector}; iuuid=nothing)
    getcostnow(x) = x.costip + x.costopt
    getratenow(x) = sum(getrate.(x.chosentransmodls); init=0.0) 
    getdistnow(x) = MINDF.getdistance(ibn, unique(getindex.(getvertexattr.([mlg], pathify(x.phypath)),2)) )

    minind = findmin(vpc -> (+getcostnow(vpc), +vpc.virtuallinks, -getratenow(vpc), +getdistnow(vpc)), vpcs)[2]
    return vpcs[minind]
end

"$(TYPEDSIGNATURES) Pick path out of candidate paths. Current implementation picks the one with lower cost"
function maximizevirtuallinks(ibn, mlg, vpcs::Vector{<:PathCostVector})
    getcostnow(x) = x.costip + x.costopt
    getratenow(x) = sum(getrate.(x.chosentransmodls); init=0.0) 
    getdistnow(x) = MINDF.getdistance(ibn, unique(getindex.(getvertexattr.([mlg], pathify(x.phypath)),2)) )

    minind = findmin(vpc -> (-vpc.virtuallinks, +getcostnow(vpc), -getratenow(vpc), +getdistnow(vpc)), vpcs)[2]
    return vpcs[minind]
end

"$(TYPEDSIGNATURES) Pick path out of candidate paths. Current implementation picks the one with lower cost"
function minlatencyoptimize(ibn, mlg, vpcs::Vector{<:PathCostVector})
    getcostnow(x) = x.costip + x.costopt
    getratenow(x) = sum(getrate.(x.chosentransmodls); init=0.0) 
    getdistnow(x) = MINDF.getdistance(ibn, unique(getindex.(getvertexattr.([mlg], pathify(x.phypath)),2)) )

    minind = findmin(vpc -> (+getdistnow(vpc), +getcostnow(vpc), +vpc.virtuallinks, -getratenow(vpc)), vpcs)[2]
    return vpcs[minind]
end


"$(TYPEDSIGNATURES) Pick path out of diverse source-destination paths"
function optimizediversepaths(ibn, mlg, vpcs)
    getdistnow(x) = MINDF.getdistance(ibn, unique(getindex.(getvertexattr.([mlg], pathify(x.phypath)),2)) )
    minind = findmin(vpc -> (+getdistnow(vpc)), vpcs)[2]
    vpcs[minind]
end
