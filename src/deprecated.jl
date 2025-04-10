"""
$(TYPEDSIGNATURES)
"""
function MINDF.compileintent!(ibnf::IBNFramework, idagnode::IntentDAGNode{<:ConnectivityIntent}, uniformrandomalg::UniformRandomCompilation)
    sourceglobalnode = getsourcenode(getintent(idagnode))
    destinationglobalnode = getdestinationnode(getintent(idagnode))

    returncode::Symbol = ReturnCodes.FAIL

    if getibnfid(ibnf) == getibnfid(sourceglobalnode) == getibnfid(destinationglobalnode)
        # intra-domain
        returncode = uniformrandom!(ibnf, idagnode, uniformrandomalg)
        if returncode === ReturnCodes.FAIL_OPTICALREACH_OPTINIT
            # uncompile
            @assert MINDF.uncompileintent!(ibnf, getidagnodeid(idagnode)) 
            # find shortest distance neighbor j

            # get a node in between the shortest paths
            splitglobalnode = getinbetweenpathnode(ibnf, idagnode, uniformrandomalg)
            intent = getintent(idagnode)
            idag = getidag(ibnf)
            firsthalfintent = ConnectivityIntent(sourceglobalnode, splitglobalnode, getrate(intent), getconstraints(intent))
            firsthalfidagnode = addidagnode!(idag, firsthalfintent; parentid = getidagnodeid(idagnode), intentissuer = MachineGenerated())
            returncode = uniformrandom!(ibnf, firsthalfidagnode, uniformrandomalg)
            updateidagnodestates!(ibnf, firsthalfidagnode)
            issuccess(returncode) || return returncode
            
            secondhalfintent = ConnectivityIntent(splitglobalnode, destinationglobalnode, getrate(intent), filter(x -> !(x isa OpticalInitiateConstraint), getconstraints(intent)))
            secondhalfidagnode = addidagnode!(idag, secondhalfintent; parentid = getidagnodeid(idagnode), intentissuer = MachineGenerated())
            returncode = uniformrandom!(ibnf, secondhalfidagnode, uniformrandomalg)
            updateidagnodestates!(ibnf, secondhalfidagnode)
        end
        updateidagnodestates!(ibnf, idagnode)
    elseif getibnfid(ibnf) == getibnfid(sourceglobalnode) && getibnfid(ibnf) !== getibnfid(destinationglobalnode)
        # source intra-domain , destination cross-domain
        # border-node
        if isbordernode(ibnf, destinationglobalnode)
            #TODO-tomorrow
            returncode = kspffcrossdomain!(ibnf, idagnode, uniformrandomalg, destinationglobalnode)
        else
            destinationglobalbordernode =  
            let 
                # randomly pick a border node
                # TODO-tomorrow
                dglobalbordernode = getfirst(shuffle(getrng(uniformrandomalg), getbordernodesasglobal(ibnf))) do globalbordernode
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
            #TODO-tomorrow
            returncode = kspffcrossdomain!(ibnf, idagnode, uniformrandomalg, destinationglobalbordernode)
        end
    end
    return returncode
end
