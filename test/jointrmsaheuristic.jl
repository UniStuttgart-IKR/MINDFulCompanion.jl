
@testset "jointrmsaheuristic.jl" begin
    globalnet = open(testdir*"/topologies/4nets.graphml") do io
        loadgraph(io, "global-network", GraphIO.GraphML.GraphMLFormat(), NestedGraphs.NestedGraphFormat())
    end
    simgraph = MINDFul.simgraph(globalnet; 
                                router_lcpool=MINDFC.defaultlinecards(), 
                                router_lccpool=MINDFC.defaultlinecardchassis(), 
                                router_lcccap=3,
                                transponderset=MINDFC.defaulttransmissionmodules())
    myibns = MINDFul.nestedGraph2IBNs!(simgraph)

    conint = ConnectivityIntent((myibns[1].id,3), (myibns[2].id, 5), 40.0)
    intentid = addintent!(myibns[1], conint)

    deploy!(myibns[1], intentid, MINDF.docompile, MINDF.SimpleIBNModus(), jointrmsagenerilizeddijkstra!; time=nexttime())
    deploy!(myibns[1], intentid, MINDF.doinstall, MINDF.SimpleIBNModus(), MINDF.directinstall!; time=nexttime());
    # check if satisfies
    @test MINDF.issatisfied(myibns[1], intentid)
    @test MINDF.anyreservations(myibns[1])

end
