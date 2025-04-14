module MINDFulCompanion

using Unitful, UnitfulData, UUIDs, DocStringExtensions
using Graphs

# is not used 
import StatsBase: sample

import AttributeGraphs as AG
import MINDFul as MINDF


import Random: MersenneTwister, rand, AbstractRNG, shuffle!, shuffle, randperm

import MINDFul: IBNFramework, IntentDAGNode, ConnectivityIntent, updateidagnodestates!, getibnfid, isbordernode, getfirst, getlocalnode, getbordernodesaslocal, getglobalnode, getsourcenode, getdestinationnode, getbordernodesasglobal, getintent, getibnag, getidag, getnodeview, getrate, getconstraints, LowLevelIntent, OpticalInitiateConstraint, OpticalTerminateConstraint, getspectrumslotsrange, getopticalreach, getpathspectrumavailabilities, gettransmissionmodulecompat, getglobalnode_input, generatelightpathoxcadddropbypassspectrumlli, MachineGenerated, addidagnode!, RouterPortLLI, getavailabletransmissionmoduleviewindex, gettransmissionmoduleviewpool, gettransmissionmode, getspectrumslotsneeded, TransmissionModuleLLI, TransmissionModuleCompatibility, getname, getrouterview, getidagnodeid, OXCAddDropBypassSpectrumLLI, getweights, NodeView, RouterView, OXCView, getreservations, getrouterportindex, getportnumber, getadddropport, getadddropportnumber, getoxcview, KMf, GBPSf, TransmissionModuleView, IntentCompilationAlgorithm, gettransmissionmodes, istransmissionmoduleandmodecompatible, requestcompileintent_init!, getopticalinitiateconstraint, remoteintent!, GlobalNode, getibnfhandler, ReturnCodes, IBNAttributeGraph, issuccess, LocalNode, GlobalNode

include("randompolicy.jl")
include("2migrate2core.jl")

end
