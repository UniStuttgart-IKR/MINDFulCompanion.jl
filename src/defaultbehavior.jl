#=
# Default data taken from:
# F. Christou, T. Enderle and A. Witt, 
# "Towards a Hybrid Architecture by Introducing Coherent Pluggable Transceivers in IP-Optical Core Networks with Optical Cross-Connects," 
# Photonic Networks; 23th ITG-Symposium, Berlin, Germany, 2022, pp. 1-8.
=#

defaultlinecards() = [MINDF.LineCardDummy(10, 100, 26.72), MINDF.LineCardDummy(2, 400, 29.36), MINDF.LineCardDummy(1, 1000, 31.99)]
defaultlinecardchassis() = [MINDF.LineCardChassisDummy(Vector{MINDF.LineCardDummy}(), 4.7, 16)]

defaulttransmissionmodules() = [MINDF.TransmissionModuleView("DummyFlexibleTransponder",
            MINDF.TransmissionModuleDummy([MINDF.TransmissionProps(5080.0u"km", 300, 8),
            MINDF.TransmissionProps(4400.0u"km", 400, 8),
            MINDF.TransmissionProps(2800.0u"km", 500, 8),
            MINDF.TransmissionProps(1200.0u"km", 600, 8),
            MINDF.TransmissionProps(700.0u"km", 700, 10),
            MINDF.TransmissionProps(400.0u"km", 800, 10)],0,20)),
                                            MINDF.TransmissionModuleView("DummyFlexiblePluggables",
            MINDF.TransmissionModuleDummy([MINDF.TransmissionProps(5840.0u"km", 100, 4),
            MINDF.TransmissionProps(2880.0u"km", 200, 6),
            MINDF.TransmissionProps(1600.0u"km", 300, 6),
            MINDF.TransmissionProps(480.0u"km", 400, 6)],0,8))
           ]
