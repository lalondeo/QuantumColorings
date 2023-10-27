##### Code for launching the proof showing that the edge-critical graphs have equal quantum and classical chromatic numbers #####
include("config.jl")
for i=1:N_threads
	run(`/bin/bash "julia check_graphs.jl $(i) > out$(i).txt"`; wait=false)
end
