### Code for running the testing pipeline on edge-critical graphs

include("config.jl")

import Pkg
Pkg.activate(koala_path)

using Koala
using LinearAlgebra, Random
BLAS.set_num_threads(1)

function test_graph(G::Matrix{Bool}, k::Int, times::Vector{Float64}, NPAs::Vector{Union{Koala.NPA, Nothing}})::Int
	test_id = 1
	times .= 0
	game, distribution = coloring_game(G, k-1)
	filterings = [same_output_reduced, same_output, probabilistic, full];

	while(test_id <= 5)
		success = false
	
		if(test_id == 1)
			times[1] = @elapsed success = xi_SDP(G; target_val = k, max_iter = 250) == k
		else
			if(NPAs[test_id-1] == nothing)
				NPAs[test_id-1] = NPASynchronous(game, 2; filtering = filterings[test_id-1])
			end
			bound = 1.0
			if(test_id >= 4)
				times[test_id] = @elapsed bound = upper_bound_game(game, distribution, NPAs[test_id-1]; verbose = false, formal = true, rho = 1e-02, target = 1 - 2e-04, kkt_solver = Koala.COSMO.
				with_options(Koala.COSMO.CGIndirectKKTSolverGPU))
			else
				times[test_id] = @elapsed bound = upper_bound_game(game, distribution, NPAs[test_id-1]; verbose = false, max_iter = 1000, formal = true, target = 1 - 2e-04, rho = 1e-02)
			end
			success = bound < 1 - 1e-04
		end
		
		if(success)
			return test_id
		else
			test_id += 1
		end
	end
	return test_id
end
					
function test_graphs(n::Int, k::Int, infile::String, outfile::String, thread_id::Int, N::Int)
	println(infile)
	@assert thread_id >= 1 && thread_id <= N
	NPAs::Vector{Union{Koala.NPA, Nothing}} = [nothing, nothing, nothing, nothing]
	times = [0.0, 0.0, 0.0, 0.0, 0.0]
	already_seen = []
	if(isfile(outfile))
		for f in readlines(outfile)
			push!(already_seen, hash(split(f)[1]))
		end
	end
	already_seen = Set(already_seen)
	out = open(outfile, "a+")
	
	for graph6_rep in shuffle(readlines(infile))
		if((hash(graph6_rep) - thread_id) % N == 0 && !(hash(graph6_rep) in already_seen))
			G = parse_graph6(graph6_rep)
			result = test_graph(G, k, times, NPAs)
			write(out, graph6_rep, " | ", join(times, " "), " | $(result)\n")
			if(times[2] != 0.0)
				flush(out)
				if(rand() < 0.3)
					GC.gc(true)
				end
			end
			#catch end
		end
	end
	close(out)
end
	
thread_id = parse(Int, ARGS[1])
for k=4:14
	for n=k:14
		if(n == 13 && k == 5) continue end
		if(isfile("$(critical_graphs_path)/output$(n)_$(k).g6"))
			test_graphs(n, k, "$(critical_graphs_path)/output$(n)_$(k).g6", "$(processed_critical_graphs_path)/output$(n)_$(k)_$(thread_id).g6", thread_id, N_threads)
		end
	end
end

