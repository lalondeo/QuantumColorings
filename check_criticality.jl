include("config.jl")
import Pkg
Pkg.activate(koala_path)
using Koala

##### Unit tests to make sure that the lists of edge-k-critical graphs are complete #####
# Will check that every graph in the list is indeed critical, and also generate a few edge-critical graphs and make sure that they are in the list


for file in readdir(critical_graphs_path)
	list = []
	n, k = broadcast((x)->parse(Int, x), split(file[7:end-3], "_"))
	println("*** TESTING n=$(n), k=$(k) ***")
	println("Running criticality check")
	for s in readlines("$(critical_graphs_path)/$(file)")
		G = parse_graph6(s)
		@assert size(G,1) == n
		@assert compute_chromatic_number(G) == k
		@assert test_edge_criticality(G)
		push!(list, compute_invariant(G))
	end
	if(n==k) continue end
	list = Set(list)
	println("Generating edge-critical graphs")
	tot = 0
	for i=1:1000000
		G = generate_edge_critical_graph(n,k)
		if(G != nothing)
			@assert compute_invariant(G) in list
			tot += 1
		end
	end
	println("Number of graphs generated: $(tot)")
end
