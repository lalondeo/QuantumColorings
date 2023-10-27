using Combinatorics
using Koala
using Random

"""
	identify_vertices(G::Matrix{Bool}, v1::Int, v2::Int)::Matrix{Bool}

Given two distinct vertices v1 and v2 of G, returns the graph in which the two vertices are identified, i.e. a vertex other than v1 and v2 is adjacent to the new vertex if and only if
it is adjacent to one of them. """
function identify_vertices(G::Matrix{Bool}, v1::Int, v2::Int)::Matrix{Bool}
	@assert v1 != v2
	if(v1 > v2)
		v1,v2 = v2,v1
	end
	
	n = size(G,1)
	G2 = zeros(Bool, n-1, n-1)
	for i=1:n-1
		for j=i+1:n-1
			i_ = (i >= v2 ? i + 1 : i)
			j_ = (j >= v2 ? j + 1 : j)
			G2[i, j] = G2[j, i] = G[i_, j_]
		end
	end
	
	for i=1:n-1
		i_ = (i >= v2 ? i + 1 : i)
		G2[v1, i] |= G[v2, i_]
		G2[i, v1] |= G[v2, i_]
	end
	return G2
end

"""
	set_identification_does_something(G::Matrix{Bool}, S::Vector{Int})::Bool
	
Given a graph G and a subset S of vertices of G, returns whether the pseudo-merge operation with respect to S will add at least one edge to G. 
"""
function set_identification_does_something(G::Matrix{Bool}, S::Vector{Int})::Bool
	@assert (length(Set(S)) == length(S)) && minimum(S) >= 1 && maximum(S) <= size(G, 1)
	
	ok = false
	for i=1:size(G,1)
		if(!(i in S))
			ok |= sum(G[i,j] for j in S) == length(S)-1
		end
	end
	return ok
end


"""
	set_identification(G::Matrix{Bool}, S::Vector{Int})::Matrix{Bool}
	
Given a graph G and a subset S of vertices of G, identifies S, i.e. a vertex not in S is made adjacent to all vertices in S in the new graph if it was adjacent to 
all vertices in S save one in the original graph. 
"""
function set_identification(G::Matrix{Bool}, S::Vector{Int})::Matrix{Bool}
	@assert (length(Set(S)) == length(S)) && minimum(S) >= 1 && maximum(S) <= size(G, 1)
	
	G2 = deepcopy(G)
	for i=1:size(G, 1)
		if(!(i in S))
			if(sum(G[i,j] for j in S) == length(S)-1)
				for j in S
					G2[i,j] = G2[j,i] = 1
				end
			end
		end
	end
	return G2
end


"""
	expand_graph(G::Matrix{Bool}, k::Int)::Vector{Matrix{Bool}
	
Given a graph G and a dimension k, attempts to construct a list of graphs such that G has an orthogonal representation in dimension k if and only if one of the graphs in the list does.
If no nontrivial such list could be built, the empty list is returned. 
"""
function expand_graph(G::Matrix{Bool}, k::Int)::Vector{Matrix{Bool}}
	@assert k >= 3
	K = zeros(Bool, k+1, k+1)
	for i=1:k-1
		for j=0:1
			K[i,k+j] = K[k+j, i] = true
		end
	end
	
	S = pick_subgraph(G,K)
	
	if(isempty(S)) # K isn't even a subgraph of G, we give up
		return S
	end
	
	ok = false
	for _=1:10000
		ok = true
		for R in powerset(S[1:k-1])
			if(length(R) >= 3 && !(set_identification_does_something(G, R))) # The expansion does nothing with respect to the vertices in R: mission abort
				S = pick_subgraph(G,K)
				ok = false
				break
			end
		end
		if(ok)
		break
		end
	end
	
	if(!ok) # Couldn't find a suitable copy of K inside G
		return [] 
	end
	
	Gs::Vector{Matrix{Bool}} = []
	
	if(!(G[S[end-1], S[end]]))
		push!(Gs, identify_vertices(G, S[end-1], S[end]))
	end
	
	for R in powerset(S[1:k-1])
		if(length(R) == 2 && !(G[R[1], R[2]])) # Two vertices: if not adjacent, merge
			push!(Gs, identify_vertices(G, R[1], R[2]))
			
		elseif(length(R) >= 3) # More than two vertices: pseudo-merge
			push!(Gs, set_identification(G, R))
		end
	end
	return Gs
end


function _prove_no_representation_exists(G::Matrix{Bool}, k::Int)::Tuple{Bool,Int}
	if(xi_SDP(G; target_val = k+1, formal = true, max_iter = 100) == k+1) # if xi_SDP(G) > k, G cannot have an orthogonal representation in C^k and so we're done
		return (true, 1)
	end
	
	new_graphs = expand_graph(G, k)
	
	if(isempty(new_graphs))
		return (false, 1)
	else
		tot = 0
		for G2 in new_graphs
			ok, _tot =  _prove_no_representation_exists(G2, k)
			tot += _tot
			if(!ok)
				return (false, tot)
			end
		end

		return (true, tot)
	end
end

"""
	prove_no_representation_exists(G::Matrix{Bool}, k::Int)::Bool

Given a graph G and a dimension k, attempts to prove that G does not have an orthogonal representation in C^k. Returns true if this succeeded. If false is returned, nothing can be 
concluded. """
function prove_no_representation_exists(G::Matrix{Bool}, k::Int)::Tuple{Bool,Int}
	Random.seed!(1)
	return _prove_no_representation_exists(G,k)
end



