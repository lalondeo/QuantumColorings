include("config.jl")

import Pkg
Pkg.activate(koala_path)
using Koala

using JLD2
using LinearAlgebra
include("no_orthogonal_representation.jl")

"""
	are_orthogonal(clump1::Vector{Vector{Vector{T}}}, clump2::Vector{Vector{Vector{T}}})::Bool

Checks whether the two clumps are orthogonal. """
function are_orthogonal(clump1::Vector{Vector{Vector{T}}}, clump2::Vector{Vector{Vector{T}}})::Bool where {T}
	m = length(clump1)
	n = length(clump1[1])
	for j1=1:m
		for j2=1:m
			tot = 0
			for i=1:n
				tot += dot(clump1[j1][i], clump2[j2][i])
			end
			if(abs(tot) > 1e-10) # In our case, T is Int, so this amounts to checking that tot is exactly zero.
				return false
			end
		end
	end
	return true
end


"""
	is_projective_measurement(E::Vector{Matrix{ComplexF64}})::Bool
	
Checks that E is a projective measurement. """
function is_projective_measurement(E::Vector{Matrix{ComplexF64}})::Bool
	n = size(E[1], 1)
	ok = norm(sum(E) - I) < 1e-10 # Sums to the identity
	for e in E
		ok &= norm(e^2 - e) < 1e-10 && norm(e-adjoint(e)) < 1e-10
	end
	return ok
end
		
	
"""
	build_measurement(clump::Vector{Vector{Vector{T}}})::Vector{Matrix{ComplexF64}}

Builds the projective measurement corresponding to the given clump. """
function build_measurement(clump::Vector{Vector{Vector{T}}})::Vector{Matrix{ComplexF64}} where {T}
	projs = Matrix{ComplexF64}[];
	r = length(clump)
	k = length(clump[1])
	E = zeros(ComplexF64, r*k^2, r*k^2)
	v = zeros(ComplexF64, r*k^2)
	zeta = exp(2*pi*im/k)
	delta = (i,n) -> [j == i ? 1.0 : 0.0 for j=1:n];
	
	for c1=1:k
		for c2=1:k
			E .*= 0
			for i=1:r
				v .*= 0
				for j=1:k					
					w = clump[i][((j+c1) % k == 0) ? k : ((j + c1) % k)]
					v .+= zeta^(c2 * j) * kron(delta(j, k), w / norm(w))
				end
				E .+= v * adjoint(v) / k
			end
			push!(projs, deepcopy(E))
		end
	end
	return projs
end
				
clumps = load_object("clumps.jld2")
measurements = broadcast(build_measurement, clumps)

println("Checking that the quantum 4-colouring of G_{21} is valid...")
@assert !(false in broadcast(is_projective_measurement, measurements)) # Sanity check: make sure that the operators returned by build_measurement are projective measurements

G21 = zeros(Bool, 21, 21)
for i=1:21
	for j=i+1:21
		G21[i,j] = G21[j,i] = are_orthogonal(clumps[i], clumps[j])
		if(G21[i,j]) # Sanity check: make sure that the measurements corresponding to the same color are indeed orthogonal if i and j are adjacent
			@assert sum(norm(measurements[i][k] * measurements[j][k]) for k=1:4) < 1e-10 
		end
	end
end
println("Checking that G_{21}'s chromatic number is indeed 5, as a sanity check...")
@assert compute_chromatic_number(G21) == 5

println("Proving that G21 has no 4-dimensional orthogonal representation...")
total_time = @elapsed ok, tot = prove_no_representation_exists(G21,4) # ok is whether the proof succeeded, tot is the number of evaluations of xi_SDP
@assert ok
println("Proof complete, evaluated xi_SDP $(tot) times, took $(total_time) seconds.")
