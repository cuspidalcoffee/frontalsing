module Frontals

using Oscar

function embed(f, R)
	
end

function double_point_space(p, q; var)
	og_ring = parent(var);
	ext_ring = PolynomialRing(
					coefficient_ring(og_ring),
					vcat(symbols(og_ring), aux_var)
					);

	hom(og_ring,
		ext_ring,
		gens(ext_ring)[1:ngens(og_ring)]
		);
									

	#=
	S, all_vars = PolynomialRing(coefficient_ring(R), vcat(symbols(R), list_of_new_vars));
	this creates a polynomial ring not over R, but over the coefficient ring of R with a flat hierarchy of all variables.
	This is probably what you want.
	In this case you have to set up the homomorphism manually via hom(R, S, gens(S)[1:ngens(R)])
	=#
	
end


end # module
