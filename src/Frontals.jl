module Frontals

#=
https://docs.julialang.org/en/v1/stdlib/Markdown/#markdown_stdlib
https://docs.julialang.org/en/v1/manual/documentation/
https://www.juliapackages.com/p/combinatorics <- partitions of an integer
https://docs.julialang.org/en/v1/manual/performance-tips/
https://oscar-system.github.io/Oscar.jl/stable
https://oscar-system.github.io/Oscar.jl/stable/AbstractAlgebra/ytabs/#Partition

=#

using Oscar

# Untested!
function jacob(I::Vector{fmpz_mpoly})
    # https://docs.julialang.org/en/v1/manual/functions/#Argument-type-declarations
    
    
    R = parent(I[1])
    n,p = length(gens(R)), length(I)
    
    MR = MatrixSpace(R, n, p)
    JI = zero(MR)
    
    for i = 1:p
        for j = 1:n
            JI[i,j] = derivative(I[i], R[j])
        end
    end
    
    return ideal(R, minors(J, min(n,p)))
    # https://oscar-system.github.io/Oscar.jl/stable/AbstractAlgebra/matrix/#minors-Tuple%7BMatElem,%20Int64%7D
end


	#=
	S, all_vars = PolynomialRing(coefficient_ring(R), vcat(symbols(R), list_of_new_vars));
	this creates a polynomial ring not over R, but over the coefficient ring of R with a flat hierarchy of all variables.
	This is probably what you want.
	In this case you have to set up the homomorphism manually via hom(R, S, gens(S)[1:ngens(R)])
	=#
	

    function prenormal(f)
    #=============================================================================
    static proc c1fsurf (map f)
    {
        /* Given a frontal f[x,y] = [x, p[x,y], q[x,y]], either p_y | q_y or
         * q_y | p_y.
         * Assuming q_y = l * p_y, the singular set S of the Nash lifting is given
         * by the zero locus of p_y and lambda_y.
         * This procedure assumes f is given as above to find x, y, gcd(p_y, q_y)
         * and the codimension of S. */

        is_frontal(f);
        string source_name = nameof(basering);
        string target_name = nameof(preimage(f));

        proc find_y (map f)
        {
            if (f[1] == var(1)) { return(var(1), var(2)); }
            if (f[1] == var(2)) { return(var(2), var(1)); }

            "// ** please introduce a map in the form";
            "// **   map f = " + target_name + ", var, p, q";
            "// ** where var is a variable in " + source_name +
            " and p, q are polynomials";

            ERROR("map not set to prenormal form");
        }

        list xy = find_y(f);          // (x,y)
        poly f2y = diff(f[2], xy[2]);
        poly f3y = diff(f[3], xy[2]);

        // return x, p, q, y, codim(Sigma) (in this order)
        if (reduce(f3y, f2y) == 0)
        {
            poly ly = diff(f3y / f2y, y);
            ideal S = f2y, ly;
            return(xy[1], f[2], f[3], xy[2], vdim(std(S)));
        }
        else
        {
            poly ly = diff(f2y / f3y, y);
            ideal S = f3y, ly;
            return(xy[1], f[3], f[2], xy[2], vdim(std(S)));
        }
    }
    ============================================================================#
    end

    function dpscheme
    #=========================================================================
        static proc dp_scheme(poly p, poly q, poly y)
{
    // add an auxiliar variable for divided differences
    ring src = basering;
    int charsrc = ringlist(src)[1];
    ring aux = charsrc, @z, ds;
    def divdifring = src + aux;     // divided differences ring
    setring divdifring;

    poly p = imap(src, p);
    poly q = imap(src, q);
    poly y = imap(src, y);

    "// Computing double point scheme...";
    poly P = (p - subst(p, y, @z)) / (y - @z);
    poly Q = (q - subst(q, y, @z)) / (y - @z);
    poly py = diff(p, y);
    poly tau = resultant(P, Q, @z) / py^2;

    poly alpha = reduce(P, py) / (y - @z);    // 2S + K
    poly beta = reduce(Q, py) / (y - @z);
    ideal PAA = py, alpha, beta;
    int dimPAA = vdim(std(PAA));
 
    setring src;
    poly @cuspedge = imap(divdifring, py);    // S
    ideal P3 = @cuspedge, diff(@cuspedge, y);
    int dimP3 = vdim(std(P3));
 
    poly @dpoints = imap(divdifring, tau);    // 2S + K + W
    ideal PT = @cuspedge, @dpoints;
    int dimPT = vdim(std(PT));

    export(@cuspedge, @dpoints);
    return(dimP3,                             // S
           dimPAA - 2 * dimP3,                // K
           dimPT - dimPAA,                    // W
           milnor(@cuspedge),                 // C
           milnor(@dpoints));                 // D
}
        ====#
    end
    
end # module
