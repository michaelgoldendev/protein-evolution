using ArgParse
include("PDBsearch.jl")

function find_random_pdb_families(t::Int)
	rng = MersenneTwister(4917104813143)
	countsequences = 0
	indices = Int[]
	pdbdatabase = "../data/pdbsequences.fasta"
	open(pdbdatabase) do file
		for (index, ln) in enumerate(eachline(file))
	    	if startswith(ln,">")
	    		push!(indices, index)
	    	end
	    end
	end
	shuffle!(rng, indices)

	numthreads = 16
	pdbdatabasecopy = "../data/pdbsequences$(t).fasta"
	if isfile(pdbdatabasecopy)
		rm(pdbdatabasecopy)
	end
	cp(pdbdatabase, pdbdatabasecopy)
 	find_random_pdb_families_helper(indices[t:numthreads:end], pdbdatabasecopy)
end

#find_random_pdb_families()

function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table settings begin
    	"thread"
    	help = ""
     	arg_type = Int
     	required = true
 	end
    return parse_args(settings)
end
find_random_pdb_families(parse_commandline()["thread"])