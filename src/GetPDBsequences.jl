using FastaIO
using BioStructures
using ArgParse

push!(LOAD_PATH,@__DIR__)
using Backbone
using Binaries
using DatasetCreator

function parse_get_pdb_sequences_commandline()
    settings = ArgParseSettings()
    #settings.prog = prog
    #settings.version = version
    #settings.add_version = true

    @add_arg_table settings begin
        "sequencelist"
        	help = "Comma-seperated list of PDB IDs"
         	arg_type = String
         	required = true
    	"outfile"
        	help = "Path to file where the sequences should be written."
         	arg_type = String
         	default = "sequences.fasta"
    end
    return parse_args(settings)
end

parsed_args  = parse_get_pdb_sequences_commandline()
pdbs = split(parsed_args["sequencelist"], ",")

pdbdir = abspath("../data/pdbs/")
fout = open(parsed_args["outfile"], "w")
for p in pdbs
	spl = split(p,"_")
	pdbname = spl[1]
	chain = ""
	if length(spl) > 1
		chain = spl[2]
	end
	structure = retrievepdb(pdbname, pdb_dir=pdbdir)
	longestsequence = ""
	for c in structure[1]
		polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(c)
		sequence = polypeptide["sequence"]
		if longestsequence == "" || length(c) > length(longestsequence)
			longestsequence =  sequence
			chain = chainid(c)
		end
	end
	if chain != ""
		polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])
		longestsequence = polypeptide["sequence"]
		println(fout, ">pdb",pdbname,"_",chain)
		println(fout, longestsequence)
	end
end
close(fout)

musclefile = string(parsed_args["outfile"], ".muscle.fas")
fastastring,cachefile = Binaries.muscle(parsed_args["outfile"])
fout = open(musclefile, "w")
print(fout,fastastring)
close(fout)

if fastastring != ""
	DatasetCreator.fromsequencealignment(musclefile, string(musclefile, ".fam"))
end

