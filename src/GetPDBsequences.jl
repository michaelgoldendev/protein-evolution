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
getpdbsequencealignment(pdbs, parsed_args["outfile"])

