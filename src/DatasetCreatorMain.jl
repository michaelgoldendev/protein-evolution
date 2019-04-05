using ArgParse

push!(LOAD_PATH,@__DIR__)
using DatasetCreator

function parse_dataset_creator_commandline()
    settings = ArgParseSettings()
    #settings.prog = prog
    #settings.version = version
    #settings.add_version = true

    @add_arg_table settings begin
        "fastafile"
        	help = "Specified a FASTA sequence alignment in order to create a .fam dataset file."
         	arg_type = String
         	required = true	        
    end
    return parse_args(settings)
end

parsed_args = parse_dataset_creator_commandline()
if parsed_args["fastafile"] != nothing	
	path = parsed_args["fastafile"]
	outpath = string(path,".fam")
	DatasetCreator.fromsequencealignment(path, outpath)
end

#=
else
	path = "../data/diverse_rna_virus_structures/norovirus_capsid_P_selection.fas"
	outpath = string(path,".fam")
	DatasetCreator.fromsequencealignment(path, outpath)
end=#





#=
using FastaIO
sequences = AbstractString[]
names = AbstractString[]
fastafile = "pdbselection.fas"
FastaIO.FastaReader(fastafile) do fr
    for (desc, seq) in fr
    	m = match(r"pdb(....)_(.).*", desc)
        #println(m)        
        DatasetCreator.fromsinglepdb(m[1], m[2], "../data/single_pdbs/")
    end
end=#

#=
using JSON
family_dir = "../data/curated_rna_virus_structures/"
fastafiles = filter(f -> endswith(f,".fam"), readdir(family_dir))
for fastafile in fastafiles
	jsondict = JSON.parse(open("$(family_dir)/$(fastafile)", "r"))
	#println(jsondict)
	println(fastafile)
	for protein in jsondict["proteins"]
		if startswith(protein["name"],"pdb")
			println(length(protein["aligned_phi_psi"]),"\t",protein["name"])
			for (a,phi_psi) in zip(protein["aligned_sequence"], protein["aligned_phi_psi"])
				println(a,"\t",phi_psi)
			end
		end
	end
	println("--------------------------------------------")
end
exit()=#

#=
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\curated_rna\\rabbit_hemorrhagic_disease_virus_VP60.select.fasta"
#fastafile = "homologs.fas.muscle.fas"
#DatasetCreator.fromsequencealignment(fastafile, string(fastafile,".fam"))
family_dir = "../data/curated_rna_virus_structures/"
#family_dir = "../data/diverse_rna_virus_structures/"
#family_dir = "../data/single_pdbs/"

fastafiles = filter(f -> endswith(f,".fas"), readdir(family_dir))
for fastafile in fastafiles
	println(fastafile)
	path = abspath(joinpath(family_dir,fastafile))
	outpath = joinpath(family_dir,string(fastafile,".fam"))
	DatasetCreator.fromsequencealignment(path, outpath)
end=#

 #DatasetCreator.fromsinglepdb("5msq", "A", "../data/single_pdbs/")
 #DatasetCreator.fromsinglepdb("9nse", "A", "../data/single_pdbs/")

#DatasetCreator.fromsinglepdb("6h4n", "y", "../data/single_pdbs/")
#DatasetCreator.fromsinglepdb("5wjz", "R", "../data/single_pdbs/")
#DatasetCreator.fromsinglepdb("4aaa", "D", "../data/single_pdbs/")
#DatasetCreator.fromsinglepdb("5v74", "35", "../data/single_pdbs/")
#DatasetCreator.fromsinglepdb("5i4l", "e1", "../data/single_pdbs/")

#=
DatasetCreator.fromsinglepdb("6b10", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3CFR", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3CF4", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("1thm", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("1thg", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("1TH2", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3pxs", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3ZMU", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5E9U", "G", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4y7x", "O", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5yr2", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5ggf", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5gep", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5gas", "N", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4aay", "G", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4aah", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4aa0", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4a9w", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3puz", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3pus", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3puo", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3puh", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2der", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("1v94", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("1v8z", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("1V8L", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("1v8b", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4uvb", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5p9t", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5p9f", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5p8l", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5f8c", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4xhb", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zyk", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zyg", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zy1", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zxw", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zxq", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5xvx", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6EGK", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5i5d", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("5i3v", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4qby", "E", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4q99", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4q8s", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4q7k", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4q79", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4q6m", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4eis", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4ei2", "J", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4egy", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4egc", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4eff", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4eed", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4edl", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4eah", "E", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4e9v", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4e5t", "H", "../data/single_pdbs/")


DatasetCreator.fromsinglepdb("4e4l", "E", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4e45", "O", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("4e3q", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6hvx", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6hw1", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6hwh", "M", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6hwz", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6hxf", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6i04", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6i07", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6i18", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("6i53", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3ktw", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kti", "G", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3ktb", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kt6", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3ks8", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3krw", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3krf", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kre", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kr7", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kr4", "I", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kqr", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kqg", "D", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kpz", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kpt", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("3kpr", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zua", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2ztm", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zt9", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zsu", "C", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zsd", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zrz", "B", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zru", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zrs", "E", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zr2", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zqm", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zpq", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zp0", "T", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zos", "A", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2zo3", "H", "../data/single_pdbs/")
DatasetCreator.fromsinglepdb("2znn", "A", "../data/single_pdbs/")=#

#DatasetCreator.fromsinglepdb("2zmw", "A", "../data/single_pdbs/")
 
