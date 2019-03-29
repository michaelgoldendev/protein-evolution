module DatasetCreator
	using FastaIO
	using JSON
	using BioStructures

	push!(LOAD_PATH,@__DIR__)
	using Backbone
	using Binaries	
	pdbdir = abspath("../data/pdbs/")

	export insertgaps
	function insertgaps(alignmentseq::AbstractString, vec, gapvalue=-1000.0)
		ret = []
		index = 1
		for c in alignmentseq
			if c == '-'
				push!(ret, gapvalue)
			else
				push!(ret, vec[index])
				index += 1
			end
		end
		return ret
	end

	function get_quality_attributes(pdbname)
		resolution = -1.0
		rvalue = -1.0
		freervalue = -1.0
		pdbfile = downloadpdb(pdbname, pdb_dir=pdbdir);
		fin = open(pdbfile, "r")
		for line in readlines(fin)
			if resolution < 0.0
				m = match(r"^REMARK\s+\d+\s+RESOLUTION.\s+((\d*\.)?\d+)\s+ANGSTROMS.*", line)
				if m != nothing
					resolution = parse(Float64,m[1])
				end
			end
			if rvalue < 0.0
				m = match(r"^REMARK\s+\d+\s+R VALUE.*:\s*((\d*\.)?\d+).*", line)
				if m != nothing
					rvalue = parse(Float64,m[1])
				end
			end 
			if freervalue < 0.0
				m = match(r"^REMARK\s+\d+\s+FREE R VALUE .*:\s*((\d*\.)?\d+).*", line)
				if m != nothing
					freervalue = parse(Float64,m[1])
				end
			end			
		end
		close(fin)
		println("$pdbname ANGSTROM:",resolution)
		println("$pdbname RVALUE:",rvalue)
		println("$pdbname FREER:", freervalue)
		if rvalue < 0.0 || rvalue > 0.25
			println("*****************************************")
		end
		return resolution,rvalue,freervalue
	end

	export fromsequencealignment
	function fromsequencealignment(fastafile,outfile, newickfile=nothing)
		jsondict = Dict{String,Any}()
		if newickfile == nothing
			newickstring, cachefile = Binaries.fasttreeaa(fastafile,branchsupport=false)
			jsondict["newick_tree"] = strip(newickstring)
		else 
			newickstring = open(newickfile) do file
                strip(read(file, String))
            end
            jsondict["newick_tree"] = newickstring
		end

		proteins = []
		FastaIO.FastaReader(fastafile) do fr
		    for (desc, alignedseq) in fr
		    	m = match(r"pdb(....)_(.).*", desc)
		    	usepdb = false
		    	if m != nothing
		    		pdbname = m[1]
		    		resolution,rvalue,freervalue = get_quality_attributes(pdbname)
					if rvalue > 0.0 && rvalue < 0.25 && freervalue < 0.25
						usepdb = true
					else
						println("PDB does not meet criteria")
					end
		    	end

		    	if usepdb 
		    		pdbname = m[1]
			     	chain = m[2]
			     	proteindict = Dict{String,Any}()
			        proteindict["name"] = desc
			        proteindict["pdb"] = pdbname			       
			        proteindict["chain"] = chain
			        resolution,rvalue,freervalue = get_quality_attributes(pdbname)
			        proteindict["resolution"] = resolution
			        proteindict["rvalue"] = rvalue
			        proteindict["freervalue"] = freervalue
			        
			        structure = retrievepdb(pdbname, pdb_dir=pdbdir)

			        polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])

			        if polypeptide["sequence"] != replace(alignedseq, "-" => "")
			        	println("MISMATCH")
			        	println(polypeptide["sequence"])
			        	println(replace(alignedseq, "-" => ""))
			        	println("EXITING")
			        	exit()
			        end

			        unalignedseq = uppercase(replace(alignedseq, "-" => ""))
			        unalignedseq = replace(unalignedseq, r"[^ACDEFGHIKLMNPQRSTVWY]" => "X")
			        
			        proteindict["sequence"] = unalignedseq
			        proteindict["aligned_sequence"] = alignedseq
			        #println(alignedseq)

			        proteindict["bond_angles"] = polypeptide["bond_angles"]
			        proteindict["aligned_bond_angles"] = insertgaps(alignedseq, polypeptide["bond_angles"], (-1000.0,-1000.0,-1000.0))

			        proteindict["bond_lengths"] = polypeptide["bond_lengths"]
			        proteindict["aligned_bond_lengths"] = insertgaps(alignedseq, polypeptide["bond_lengths"], (-1000.0,-1000.0,-1000.0))
			        
			        proteindict["omega"] = polypeptide["omega"]
			        proteindict["aligned_omega"] = insertgaps(alignedseq, polypeptide["omega"], -1000.0)
			       	
			       	proteindict["phi_psi"] =  polypeptide["phi_psi"]
			       	proteindict["aligned_phi_psi"] = insertgaps(alignedseq, polypeptide["phi_psi"], (-1000.0,-1000.0))
			       	push!(proteins, proteindict)
		    	else
		    		unalignedseq = uppercase(replace(alignedseq, "-" => ""))
			        unalignedseq = replace(unalignedseq, r"[^ACDEFGHIKLMNPQRSTVWY]" => "X")

			        proteindict = Dict{String,Any}()
			        proteindict["name"] = desc
			        
			        proteindict["sequence"] = unalignedseq
			        proteindict["aligned_sequence"] = alignedseq

			        proteindict["bond_angles"] = [(-1000.0,-1000.0,-1000.0) for c in unalignedseq]
			        proteindict["aligned_bond_angles"] = [(-1000.0,-1000.0,-1000.0) for c in alignedseq]

			        proteindict["bond_lengths"] = [(-1000.0,-1000.0,-1000.0) for c in unalignedseq]
			        proteindict["aligned_bond_lengths"] = [(-1000.0,-1000.0,-1000.0) for c in alignedseq]
			        
			        proteindict["omega"] = [-1000.0 for c in unalignedseq]
			        proteindict["aligned_omega"] = [-1000.0 for c in alignedseq]
			       	
			       	proteindict["phi_psi"] = [(-1000.0,-1000.0) for c in unalignedseq]
			       	proteindict["aligned_phi_psi"] = [(-1000.0,-1000.0) for c in alignedseq]
			       	push!(proteins, proteindict)
		    	end
		    end
		end
		jsondict["proteins"] = proteins	

		fout = open(outfile,"w")
		JSON.print(fout, jsondict)
		close(fout)
	end

	export fromsinglepdb
	function fromsinglepdb(pdbname, chain, outdir)		
		#println(pdbname,"\t",chain)
		jsondict = Dict{String,Any}()
		desc = string("pdb",pdbname,"_",chain)
		outfile = joinpath(outdir, string(desc,".fam"))
		jsondict["newick_tree"] = "($(desc):0.00000);"

		proteins = []	
	 	proteindict = Dict{String,Any}()
	    proteindict["name"] = desc
	    proteindict["pdb"] = pdbname
	    proteindict["chain"] = chain
	    structure = retrievepdb(pdbname, pdb_dir=pdbdir)

	    polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])

	    unalignedseq =  polypeptide["sequence"]
	    alignedseq =  polypeptide["sequence"]
	    
	    proteindict["sequence"] = polypeptide["sequence"]
	    proteindict["aligned_sequence"] = polypeptide["sequence"]

	    proteindict["bond_angles"] = polypeptide["bond_angles"]
	    proteindict["aligned_bond_angles"] = insertgaps(alignedseq, polypeptide["bond_angles"], (-1000.0,-1000.0,-1000.0))

	    proteindict["bond_lengths"] = polypeptide["bond_lengths"]
	    proteindict["aligned_bond_lengths"] = insertgaps(alignedseq, polypeptide["bond_lengths"], (-1000.0,-1000.0,-1000.0))
	    
	    proteindict["omega"] = polypeptide["omega"]
	    proteindict["aligned_omega"] = insertgaps(alignedseq, polypeptide["omega"], -1000.0)
	   	
	   	proteindict["phi_psi"] =  polypeptide["phi_psi"]
	   	proteindict["aligned_phi_psi"] = insertgaps(alignedseq, polypeptide["phi_psi"], (-1000.0,-1000.0))
	   	push!(proteins, proteindict)

		jsondict["proteins"] = proteins

		resolution,rvalue,freervalue = get_quality_attributes(pdbname)
		if rvalue > 0.0 && rvalue < 0.25 && freervalue < 0.25
			fout = open(outfile,"w")
			JSON.print(fout, jsondict)
			close(fout)
		end
	end
end


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
end

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
 
