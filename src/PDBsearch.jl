using FastaIO
using Statistics
using BioAlignments
using BioStructures

push!(LOAD_PATH,@__DIR__)
using Backbone
using Binaries

pdbdir = abspath("../data/pdbs/")

include("KmerAACommon.jl")

function issubsequence(query::String, target::String)
	if length(query) > length(target)
		return false,0
	end
	b = 0
	len = length(target)
	numinserts = 0
	for q in query
		b += 1
		if b > len
			return false, numinserts
		end
		isinsert = false
		while q != target[b]
			isinsert = true
			b += 1
			if b > len
				return false, numinserts
			end
		end

		if isinsert
			numinserts += 1
		end
	end
	return true, numinserts
end

function find_pdb_homologs(fastafile, outfile="homologs.fas")
	kmers2 = Array{Int,1}[]
	kmers3 = Array{Int,1}[]
	sequences = AbstractString[]
	names = AbstractString[]
	FastaIO.FastaReader(fastafile) do fr
	    for (desc, seq) in fr
	        push!(names,desc)
	        degapped = uppercase(replace(seq, "-" => ""))
	        degapped = replace(degapped, r"[^ACDEFGHIKLMNPQRSTVWY]" => "X")
	        push!(sequences, degapped)
	       	push!(kmers2, freqvectorint(degapped, 2))
	       	push!(kmers3, freqvectorint(degapped, 3))
	    end
	end

	fout = open(outfile, "w")
	for (name,seq) in zip(names,sequences)
		println(fout,">",name)
		println(fout,seq)
	end


	k2cutoff = 0.2
	k3cutoff = 0.2
	alignmentscorecutoff = 0.3


	pdbmatches = Dict{AbstractString,AbstractString}()
	open("../data/pdbsequences.fasta") do file
		name = ""
		seq = ""
		count = 0	
		k2counts = 0
		k3counts = 0
	    for ln in eachline(file)
	    	if startswith(ln,">")
	    		name = ln
	    		seq = ""
	    	else
	    		count += 1
	    		if count % 10000 == 0
	    			#println(count)
	    		end
	    		seq = ln

	    		cutoff = 0.5
				uppercutoff = 1.0/cutoff
				seqindex = 1
				tv2 = nothing
				tv3 = nothing
				maxmatch2 = 0.0
				maxmatch3 = 0.0
				maxalignmentscore = 0.0
				k2match = false
				k3match = false
				bestquery = ""
	    		for (query,queryname) in zip(sequences, names)
	    			lengthratio = length(seq)/length(query)

	    			if lengthratio >= cutoff && lengthratio <= uppercutoff
	    				if tv2 == nothing
		    				tv2 = freqvectorint(seq, 2)
		    			end
		    			qv2 = kmers2[seqindex]
		    			k2score = 1.0 - (sum(abs.(tv2 - qv2)) / sum(max.(tv2,qv2)))
		    			if k2score > k2cutoff
		    				k2match = true
		    				if tv3 == nothing
			    				tv3 = freqvectorint(seq, 3)
			    			end
			    			qv3 =  kmers3[seqindex]
			    			k3score = 1.0 - (sum(abs.(tv3 - qv3)) / sum(max.(tv3,qv3)))
			    			if k3score > k3cutoff
			    				k3match = true

			    				res = pairalign(GlobalAlignment(), seq, query, AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1))
			    				aln = res.aln
			    				#println(count_matches(aln),"\t",count_mismatches(aln),"\t",count_aligned(aln),"\t", length(aln))
				    			if k3score > maxmatch3
				    				maxmatch2 = k2score
				    				maxmatch3 = k3score
				    				maxalignmentscore = count_matches(aln)/(count_matches(aln)+count_mismatches(aln))
				    				bestquery = query
				    			end
					        end
		    			end
		    		end
			        seqindex += 1
			    end
			    if k2match
					k2counts += 1				
				end
				if k3match 
					k3counts += 1
				end
				if count % 10000 == 0
					#println("$(k2counts) / $(count) ($(k2counts/count)) and $(k3counts) / $(count) ($(k3counts/count)) ")
				end
			    if maxmatch2 > 0.0
			    	nameandchain = split(name[2:end])[1]
			    	spl2 = split(nameandchain,"_")
			    	pdbname = spl2[1]
			    	chain = spl2[2]
			    	key = string(">",pdbname,"\n",seq)
			    	if !haskey(pdbmatches, key)
				    	println(count,"\t",pdbname,"\t",chain,"\t",maxmatch2, "\t", maxmatch3,"\t",maxalignmentscore)
				    	res = pairalign(GlobalAlignment(), seq, bestquery, AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1))
	    				aln = res.aln
	    				if maxalignmentscore > alignmentscorecutoff
	    					#println(aln)	
	    					try
		    					structure = retrievepdb(pdbname, pdb_dir=pdbdir)
		    					polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])
		    					#println(name)
		    					#println(polypeptide["sequence"])
		    					println(fout, string(">pdb",pdbname,"_",chain))
		    					println(fout, polypeptide["sequence"])
		    				catch
		    					println("ERROR DOWNLOADING ", pdbname)
		    				end
	    				end
				    	pdbmatches[key] = nameandchain
				    end
			    end
	    	end
	       
	    end
	end
	close(fout)


	if length(pdbmatches) > 0
		musclefile = string(outfile,".muscle.fas")
		muscle_alignment, cachefile = Binaries.muscle(outfile)
		fout = open(musclefile,"w")
		println(fout,muscle_alignment)
		close(fout)

		newickstring, cachefile = Binaries.fasttreeaa(musclefile,branchsupport=true)
		fout = open(string(outfile,".nwk"),"w")
		println(fout,newickstring)
		close(fout)
	end
end


#fastafile = abspath("../data/influenza_a/HA/selection3.fasta")
#fastafile = abspath("../data/hiv/curated6.fasta")
#fastafile = abspath("../data/test_data/hiv_pol_selection.fasta")
#fastafile = abspath("../data/test_data/maise_streak_virus_coat_protein_selection.fasta")
#fastafile = abspath("../data/test_data/westnile_dengue_selection.fasta")
#find_pdb_homologs(fastafile)


#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\curated_rna\\avian_coronavirus_NSP3.select.fasta"
#family_directories = ["../data/curated/curated_rna/", "../data/selected_families/"]
family_directories = ["../data/curated_selection"]
for family_dir in family_directories
	fastafiles = filter(f -> endswith(f,".fasta"), readdir(family_dir))
	for fastafile in fastafiles
		path = abspath(joinpath(family_dir,fastafile))
		outpath = string("../data/diverse_rna_virus_structures/",fastafile)
		find_pdb_homologs(path, outpath)
	end
end


#=
pdbname = "1rvx"
chain = "A"
structure = retrievepdb(pdbname, pdb_dir=pdbdir)
polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])
=#