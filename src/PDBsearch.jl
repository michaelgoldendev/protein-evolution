using FastaIO
using Statistics
using BioAlignments
using BioStructures
using Random
using Printf

push!(LOAD_PATH,@__DIR__)
using Backbone
using Binaries
using DatasetCreator
using CommonUtils

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

function find_pdb_homologs(fastafile, outfile="homologs.fas", minseqs::Int=2)
	println("Find homologs: ",fastafile)
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
	alignmentscorecutoff = 0.2


	pdbmatches = Dict{AbstractString,AbstractString}()
	try
		file = open("../data/pdbsequences.fasta")
		name = ""
		seq = ""
		count = 0	
		k2counts = 0
		k3counts = 0
	    for (lineno,ln) in enumerate(eachline(file))
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
	    			namelower = lowercase(name)
	    			if !occursin("mutant", namelower) && !occursin("recombinant", namelower) && !occursin("synthetic", namelower) && !occursin("humanized", namelower) && !occursin("humanised", namelower)
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
			    	#key = string(">",pdbname,"\n",seq)
			    	key = string(pdbname)
			    	if !haskey(pdbmatches, key)
			    		resolution,rvalue,freervalue = get_quality_attributes(pdbname)
						if rvalue > 0.0 && rvalue < 0.25 && freervalue < 0.25
					    	println(lineno,"\t",count,"\t",pdbname,"\t",chain,"\t",maxmatch2, "\t", maxmatch3,"\t",maxalignmentscore)
					    	res = pairalign(GlobalAlignment(), seq, bestquery, AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1))
		    				aln = res.aln
		    				if maxalignmentscore > alignmentscorecutoff
		    					try
			    					structure = retrievepdb(pdbname, pdb_dir=pdbdir)
			    					polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])
			    					println(fout, string(">pdb",pdbname,"_",chain), "_r", @sprintf("%0.3f",rvalue) ,"_freer", @sprintf("%0.3f",freervalue))
			    					println(fout, polypeptide["sequence"])
					    			pdbmatches[key] = nameandchain
			    				catch
			    					println("ERROR DOWNLOADING ", pdbname)
			    				end
		    				end
					    end
				    end
			    end
	    	end
	       
	    end
	finally 
		println("FINISHED READING FILE")
	end
	close(fout)


	if length(pdbmatches) >= minseqs
		musclefile = string(outfile,".muscle.fas")
		muscle_alignment, cachefile = Binaries.muscle(outfile)
		fout = open(musclefile,"w")
		println(fout,muscle_alignment)
		close(fout)

		sequences = AbstractString[]
		names = AbstractString[]
		FastaIO.FastaReader(musclefile) do fr
			for (desc, seq) in fr
				push!(names,desc)
				push!(sequences,seq)
			end
		end
		println(sequences)

		#=
		newickstring, cachefile = Binaries.fasttreeaa(musclefile,branchsupport=true)
		fout = open(string(outfile,".nwk"),"w")
		println(fout,newickstring)
		close(fout)=#
	end

	println(length(pdbmatches))
	exit()
	return length(pdbmatches)
end

function find_random_pdbs(outfile="pdbselection.fas")
	rng = MersenneTwister(4917104813143)
	fout = open(outfile,"w")
	pdbmatches = Dict{AbstractString,AbstractString}()
	open("../data/pdbsequences.fasta") do file
		name = ""
		seq = ""
		count = 0
	    for ln in eachline(file)
	    	if startswith(ln,">")
	    		name = ln
	    		seq = ""
	    	else
	    		count += 1
	    		seq = ln

	    		if rand(rng,1:250) == 1 && match(r".*mol:protein.*", name) != nothing && length(seq) >= 100 
	    			namelower = lowercase(name)

	    			if !occursin("mutant", namelower) && !occursin("recombinant", namelower) && !occursin("synthetic", namelower) && !occursin("humanized", namelower) && !occursin("humanised", namelower)
				    	nameandchain = split(name[2:end])[1]
				    	spl2 = split(nameandchain,"_")
				    	pdbname = spl2[1]
				    	chain = spl2[2]
					    
					    if !haskey(pdbmatches, pdbname)
						    try
					    		resolution,rvalue,freervalue = DatasetCreator.get_quality_attributes(pdbname)
					    		if rvalue > 0.0 && rvalue < 0.25 && freervalue > 0.0 && freervalue < 0.25
						    		pdbmatches[pdbname] = pdbname
						    		println(fout, ">pdb$(pdbname)_$(chain) $(name[2:end])")
						    		println(fout, seq)
						    		flush(fout)
						    	end
						   	catch e

						   	end
					   end
				   end
		    	end
	    	end	       
	    end
	end
	close(fout)
end

function find_non_homologous_random_pdbs(outfile="nonhomologous.fas",homologfasta=nothing)
	kmers2 = Array{Int,1}[]
	kmers3 = Array{Int,1}[]
	sequences = AbstractString[]
	names = AbstractString[]
	if homologfasta != nothing
		FastaIO.FastaReader(homologfasta) do fr
		    for (desc, seq) in fr
		        push!(names,desc)
		        degapped = uppercase(replace(seq, "-" => ""))
		        degapped = replace(degapped, r"[^ACDEFGHIKLMNPQRSTVWY]" => "X")		        
		        degapped = replace(degapped, "X" => "")
		        push!(sequences, degapped)
		       	push!(kmers2, freqvectorint(degapped, 2))
		       	push!(kmers3, freqvectorint(degapped, 3))
		    end
		end
	end

	rng = MersenneTwister(4917104813143)
	fout = open(outfile,"w")
	pdbmatches = Dict{AbstractString,AbstractString}()
	open("../data/pdbsequences.fasta") do file
		name = ""
		seq = ""
		count = 0
	    for ln in eachline(file)
	    	if startswith(ln,">")
	    		name = ln
	    		seq = ""
	    	else
	    		count += 1
	    		seq = ln
	    		r = rand(rng,1:125)
	    		if (r == 1 || r == 16 || r == 65 || r == 32 || r == 48  || r == 72 || r == 98 || r == 114) && match(r".*mol:protein.*", name) != nothing && length(seq) >= 100 
	    			namelower = lowercase(name)

	    			if !occursin("mutant", namelower) && !occursin("recombinant", namelower) && !occursin("synthetic", namelower) && !occursin("humanized", namelower) && !occursin("humanised", namelower)
				    	nameandchain = split(name[2:end])[1]
				    	spl2 = split(nameandchain,"_")
				    	pdbname = spl2[1]
				    	chain = spl2[2]
					    
					    if !haskey(pdbmatches, pdbname)
					    	tv3 = freqvectorint(seq, 3)
					    	homologous = false
					    	k3cutoff = 0.1
					    	for qv3 in kmers3
						    	k3score = 1.0 - (sum(abs.(tv3 - qv3)) / sum(max.(tv3,qv3)))
					    		if k3score > k3cutoff
					    			homologous = true
					    			break
					    		end
					    	end
					    	if !homologous
							    try
						    		resolution,rvalue,freervalue = DatasetCreator.get_quality_attributes(pdbname)
						    		if rvalue > 0.0 && rvalue < 0.25 && freervalue > 0.0 && freervalue < 0.25
						    			push!(kmers3, tv3)
							    		pdbmatches[pdbname] = pdbname
							    		println(fout, ">pdb$(pdbname)_$(chain) $(name[2:end])")
							    		println(fout, seq)
							    		flush(fout)
							    		DatasetCreator.fromsinglepdb(pdbname, chain, "../data/nonhomologous_singles_xlarge/")
							    	end
							   	catch e

							   	end
						   end
					   end
				   end
		    	end
	    	end	       
	    end
	end
	close(fout)
end



function percentaligned(seq1::AbstractString, seq2::AbstractString)
	gapcounts = 0.0
	total = 0.0
	for (a,b) in zip(seq1,seq2)
		if a in aminoacids || b in aminoacids || a == 'X' || b == 'X'
			if a == '-' || b == '-'
				gapcounts += 1.0				
			end
			total += 1.0
		end
	end
	return (total-gapcounts)/total
end

function percentmatch(seq1::AbstractString, seq2::AbstractString)
	matchcount = 0.0
	total = 0.0
	for (a,b) in zip(seq1,seq2)
		if a in aminoacids || b in aminoacids
			if a == b
				matchcount += 1.0
			end
			total += 1.0
		end
	end
	return matchcounttotal
end

mutable struct ProteinAttributes
	rvalue::Float64
	freervalue::Float64

	function ProteinAttributes(rvalue::Float64=0.0, freervalue::Float64=0.0)
		new(rvalue, freervalue)
	end
end

function evaluate2(sequences::Array{AbstractString,1}, seqattr::Array{ProteinAttributes,1}, bitvector::Array{Int,1})
	minpercaligned = 1.0
	for i=1:length(sequences)
		if bitvector[i] > 0
			for j=i+1:length(sequences)
				if bitvector[j] > 0
					minpercaligned = min(minpercaligned, percentaligned(sequences[i], sequences[j]))
				end
			end
		end
	end

	minpairwisematch = 1.0
	maxpairwisematch = 0.0
	for i=1:length(sequences)
		if bitvector[i] > 0
			for j=i+1:length(sequences)
				if bitvector[j] > 0
	 				percmatch = percentmatch(sequences[i], sequences[j])
	 				minpairwisematch = min(minpairwisematch, percmatch)
	 				maxpairwisematch = max(maxpairwisematch, percmatch)
	 			end
	 		end
	 	end
	 end
	
	numpdbs = 0
	numseqs = 0
	if minpercaligned < 0.90 || minpairwisematch < 0.1 || maxpairwisematch > 0.90
		score = (0,0)
	else
		for i=1:length(sequences)
			if bitvector[i] > 0
				if seqattr[i].rvalue > 0.0 && seqattr[i].rvalue < 0.25 && seqattr[i].freervalue > 0.0 && seqattr[i].freervalue < 0.25
					numpdbs += 1
				else
					numseqs += 1
				end
			end
		end

		score = (numpdbs,numseqs)
	end
	return score[1]*length(sequences) + score[2]
end

function findmax2(sequences::Array{AbstractString,1}, seqattr::Array{ProteinAttributes,1}, evalualationfunc=evaluate2)
	rng = MersenneTwister(91093194923081)
	population = []
	popsize = 200

	v = ones(Int,length(sequences))
	push!(population, (v, evalualationfunc(sequences, seqattr, v)))
	for i=1:length(sequences)
		v = zeros(Int, length(sequences))
		v[i] = 1
		push!(population, (v, evalualationfunc(sequences, seqattr,v)))

		for j=i+1:length(sequences)
			v = zeros(Int, length(sequences))
			v[i] = 1
			v[j] = 1
			push!(population, (v, evaluate(sequences,seqattr,v)))
		end
	end

	while length(population) < popsize 
		v = zeros(Int, length(sequences))
		randorder = Int[i for i=1:length(sequences)]
		shuffle!(rng, randorder)
		for k=1:min(3, length(sequences))
			v[randorder[k]] = 1
		end
		push!(population, (v, evaluate(sequences,seqattr,v)))
	end

	for z=1:50
		scores = zeros(Float64, length(population))
		for i=1:length(population)
			scores[i] = population[i][2]
		end

		maxindex = 1
		for i=1:length(population)
			if population[i][2] > population[maxindex][2]
				maxindex = i
			end
		end
		
		newpopulation = []
		push!(newpopulation, population[maxindex])
		numsel = div(popsize,2)
		selection = zeros(Int, length(population))
		iter = 0
		while length(newpopulation) < numsel
			ind = CommonUtils.sample(rng,scores)
			if selection[ind] == 0 || iter > 1000
				push!(newpopulation, deepcopy(population[ind]))
				selection[ind] = 1
			end
			iter += 1
		end
		upper = div(popsize,3)*2
		while length(newpopulation) < upper
			ind = rand(rng, 1:length(newpopulation))
			newbitarray = copy(newpopulation[ind][1])
			randindex =  rand(rng, 1:length(newbitarray))
			newbitarray[randindex]= 1 - newbitarray[randindex]
			push!(newpopulation, (newbitarray, evaluate(sequences,seqattr,newbitarray)))
		end
		while length(newpopulation) < popsize
			ind1 = rand(rng, 1:length(newpopulation))
			ind2 = rand(rng, 1:length(newpopulation))
			newbitarray = copy(newpopulation[ind1][1])
			for i=1:length(newbitarray)
				if rand(rng) < 0.5
					newbitarray[i] = newpopulation[ind1][1][i]
				else
					newbitarray[i] = newpopulation[ind2][1][i]
				end
			end
			push!(newpopulation, (newbitarray, evaluate(sequences,seqattr,newbitarray)))
		end
		population = newpopulation
	end
	return population[1]
end

function evaluate(sequences::Array{AbstractString,1}, seqattr::Array{Any,1}, bitvector::Array{Int,1})
	minpercaligned = 1.0
	for i=1:length(sequences)
		if bitvector[i] > 0
			for j=i+1:length(sequences)
				if bitvector[j] > 0
					minpercaligned = min(minpercaligned, percentaligned(sequences[i], sequences[j]))
				end
			end
		end
	end
	
	numpdbs = 0
	numseqs = 0
	if minpercaligned < 0.90
		score = (0,0)
	else
		for i=1:length(sequences)
			if bitvector[i] > 0
				rvalue = seqattr[i][4]
				freervalue = seqattr[i][5]
				if rvalue > 0.0 && rvalue < 0.25 && freervalue > 0.0 && freervalue < 0.25
					numpdbs += 1
				else
					numseqs += 1
				end
			end
		end

		score = (numpdbs,numseqs)
	end
	return score[1]*length(sequences) + score[2]
end

function findmax(sequences::Array{AbstractString,1}, seqattr::Array{Any,1}, evalualationfunc=evaluate)
	rng = MersenneTwister(91093194923081)
	population = []
	popsize = 200

	v = ones(Int,length(sequences))
	push!(population, (v, evalualationfunc(sequences, seqattr, v)))
	for i=1:length(sequences)
		v = zeros(Int, length(sequences))
		v[i] = 1
		push!(population, (v, evalualationfunc(sequences, seqattr,v)))

		for j=i+1:length(sequences)
			v = zeros(Int, length(sequences))
			v[i] = 1
			v[j] = 1
			push!(population, (v, evaluate(sequences,seqattr,v)))
		end
	end

	while length(population) < popsize 
		v = zeros(Int, length(sequences))
		randorder = Int[i for i=1:length(sequences)]
		shuffle!(rng, randorder)
		for k=1:min(3, length(sequences))
			v[randorder[k]] = 1
		end
		push!(population, (v, evaluate(sequences,seqattr,v)))
	end

	for z=1:50
		scores = zeros(Float64, length(population))
		for i=1:length(population)
			scores[i] = population[i][2]
		end

		maxindex = 1
		for i=1:length(population)
			if population[i][2] > population[maxindex][2]
				maxindex = i
			end
		end
		
		newpopulation = []
		push!(newpopulation, population[maxindex])
		numsel = div(popsize,2)
		selection = zeros(Int, length(population))
		iter = 0
		while length(newpopulation) < numsel
			ind = CommonUtils.sample(rng,scores)
			if selection[ind] == 0 || iter > 1000
				push!(newpopulation, deepcopy(population[ind]))
				selection[ind] = 1
			end
			iter += 1
		end
		upper = div(popsize,3)*2
		while length(newpopulation) < upper
			ind = rand(rng, 1:length(newpopulation))
			newbitarray = copy(newpopulation[ind][1])
			randindex =  rand(rng, 1:length(newbitarray))
			newbitarray[randindex]= 1 - newbitarray[randindex]
			push!(newpopulation, (newbitarray, evaluate(sequences,seqattr,newbitarray)))
		end
		while length(newpopulation) < popsize
			ind1 = rand(rng, 1:length(newpopulation))
			ind2 = rand(rng, 1:length(newpopulation))
			newbitarray = copy(newpopulation[ind1][1])
			for i=1:length(newbitarray)
				if rand(rng) < 0.5
					newbitarray[i] = newpopulation[ind1][1][i]
				else
					newbitarray[i] = newpopulation[ind2][1][i]
				end
			end
			push!(newpopulation, (newbitarray, evaluate(sequences,seqattr,newbitarray)))
		end
		population = newpopulation
	end
	return population[1]
end

						

function create_from_homstrad(homstraddir="../data/homstrad_with_PDB_2019_Apr_1/", outputdir="../data/homstrad_curated_highquality/")
	if !isdir(outputdir)
		mkdir(outputdir)
	end
	countdatasets = 0
	folders = readdir(homstraddir)
	for folder in folders
		familyfolder = joinpath(homstraddir,folder)
		if isdir(familyfolder)
			familyfiles = filter(x -> endswith(x, ".ali"), readdir(familyfolder))
			if length(familyfiles) > 0
				namesandchains = []
				alifile	= joinpath(familyfolder, familyfiles[1])
				open(alifile) do file
					for ln in eachline(file)
			    		if startswith(ln,">")
			    			m = match(r".+;(....)(.)?.*", ln)
			    			pdbname = ""
			    			chain = ""
			    			if m[1] != nothing
			    				pdbname = m[1]
			    				if m[2] != nothing
			    					chain = uppercase(m[2])
			    				end
			    				push!(namesandchains, (pdbname,chain))
			    			end
			    		end
			    	end
			    end
			    if length(namesandchains) >= 2
			    	familyname = familyfiles[1][1:end-4]
			    	println("Family: ", familyname)
			    	selection = []
			    	usedpdbs = Dict{String,String}()
			    	for (pdbname, chain) in namesandchains
			    		println(pdbname)
			    		structure = retrievepdb(pdbname, pdb_dir=pdbdir)
			    		#=
			    		if !haskey(structure[1], chain)
			    			chain = ""
			    		end=#
			    		if chain == "" && length(structure) > 0
				    		longestsequence = ""
							for c in structure[1]
								polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(c)
								sequence = polypeptide["sequence"]
								if longestsequence == "" || length(c) > length(longestsequence)
									longestsequence =  sequence
									chain = chainid(c)
								end
							end
						end
						try 
				    		polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])
				    		resolution,rvalue,freervalue = DatasetCreator.get_quality_attributes(pdbname)
				    		if !haskey(usedpdbs,pdbname)
				    			push!(selection, (pdbname,chain,polypeptide["sequence"],rvalue,freervalue))
				    			usedpdbs[pdbname] = pdbname
				    		end
				    	catch 

				    	end
				    	#pdbname 
			    	  	#structure = retrievepdb(pdbname, pdb_dir=pdbdir)
		   				#polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])
				    	#println(namesandchains)
				    end
				    outpath = joinpath(outputdir, string(familyname,".fas"))
				    
				    if length(selection) >= 1
				    	selectiondict = Dict{AbstractString,Any}()
				    	fout = open(outpath, "w")
				    	for s in selection
				    		seqname = ""
				    		resolution,rvalue,freervalue = DatasetCreator.get_quality_attributes(s[1])
			    			if rvalue > 0.0 && rvalue < 0.25 && freervalue > 0.0 && freervalue < 0.25
			    			#if rvalue > 0.0 && rvalue < 0.25
			    				seqname = string(">pdb",s[1],"_",s[2])
						    	println(fout, seqname)
						    	println(fout,s[3])
						    else
						    	seqname = string(">seq",s[1],"_",s[2])
						    	println(fout,seqname)
						    	println(fout,s[3])
						    end
						    selectiondict[seqname[2:end]] = s
					    end
				    	close(fout)

				    	musclefile = joinpath(outputdir, string(familyname,".muscle.fas"))
						muscle_alignment, cachefile = Binaries.muscle(outpath)
						fout = open(musclefile,"w")
						println(fout,muscle_alignment)
						close(fout)

						orderedselection = []
						sequences = AbstractString[]
						names = AbstractString[]
						FastaIO.FastaReader(musclefile) do fr
						    for (desc, seq) in fr
						        push!(names, desc)
						        push!(orderedselection, selectiondict[desc])
						        push!(sequences, seq)
						    end
						end
						
						bitarray = findmax(sequences, orderedselection)[1]
						selectionfile = joinpath(outputdir, string(familyname,".selection.fas"))
						fout = open(selectionfile,"w")
						for i=1:length(bitarray)
							if bitarray[i] > 0
								println(fout,">$(names[i])")
								println(fout,sequences[i])
							end
						end
						close(fout)

						musclefile = joinpath(outputdir, string(familyname,".selection.muscle.fasta"))
						muscle_alignment, cachefile = Binaries.muscle(selectionfile)
						fout = open(musclefile,"w")
						println(fout,muscle_alignment)
						close(fout)

						familyoutfile = joinpath(outputdir, string(familyname,".selection.fam"))
						DatasetCreator.fromsequencealignment(musclefile, familyoutfile)
				    end
			    end
			end
		end
		#println(familyfolder)
	end
end

function find_random_pdb_families()
	rng = MersenneTwister(4917104813143)
	pdbmatches = Dict{AbstractString,AbstractString}()
	open("../data/pdbsequences.fasta") do file
		name = ""
		seq = ""
		count = 0
	    for ln in eachline(file)
	    	if startswith(ln,">")
	    		name = ln
	    		seq = ""
	    	else
	    		count += 1
	    		seq = ln

	    		if rand(rng,1:500) == 125 && match(r".*mol:protein.*", name) != nothing && length(seq) >= 100 
	    			namelower = lowercase(name)

	    			if !occursin("mutant", namelower) && !occursin("recombinant", namelower) && !occursin("synthetic", namelower) && !occursin("humanized", namelower) && !occursin("humanised", namelower)
				    	nameandchain = split(name[2:end])[1]
				    	spl2 = split(nameandchain,"_")
				    	pdbname = spl2[1]
				    	chain = spl2[2]
					    
					    if !haskey(pdbmatches, pdbname)
						    try
					    		resolution,rvalue,freervalue = DatasetCreator.get_quality_attributes(pdbname)
					    		if rvalue > 0.0 && rvalue < 0.25 && freervalue > 0.0 && freervalue < 0.25
						    		pdbmatches[pdbname] = pdbname

						    		#=
						    		println(fout, ">pdb$(pdbname)_$(chain) $(name[2:end])")
						    		println(fout, seq)
						    		flush(fout)=#
						    		tempfile = string("../data/scratch/",pdbname,".temp.fas")
						    		if !Base.Filesystem.isfile(tempfile)
							    		outfile = string("../data/scratch/",pdbname,".family.fas")
							    		fastafile = DatasetCreator.getpdbsequencealignment(String[pdbname], tempfile, false)
										numseqs = find_pdb_homologs(fastafile, outfile)
										println("FINISHED ", pdbname,"\t",numseqs)
										
										if numseqs <= 1
											Base.Filesystem.rm(fastafile)
											Base.Filesystem.rm(tempfile)
											Base.Filesystem.rm(outfile)
										end
									end
						    	end
						   	catch e

						   	end
					   end
				   end
		    	end
	    	end	       
	    end
	end
end

#find_non_homologous_random_pdbs("nonhomologous.fasta", abspath("selectedsequences.fasta"))
#create_from_homstrad()
find_random_pdb_families()
#find_random_pdbs()

#fastafile = abspath("../data/influenza_a/HA/selection3.fasta")
#fastafile = abspath("../data/hiv/curated6.fasta")
#fastafile = abspath("../data/test_data/hiv_pol_selection.fasta")
#fastafile = abspath("../data/test_data/maise_streak_virus_coat_protein_selection.fasta")
#fastafile = abspath("../data/test_data/westnile_dengue_selection.fasta")
#
#fastafile = "../data/diverse_rna_virus_structures/norovirus_capsid.select.fasta"

#fastafile = DatasetCreator.getpdbsequencealignment(String["3u4f"], "../data/scratch/temp.fas")
#find_pdb_homologs(fastafile, "../data/scratch/output.fas")


#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\curated_rna\\avian_coronavirus_NSP3.select.fasta"
#family_directories = ["../data/curated/curated_rna/", "../data/selected_families/"]

#=
family_directories = ["../data/curated_selection"]
for family_dir in family_directories
	fastafiles = filter(f -> endswith(f,".fasta"), readdir(family_dir))
	for fastafile in fastafiles
		path = abspath(joinpath(family_dir,fastafile))
		outpath = string("../data/diverse_rna_virus_structures/",fastafile)
		find_pdb_homologs(path, outpath)
	end
end
=#


#=
pdbname = "1rvx"
chain = "A"
structure = retrievepdb(pdbname, pdb_dir=pdbdir)
polypeptide = Backbone.backbone_angles_and_bond_lengths_from_pdb(structure[1][chain])
=#