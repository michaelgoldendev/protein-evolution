include("Main.jl")

function estimatehmm(rng::AbstractRNG, trainingexamples, modelparams::ModelParams, numiters::Int, parsed_args)
	loglikelihoods = Float64[]
	for iter=1:numiters
		if !(parsed_args["angles-cond-aa"] <= 0) && iter >= parsed_args["angles-cond-aa"]
			modelparams.hidden_conditional_on_aa = true
		end
		println("PRECLUSTER ",iter)

		loglikelihood = 0.0
		for (index, (proteins,nodelist,json_family,sequences)) in enumerate(trainingexamples) 
			for node in nodelist
				if isleafnode(node)
					#backwardsamplefast(rng, node, modelparams)				
					loglikelihood += backwardsamplesingle(rng, node, modelparams)
				end
			end
		end
		push!(loglikelihoods, loglikelihood)
		println("Log-likelihoods: ", loglikelihoods)

		transitionrate_counts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.01
		transitionrate_totals = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.001
		for h=1:modelparams.numhiddenstates
			transitionrate_counts[h,h] = 0.0
		end
		transitionrate_events = ones(Float64, modelparams.numhiddenstates)*0.01
		transitionrate_times = ones(Float64, modelparams.numhiddenstates)*0.001

		aatransitionrate_events = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		aatransitionrate_times = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.001

		aatransitionrate_counts = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		aatransitionrate_totals = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		for aa=1:modelparams.alphabet
			aatransitionrate_counts[aa,aa] = 0.0
		end		
		accepted_hidden = 0.0
		accepted_hidden_total = 0.0
		accepted_aa = 0.0
		accepted_aa_total = 0.0
		for (trainingindex, (proteins,nodelist,json_family,sequences)) in enumerate(trainingexamples)
			numcols = length(proteins[1])

			for node in nodelist
				if isleafnode(node)
					for col=1:numcols-1					
						modelparams.transitioncounts[node.data.branchpath.paths[col][end], node.data.branchpath.paths[col+1][end]] += 1.0					
					end					
					for col=1:numcols
						modelparams.hiddennodes[node.data.branchpath.paths[col][end]].aa_node.counts[node.data.aabranchpath.paths[col][end]] += 1.0	
					end
				end
			end

			for h=1:modelparams.numhiddenstates
				modelparams.hiddennodes[h].bond_lengths_node.data = Array{Float64,1}[]
			end

			for col=1:numcols
				for node in nodelist				
					if isleafnode(node)
						h = node.data.branchpath.paths[col][end]
						site = proteins[node.seqindex].sites[col]
						if site.aa > 0
							#modelparams.hiddennodes[h].aa_node.counts[site.aa] += 1.0
							push!(modelparams.hiddennodes[h].phi_node.data, site.phi)
							push!(modelparams.hiddennodes[h].omega_node.data, site.omega)
							push!(modelparams.hiddennodes[h].psi_node.data, site.psi)
							BivariateVonMises.add_bvm_point(modelparams.hiddennodes[h].phipsi_node, Float64[site.phi, site.psi])
							push!(modelparams.hiddennodes[h].phi_nodes[site.aa].data, site.phi)
							push!(modelparams.hiddennodes[h].omega_nodes[site.aa].data, site.omega)
							push!(modelparams.hiddennodes[h].psi_nodes[site.aa].data, site.psi)
							BivariateVonMises.add_bvm_point(modelparams.hiddennodes[h].phipsi_nodes[site.aa], Float64[site.phi, site.psi])
							push!(modelparams.hiddennodes[h].bond_angle1_node.data, site.bond_angle1)
							push!(modelparams.hiddennodes[h].bond_angle2_node.data, site.bond_angle2)
							push!(modelparams.hiddennodes[h].bond_angle3_node.data, site.bond_angle3)
							add_point(modelparams.hiddennodes[h].bond_lengths_node, Float64[site.bond_length1, site.bond_length2, site.bond_length3])
						end
					end
				end
			end			
		end
		estimate_hidden_transition_probs(modelparams)

		for h=1:modelparams.numhiddenstates			
			estimate_categorical(modelparams.hiddennodes[h].aa_node, 1.0)		
			estimatevonmises(modelparams.hiddennodes[h].phi_node)
			estimatevonmises(modelparams.hiddennodes[h].psi_node)
			estimatevonmises(modelparams.hiddennodes[h].omega_node)
			estimate_bvm(modelparams.hiddennodes[h].phipsi_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle1_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle2_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle3_node)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].phipsi_node.mu,"\t",modelparams.hiddennodes[h].phipsi_node.k, "\t", modelparams.hiddennodes[h].phi_node.mu, "\t", modelparams.hiddennodes[h].psi_node.mu, "\t", modelparams.hiddennodes[h].phi_node.kappa, "\t", modelparams.hiddennodes[h].psi_node.kappa)		
			println(modelparams.hiddennodes[h].phipsi_node.count,"\t",modelparams.hiddennodes[h].phi_node.N,"\t",modelparams.hiddennodes[h].psi_node.N)
			#println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle1_node.mu,"\t",modelparams.hiddennodes[h].bond_angle1_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle1_node.N)
			#println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle2_node.mu,"\t",modelparams.hiddennodes[h].bond_angle2_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle2_node.N)
			#println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle3_node.mu,"\t",modelparams.hiddennodes[h].bond_angle3_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle3_node.N)
			#println(iter,"\t",h,"\t", modelparams.hiddennodes[h].bond_lengths_node.mvn.μ)
			estimate_multivariate_node(modelparams.hiddennodes[h].bond_lengths_node)
			for aa=1:modelparams.alphabet
				estimatevonmises(modelparams.hiddennodes[h].phi_nodes[aa])
				estimatevonmises(modelparams.hiddennodes[h].psi_nodes[aa])
				estimatevonmises(modelparams.hiddennodes[h].omega_nodes[aa])
				estimate_bvm(modelparams.hiddennodes[h].phipsi_nodes[aa])
				if modelparams.hidden_conditional_on_aa
					print(iter,"\t",h,"\t",aminoacids[aa],"\t",@sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]),"\t",modelparams.hiddennodes[h].phi_nodes[aa].mu,"\t",modelparams.hiddennodes[h].phi_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].phi_nodes[aa].N)		
					println("\t",modelparams.hiddennodes[h].psi_nodes[aa].mu,"\t",modelparams.hiddennodes[h].psi_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].psi_nodes[aa].N)
					#println(iter,"\t",h,"\t",aminoacids[aa],"\t",@sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]),"\t",modelparams.hiddennodes[h].omega_nodes[aa].mu,"\t",modelparams.hiddennodes[h].omega_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].omega_nodes[aa].N)						end
				end
			end
			for aa=1:20
				println(iter,"\t",h,"\t",aminoacids[aa],"\t", @sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]))
			end			
		end
	end

	reset_matrix_cache(modelparams)
end

function datasummary(countdict, keylist)
	totalaminoacidcounts = 0
	totalstructuralcounts = 0
	totalsequences = 0
	totalstructures = 0
	for key in keylist
		combinationcounts, aminoacidcounts, structuralcounts = countdict[key]
		totalaminoacidcounts += aminoacidcounts
		totalstructuralcounts += structuralcounts
		totalsequences += key[1]*combinationcounts
		totalstructures += key[2]*combinationcounts
	end
	println("Number in category: ", sum(Int[countdict[key][1] for key in keylist]))
	println("Amino acid site observations: $(totalaminoacidcounts) (avg. $(@sprintf("%0.1f", totalaminoacidcounts/totalsequences)))")
	println("Structural site observations: $(totalstructuralcounts) (avg. $(@sprintf("%0.1f", totalstructuralcounts/totalstructures)))")
	println("Total sequences: ", totalsequences)
	println("Total structures: ", totalstructures)
end

function printsummary(countdict)
	println("Exactly one sequence")
	datasummary(countdict, Tuple{Int,Int}[key for key in filter(x -> x[1] == 1 && x[2] == 0, keys(countdict))])
	println("")

	println("Exactly two sequences")
	datasummary(countdict, Tuple{Int,Int}[key for key in filter(x -> x[1] == 2 && x[2] == 0, keys(countdict))])
	println("")

	println("Three or more sequences")
	datasummary(countdict, Tuple{Int,Int}[key for key in filter(x -> x[1] >= 3 && x[2] == 0, keys(countdict))])
	println("")

	println("Exactly one sequence and structure")
	datasummary(countdict, Tuple{Int,Int}[key for key in filter(x -> x[1] == 1 && x[2] == 1, keys(countdict))])
	println("")

	println("Exactly one structure (and two or more sequences)")
	datasummary(countdict, Tuple{Int,Int}[key for key in filter(x -> x[1] >= 2 && x[2] == 1, keys(countdict))])
	println("")

	println("Exactly two structures")
	datasummary(countdict, Tuple{Int,Int}[key for key in filter(x -> x[2] == 2, keys(countdict))])
	println("")

	println("Three or more structures")
	datasummary(countdict, Tuple{Int,Int}[key for key in filter(x -> x[2] >= 3, keys(countdict))])
	println("")

	println("Total")
	datasummary(countdict, Tuple{Int,Int}[key for key in keys(countdict)])
	println("")
	println("================================================")
end

function stripfamily(json_family)
	startindex = typemax(Int)
	endindex = 1
	for protein in json_family["proteins"]
		alignedsequence = protein["aligned_sequence"]
		#println(length(alignedsequence))
		#println(alignedsequence)
		for i=1:length(alignedsequence)
			if !(alignedsequence[i] == '-' || alignedsequence[i] == 'X')
				startindex = min(startindex, i)
				break
			end
		end
		for i=length(alignedsequence):-1:1
			if !(alignedsequence[i] == '-' || alignedsequence[i] == 'X')
				endindex = max(endindex, i)
				break
			end
		end
	end	
	println("startindex: ", startindex)
	println("endindex: ", endindex)	
	startgaps = Int[]
	endgaps = Int[]
	for protein in json_family["proteins"]
		alignedsequence = protein["aligned_sequence"]
		push!(startgaps, count(x -> x == 'X', alignedsequence[1:startindex-1]))
		push!(endgaps, count(x -> x == 'X', alignedsequence[endindex:end]))
	end
	for (index, protein) in enumerate(json_family["proteins"])
		protein["sequence"] = protein["sequence"][startgaps[index]+1:end-endgaps[index]]
		protein["aligned_sequence"] = protein["aligned_sequence"][startindex:endindex]

		protein["bond_angles"] = protein["bond_angles"][startgaps[index]+1:end-endgaps[index]]
		protein["aligned_bond_angles"] = protein["aligned_bond_angles"][startindex:endindex]

		protein["bond_lengths"] = protein["bond_lengths"][startgaps[index]+1:end-endgaps[index]]
		protein["aligned_bond_lengths"] = protein["aligned_bond_lengths"][startindex:endindex]
		
		protein["omega"] = protein["omega"][startgaps[index]+1:end-endgaps[index]]
		protein["aligned_omega"] = protein["aligned_omega"][startindex:endindex]

		protein["phi_psi"] = protein["phi_psi"][startgaps[index]+1:end-endgaps[index]]
		protein["aligned_phi_psi"] = protein["aligned_phi_psi"][startindex:endindex]

		if haskey(protein, "Ntempfactor")
			protein["Ntempfactor"] = protein["Ntempfactor"][startgaps[index]+1:end-endgaps[index]]
			protein["aligned_Ntempfactor"] = protein["aligned_Ntempfactor"][startindex:endindex]

			protein["CAtempfactor"] = protein["CAtempfactor"][startgaps[index]+1:end-endgaps[index]]
			protein["aligned_CAtempfactor"] = protein["aligned_CAtempfactor"][startindex:endindex]

			protein["Ctempfactor"] = protein["Ctempfactor"][startgaps[index]+1:end-endgaps[index]]
			protein["aligned_Ctempfactor"] = protein["aligned_Ctempfactor"][startindex:endindex]
		end
	end 
	#println(startgaps,"\t",endgaps)

end

function loadtrainingexamples(rng::AbstractRNG, parsed_args, family_directories, modelparams::ModelParams)
	countdict = Dict{Tuple{Int,Int},Tuple{Int,Int,Int}}()
	family_names = String[]	
	trainingexamples = Tuple[]
	traininghash = zero(UInt)
	for family_dir in family_directories
		family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))
		for family_file in family_files
			full_path = abspath(joinpath(family_dir, family_file))
			json_family = JSON.parse(open(full_path, "r"))
			stripfamily(json_family)
			traininghash = hash(json_family, traininghash)
			if 1 <= length(json_family["proteins"]) <= 1e10
				training_example = training_example_from_json_family(rng, modelparams, json_family)
				if rand(rng) < 0.01
					println(json_family["newick_tree"])
				end
				push!(trainingexamples, training_example)
				#println(getnewick(training_example[2][1]))
				push!(family_names, family_file)

				numsequences = 0
				numstructures = 0
				aminoacidcounts = 0
				structuralcounts = 0
				for p in json_family["proteins"]
					hasstructure = false
					for (phi,psi) in p["phi_psi"]
						if phi > -100.0 || psi > -100.0
							hasstructure = true
							structuralcounts += 1
						end
					end
					for aa in p["sequence"]
						if aa != 'X' && aa != '-'
							aminoacidcounts += 1
						end
					end

					numsequences += 1
					if hasstructure
						numstructures += 1
					end
				end
				key = (numsequences, numstructures)
				vals = get(countdict, key, (0,0,0))
				vals = (vals[1]+1, vals[2]+aminoacidcounts, vals[3]+structuralcounts)
				countdict[key] = vals

				if rand(rng) < 0.01
					#println(countdict)
					#println("Exactly one sequence: ", sum(Int[countdict[key][1] for key in filter(x -> x[1] == 1 && x[2] == 0, keys(countdict))]))
					#println("Exactly two sequences: ", sum(Int[countdict[key][1] for key in filter(x -> x[1] == 2 && x[2] == 0, keys(countdict))]))
					#println("Three or more sequences: ", sum(Int[countdict[key][1] for key in filter(x -> x[1] >= 3 && x[2] == 0, keys(countdict))]))
					#println("Exactly one sequence and structure: ", sum(Int[countdict[key][1] for key in filter(x -> x[1] == 1 && x[2] == 1, keys(countdict))]))
					#println("Exactly one structure (and two or more sequences): ", sum(Int[countdict[key][1] for key in filter(x -> x[1] >= 2 && x[2] == 1, keys(countdict))]))
					#println("Exactly two structures: ", sum(Int[countdict[key][1] for key in filter(x -> x[2] == 2, keys(countdict))]))
					#println("Three or more structures: ", sum(Int[countdict[key][1] for key in filter(x -> x[2] >= 3, keys(countdict))]))
					printsummary(countdict)
				end

				if parsed_args["traininginstances"] != nothing && length(trainingexamples) == parsed_args["traininginstances"]
					break
				end
			end
		end
		if parsed_args["traininginstances"] != nothing && length(trainingexamples) == parsed_args["traininginstances"]
			break
		end
	end
	printsummary(countdict)

	traininghashbase36 = string(traininghash, base=36)
	return family_names,trainingexamples,traininghashbase36
end

function sampletraininginstances(iter::Int, rng::AbstractRNG, trainingexamples, modelparams::ModelParams; maxsamplesperiter::Int=500, familyiter::Int=5, sitethreshold::Int=2, dosamplesiterates::Bool=true, samplebranchlengths::Bool=true, family_names::Array{String,1})
	totalbranchlength_output = 0.0
	accepted_hidden = 0.0
	accepted_hidden_total = 0.0
	accepted_aa = 0.0
	accepted_aa_total = 0.0
	totalhiddentime = 0.0
	totalaatime = 0.0		
	starttime = time()
	for (trainingindex, (proteins,nodelist,json_family,sequences)) in enumerate(trainingexamples)
		numcols = length(proteins[1])

		if length(proteins) == 1
			backwardsamplesingle(rng, nodelist[1], modelparams)
		else
			maxsamplesthisiter = maxsamplesperiter
			if (trainingindex+iter) % familyiter != 0
				maxsamplesthisiter = 1
			end
			samplehiddenstates = true
							#=
			maxsamplesthisiter = maxsamplesperiter
			samplehiddenstates = (trainingindex+iter) % familyiter == 0
			if !samplehiddenstates
				maxsamplesthisiter = 1
			end=#
			accepted = zeros(Int, numcols)			
			for i=1:maxsamplesthisiter					
				randcols = shuffle(rng, Int[i for i=1:numcols])
				for col in randcols
					if accepted[col] < sitethreshold || i % 20 == 0 || (col > 1 && accepted[col-1] < sitethreshold) || (col < numcols && accepted[col+1] < sitethreshold) 
						a1,a2,a3,a4, hidden_accepted, aa_accepted, hiddentime, aatime = samplepaths_seperate_new(rng,col,proteins,nodelist, modelparams, samplehiddenstates=samplehiddenstates, dosamplesiterates=dosamplesiterates)
						totalhiddentime += hiddentime
						totalaatime += aatime
						accepted_hidden += a1
						accepted_hidden_total += a2
						accepted_aa += a3
						accepted_aa_total += a4
						if hidden_accepted
							accepted[col] += 1
						end						
					end
				end
				min_accepted = minimum(accepted)
				if min_accepted >= sitethreshold
					elapsedtime = time()-starttime					
					println("min_accepted ", min_accepted," out of ", i, " mean is ", mean(accepted))
					println(totalhiddentime,"\t",totalaatime,"\t",family_names[trainingindex],"\t", elapsedtime)
					break
				end


				if samplebranchlengths && (i <= sitethreshold || i % 10 == 0)
					for node in nodelist
						if !isroot(node)
							t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
							node.branchlength = t
							#events = mean([length(p)-1.0 for p in node.data.branchpath.paths])
							#aa_events = mean([length(p)-1.0 for p in node.data.aabranchpath.paths])
							#println(node.nodeindex,"\t",node.data.inputbranchlength,"\t",node.branchlength,"\t",events,"\t",aa_events)
						end
					end
				end
			end
		end

		for col=1:numcols-1
			for node in nodelist
				if isroot(node)
					h1 = node.data.branchpath.paths[col][end]
					h2 = node.data.branchpath.paths[col+1][end]
					modelparams.transitioncounts[h1,h2] += 1.0
				else
					branchiterator = BranchPathIterator(node.data.branchpath, Int[col,col+1])
					for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
						#dt = (currtime-prevtime)*node.branchlength
						if changecol == 2
							modelparams.transitioncounts[currstates[1],currstates[2]] += 1.0
						end
						#modelparams.transitioncounts[prevstates[1],prevstates[2]] += dt
					end
				end
			end
		end

		for col=1:numcols
			for node in nodelist
				if isroot(node)						
					modelparams.hiddennodes[node.data.branchpath.paths[col][end]].aa_node.counts[node.data.aabranchpath.paths[col][end]] += 1.0
				else
					multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,Int[col]), BranchPathIterator(node.data.aabranchpath,Int[col])])
					hiddeniter = multi_iter.branchpathiterators[1]
					aaiter = multi_iter.branchpathiterators[2]
					for it in multi_iter
						if multi_iter.branchpathindex == 2 && aaiter.mincol == 1
							modelparams.hiddennodes[hiddeniter.prevstates[1]].aa_node.counts[aaiter.currstates[1]] += 1.0								
						end
					end
				end
			end
		end

		for h=1:modelparams.numhiddenstates
			modelparams.hiddennodes[h].bond_lengths_node.data = Array{Float64,1}[]
		end

		for col=1:numcols
			for node in nodelist				
				if isleafnode(node)
					h = node.data.branchpath.paths[col][end]
					site = proteins[node.seqindex].sites[col]
					if site.aa > 0
						#modelparams.hiddennodes[h].aa_node.counts[site.aa] += 1.0
						push!(modelparams.hiddennodes[h].phi_node.data, site.phi)
						push!(modelparams.hiddennodes[h].omega_node.data, site.omega)
						push!(modelparams.hiddennodes[h].psi_node.data, site.psi)
						BivariateVonMises.add_bvm_point(modelparams.hiddennodes[h].phipsi_node, Float64[site.phi, site.psi])
						push!(modelparams.hiddennodes[h].phi_nodes[site.aa].data, site.phi)
						push!(modelparams.hiddennodes[h].omega_nodes[site.aa].data, site.omega)
						push!(modelparams.hiddennodes[h].psi_nodes[site.aa].data, site.psi)
						BivariateVonMises.add_bvm_point(modelparams.hiddennodes[h].phipsi_nodes[site.aa], Float64[site.phi, site.psi])
						push!(modelparams.hiddennodes[h].bond_angle1_node.data, site.bond_angle1)
						push!(modelparams.hiddennodes[h].bond_angle2_node.data, site.bond_angle2)
						push!(modelparams.hiddennodes[h].bond_angle3_node.data, site.bond_angle3)
						add_point(modelparams.hiddennodes[h].bond_lengths_node, Float64[site.bond_length1, site.bond_length2, site.bond_length3])
					end
				end
			end
		end			
	end
	return totalbranchlength_output,accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total,totalhiddentime,totalaatime
end

function train(parsed_args=Dict{String,Any}()) 
	rng = MersenneTwister(85649871544631)
	Random.seed!(rand(rng,UInt))

	numhiddenstates = parsed_args["numhiddenstates"]
	maxsamplesperiter = parsed_args["maxsamplesperiter"]
	sitethreshold = parsed_args["sitethreshold"]


	learnrates = true
	samplebranchlengths = parsed_args["samplebranchlengths"]
	independentsites = parsed_args["independentsites"]
	usestructureobservations = !parsed_args["sequenceonly"]
	#family_dir = "../data/families/"
	dosamplesiterates = parsed_args["samplesiterates"]
	uselgrates = parsed_args["uselgrates"]
	hiddenaascaling = parsed_args["hiddenaascaling"]

	modelparams = ModelParams(LGmatrix,numhiddenstates,1.0)	
	maxloglikelihood = -Inf
	noimprovement = 0


	if dosamplesiterates
		modelparams.rate_alpha = 0.5
		modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)		
	end	

	modelparams.usestructureobservations = usestructureobservations
	modelparams.hidden_conditional_on_aa = false	
	modelparams.use_bivariate_von_mises = true	 
	modelparams.ratemode = parsed_args["ratemode"]
	

	#family_directories = ["../data/curated/curated_rna/", "../data/selected_families/"]
	#family_directories = ["../data/homstrad_families/", "../data/curated_rna_virus_structures/"]	
	#family_directories = ["../data/single_pdbs/", "../data/homstrad_families/", "../data/curated_rna_virus_structures/"]
	#family_directories = ["../data/single_pdbs/", "../data/homstrad_families/"]
	#family_directories = ["../data/single_pdbs/", "../data/homstrad_families/", "../data/curated/curated_rna/"]
	#family_directories = ["../data/single_pdbs/", "../data/homstrad_families/", "../data/curated_rna_virus_structures/"] 
	#family_directories = ["../data/single_pdbs/", "../data/homstrad_families/"] 
	#family_directories = ["../data/homstrad_curated/", "../data/curated_rna_virus_structures/"]	
	family_directories = ["../data/homstrad_curated_highquality/", "../data/curated_rna_virus_structures/", "../data/nonhomologous_singles_xlarge/"]
	#family_directories = ["../data/homstrad_curated_highquality/", "../data/curated_rna_virus_structures/"]	
	familyiter = 1

	family_names,trainingexamples,traininghashbase36 = loadtrainingexamples(rng, parsed_args, family_directories, modelparams)

	#=
	countdict = Dict{Tuple{Int,Int},Int}()
	for (index, (proteins,nodelist,json_family,sequences)) in enumerate(trainingexamples)
		numsequences = 0
		numstructures = 0
		for p in json_family["proteins"]
			hasstructure = false
			for (phi,psi) in p["phi_psi"]
				if phi > -100.0 && psi > -100.0
					hasstructure = true
					break
				end
			end

			numsequences += 1
			if hasstructure
				numstructures += 1
			end
		end
		key = (numsequences, numstructures)
		numinstances = get(countdict, key, 0)
		numinstances += 1
		countdict[key] = numinstances
		println(countdict)
		println("Sequences only: ", sum(Int[countdict[key] for key in filter(x -> x[2] == 0, keys(countdict))]))
		println("Exactly one structure: ", sum(Int[countdict[key] for key in filter(x -> x[2] == 1, keys(countdict))]))
		println("Exactly two structures: ", sum(Int[countdict[key] for key in filter(x -> x[2] == 2, keys(countdict))]))
		println("Three or more structures: ", sum(Int[countdict[key] for key in filter(x -> x[2] >= 3, keys(countdict))]))
	end
	exit()
	=#

	#=
	selectionout = open("selectedsequences.fasta","w")
	count = 1
	for (proteins,nodelist,json_family,sequences) in trainingexamples
		for protein in json_family["proteins"]
			println(selectionout,">",protein["name"])
			println(selectionout,protein["sequence"])
			count += 1
		end
	end
	close(selectionout)=#

	outputmodelname = string(".h",modelparams.numhiddenstates, ".",traininghashbase36)
	outputmodelname = string(outputmodelname,".thresh",sitethreshold)
	if dosamplesiterates
		outputmodelname = string(outputmodelname,".samplesiterates")
	end
	if uselgrates
		outputmodelname = string(outputmodelname,".lgrates")
	end
	if hiddenaascaling
		outputmodelname = string(outputmodelname,".hiddenaascaling")
	else
		outputmodelname = string(outputmodelname,".nohiddenaascaling")
	end
	if !samplebranchlengths
		outputmodelname = string(outputmodelname,".nobranchsampling")
	end
	if independentsites
		outputmodelname = string(outputmodelname,".independentsites")
	end
	if parsed_args["sequenceonly"]
		outputmodelname = string(outputmodelname,".seqonly")
	end
	if parsed_args["angles-cond-aa"] > 0
		outputmodelname = string(outputmodelname,".anglescondaa", parsed_args["angles-cond-aa"])
	end
	if parsed_args["precluster"]
		outputmodelname = string(outputmodelname,".precluster")
	end
	outputmodelname = string(outputmodelname,".ratemode", modelparams.ratemode)
	modelfile = string("models/model", outputmodelname, ".model")
	tracefile = string("models/trace", outputmodelname,".log")
	familytracefile = string("models/family.trace", outputmodelname,".log")

	if parsed_args["loadfromcache"]
		modelparams = Serialization.deserialize(open(modelfile, "r"))
	end

	if modelparams.numhiddenstates == 1
		fout = open(modelfile, "w")
		reset_matrix_cache(modelparams)
		Serialization.serialize(fout, modelparams)
		close(fout)
		exit()
	end

	if parsed_args["precluster"]
		estimatehmm(rng, trainingexamples, modelparams, 30, parsed_args)
		family_names,trainingexamples,traininghashbase36 = loadtrainingexamples(rng, parsed_args, family_directories, modelparams)
	end
	
	#=
	usedindices = Int[]
	selected_trainingexamples = Tuple[]
	selected_family_names = String[]
	for (index, (family_name,training_example)) in enumerate(zip(family_names,trainingexamples))
		if length(training_example[2]) == 1
			push!(selected_trainingexamples, training_example)
			push!(selected_family_names, family_name)
			push!(usedindices, index)
		end
	end
	
	permutedindices = randperm(rng,length(family_names))
	permutedindices	= permutedindices[1:200]
	for index in permutedindices
		if !(index in usedindices)
			push!(selected_trainingexamples, trainingexamples[index])
			push!(selected_family_names, family_names[index])
		end
	end
	trainingexamples = selected_trainingexamples
	family_names = selected_family_names
	=#


	familytracefile = string("models/family.trace", outputmodelname,".log")
	familyoutfile = open(familytracefile,"w")
	print(familyoutfile,"iter")
	print(familyoutfile,"\tll\taall\tstructurell\tpathll")
	for name in family_names
		print(familyoutfile,"\t",name,"_ll","\t",name,"_obsll","\t",name,"_aall")
	end
	println(familyoutfile)

	println("Number of training examples: ", length(trainingexamples))

	totalbranchlength_input = 0.0
	for training_example in trainingexamples
		for node in training_example[2]
			if !isroot(node)
				totalbranchlength_input += node.branchlength
			end
		end
	end

	for (index, (proteins,nodelist,json_family,sequences)) in enumerate(trainingexamples)
		if index % 50 == 0
			println(index)
		end 
		if length(proteins) == 1
			backwardsamplesingle(rng, nodelist[1], modelparams)
		else
			numcols = length(proteins[1])
			for i=1:2
				randcols = shuffle(rng, Int[i for i=1:numcols])
				for col in randcols
					samplepaths_seperate_new(rng,col,proteins,nodelist, modelparams, dosamplesiterates=dosamplesiterates, accept_everything=true)
				end
			end

			randcols = shuffle(rng, Int[i for i=1:numcols])
			for col in randcols
				samplepaths_seperate_new(rng,col,proteins,nodelist, modelparams, dosamplesiterates=dosamplesiterates)
			end
		end
	end 

	logwriter = open(tracefile, "w")
	println(logwriter, "iter\tll\taall\tstructurell\tpathll")
	for iter=1:10000
		if !(parsed_args["angles-cond-aa"] <= 0) && iter >= parsed_args["angles-cond-aa"]
			modelparams.hidden_conditional_on_aa = true
		end 

		reset_matrix_cache(modelparams)
		transitionrate_counts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.01
		transitionrate_totals = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.001
		for h=1:modelparams.numhiddenstates
			transitionrate_counts[h,h] = 0.0
		end
		transitionrate_events = ones(Float64, modelparams.numhiddenstates)*0.01
		transitionrate_times = ones(Float64, modelparams.numhiddenstates)*0.001

		aatransitionrate_events = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		aatransitionrate_times = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.001

		aatransitionrate_counts = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		aatransitionrate_totals = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		for aa=1:modelparams.alphabet
			aatransitionrate_counts[aa,aa] = 0.0
		end

		totalbranchlength_output,accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total,totalhiddentime,totalaatime = sampletraininginstances(iter, rng, trainingexamples, modelparams, maxsamplesperiter=maxsamplesperiter, familyiter=familyiter, sitethreshold=sitethreshold, dosamplesiterates=dosamplesiterates, samplebranchlengths=samplebranchlengths, family_names=family_names)

		if !independentsites
			estimate_hidden_transition_probs(modelparams)
		else
			estimate_hidden_mixture(modelparams)
		end

		for h=1:modelparams.numhiddenstates			
			estimate_categorical(modelparams.hiddennodes[h].aa_node, 1.0)		
			estimatevonmises(modelparams.hiddennodes[h].phi_node)
			estimatevonmises(modelparams.hiddennodes[h].psi_node)
			estimatevonmises(modelparams.hiddennodes[h].omega_node)
			estimate_bvm(modelparams.hiddennodes[h].phipsi_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle1_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle2_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle3_node)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].phipsi_node.mu,"\t",modelparams.hiddennodes[h].phipsi_node.k, "\t", modelparams.hiddennodes[h].phi_node.mu, "\t", modelparams.hiddennodes[h].psi_node.mu, "\t", modelparams.hiddennodes[h].phi_node.kappa, "\t", modelparams.hiddennodes[h].psi_node.kappa)		
			println(modelparams.hiddennodes[h].phipsi_node.count,"\t",modelparams.hiddennodes[h].phi_node.N,"\t",modelparams.hiddennodes[h].psi_node.N)
			#println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle1_node.mu,"\t",modelparams.hiddennodes[h].bond_angle1_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle1_node.N)
			#println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle2_node.mu,"\t",modelparams.hiddennodes[h].bond_angle2_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle2_node.N)
			#println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle3_node.mu,"\t",modelparams.hiddennodes[h].bond_angle3_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle3_node.N)
			#println(iter,"\t",h,"\t", modelparams.hiddennodes[h].bond_lengths_node.mvn.μ)
			estimate_multivariate_node(modelparams.hiddennodes[h].bond_lengths_node)
			for aa=1:modelparams.alphabet
				estimatevonmises(modelparams.hiddennodes[h].phi_nodes[aa])
				estimatevonmises(modelparams.hiddennodes[h].psi_nodes[aa])
				estimatevonmises(modelparams.hiddennodes[h].omega_nodes[aa])
				estimate_bvm(modelparams.hiddennodes[h].phipsi_nodes[aa])
				if modelparams.hidden_conditional_on_aa
					print(iter,"\t",h,"\t",aminoacids[aa],"\t",@sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]),"\t",modelparams.hiddennodes[h].phi_nodes[aa].mu,"\t",modelparams.hiddennodes[h].phi_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].phi_nodes[aa].N)		
					println("\t",modelparams.hiddennodes[h].psi_nodes[aa].mu,"\t",modelparams.hiddennodes[h].psi_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].psi_nodes[aa].N)
					#println(iter,"\t",h,"\t",aminoacids[aa],"\t",@sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]),"\t",modelparams.hiddennodes[h].omega_nodes[aa].mu,"\t",modelparams.hiddennodes[h].omega_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].omega_nodes[aa].N)						end
				end
			end
			for aa=1:20
				println(iter,"\t",h,"\t",aminoacids[aa],"\t", @sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]))
			end			
		end
		

		#=
		if learnrates
			rate_cat_events = ones(Float64, modelparams.numrates)*0.01
			rate_cat_times = ones(Float64, modelparams.numrates)*0.001	
			for (proteins,nodelist,json_family,sequences) in trainingexamples
				numcols = length(proteins[1])
				inputcols =  Int[col for col=1:numcols]
				for col in inputcols
					cols = Int[]
					selcol = 1
					if col > 1
						push!(cols,col-1)
						selcol = 2
					end
					push!(cols,col)
					if col < numcols
						push!(cols,col+1)
					end

					for node in nodelist
						if !isroot(node)
							
							multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col])])
							hiddeniter = multi_iter.branchpathiterators[1]
							aaiter = multi_iter.branchpathiterators[2]
							for it in multi_iter
								dt = (multi_iter.currtime-multi_iter.prevtime)
								
								changecol = selcol
								rate_entry = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[1], aaiter.prevstates[1], node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
								rate_cat_times[node.data.ratesbranchpath.paths[col][end]] += -rate_entry*dt*node.branchlength/modelparams.rates[node.data.ratesbranchpath.paths[col][end]]

								changecol = 0
								if multi_iter.branchpathindex == 1
									changecol = hiddeniter.mincol
								else
									changecol = aaiter.mincol
								end
								if multi_iter.branchpathindex == 1 && hiddeniter.mincol == selcol
									rate_cat_events[node.data.ratesbranchpath.paths[col][end]] += 1.0	
								end
								if multi_iter.branchpathindex == 2
									rate_cat_events[node.data.ratesbranchpath.paths[col][end]] += 1.0
								end
							end
						end
					end
				end
			end
			#modelparams.rates = rate_cat_events./rate_cat_times
			for r=1:modelparams.numrates
				alpha = rate_cat_events[r]+1.0
				beta = rate_cat_times[r]
				dist = Gamma(alpha, 1.0/beta)
				modelparams.rates[r] = rand(dist)
			end
			println("RATES ", modelparams.rates)
		end=#

		println("Acceptance:\t",accepted_hidden/accepted_hidden_total,"\t",accepted_aa/accepted_aa_total)

		for training_example in trainingexamples
			for node in training_example[2]
				if !isroot(node)
					totalbranchlength_output += node.branchlength
				end
			end
		end
		modelparams.branchscalingfactor = totalbranchlength_output/totalbranchlength_input

		aahiddentransitionsrates = zeros(Float64, modelparams.numhiddenstates)
		aahiddentransition_events = zeros(Float64, modelparams.numhiddenstates)
		if learnrates
			if hiddenaascaling
				for (proteins,nodelist,json_family,sequences) in trainingexamples
					numcols = length(proteins[1])
					inputcols =  Int[col for col=1:numcols]
					for col in inputcols
						cols = Int[col]
						selcol = 1

						for node in nodelist
							if !isroot(node)
								
								multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col])])
								hiddeniter = multi_iter.branchpathiterators[1]
								aaiter = multi_iter.branchpathiterators[2]
								for it in multi_iter
									dt = (multi_iter.currtime-multi_iter.prevtime)
									
									changecol = selcol

									for aa=1:modelparams.alphabet
										thisaa = aaiter.prevstates[1]
										if thisaa != aa
											aaentry = getaaentry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], thisaa, aa, node.data.ratesbranchpath.paths[col][end])
											h = hiddeniter.prevstates[changecol]
											aahiddentransitionsrates[h] += aaentry*dt*node.branchlength/modelparams.aarates[h]
										end
									end

									changecol = 0
									if multi_iter.branchpathindex == 1
										changecol = hiddeniter.mincol
									else
										changecol = aaiter.mincol
									end
									if multi_iter.branchpathindex == 1 && hiddeniter.mincol == selcol
										
									else multi_iter.branchpathindex == 2
										h = hiddeniter.prevstates[1]
										aahiddentransition_events[h] += 1.0
									end
								end
							end
						end
					end
				end
				modelparams.aarates = aahiddentransitionsrates./aahiddentransition_events
				modelparams.aarates /= mean(modelparams.aarates)
			end
			println("AARATES", modelparams.aarates)

			aatransitions = 0
			hiddentransitions = 0
			for (proteins,nodelist,json_family,sequences) in trainingexamples
				numcols = length(proteins[1])
				inputcols =  Int[col for col=1:numcols]
				for col in inputcols
					cols = Int[]
					selcol = 1
					if col > 1
						push!(cols,col-1)
						selcol = 2
					end
					push!(cols,col)
					if col < numcols
						push!(cols,col+1)
					end

					for node in nodelist
						if !isroot(node)
							
							multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col])])
							hiddeniter = multi_iter.branchpathiterators[1]
							aaiter = multi_iter.branchpathiterators[2]
							for it in multi_iter
								dt = (multi_iter.currtime-multi_iter.prevtime)
								
								changecol = selcol
								for h=1:modelparams.numhiddenstates
									thish = hiddeniter.prevstates[changecol]
									if thish != h
										hiddenrateentry = gethiddenentry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), thish, h, aaiter.prevstates[1], node.data.ratesbranchpath.paths[col][end])
										#transitionrate_times[thish] += hiddenrateentry*dt*node.branchlength/modelparams.transitionrates[thish,h]
										transitionrate_totals[thish, h] += hiddenrateentry*dt*node.branchlength/modelparams.transitionrates[thish,h]
									end
								end

								for aa=1:modelparams.alphabet
									thisaa = aaiter.prevstates[1]
									if thisaa != aa
										aaentry = getaaentry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], thisaa, aa, node.data.ratesbranchpath.paths[col][end])
										aatransitionrate_times[thisaa, aa] += aaentry*dt*node.branchlength/modelparams.aa_exchangeablities[thisaa, aa]
									end
								end

								changecol = 0
								if multi_iter.branchpathindex == 1
									changecol = hiddeniter.mincol
								else
									changecol = aaiter.mincol
								end
								if multi_iter.branchpathindex == 1 && hiddeniter.mincol == selcol
									transitionrate_counts[hiddeniter.prevstates[changecol], hiddeniter.currstates[changecol]] += 1.0
									transitionrate_events[hiddeniter.prevstates[changecol]] += 1.0
									hiddentransitions += 1
								else multi_iter.branchpathindex == 2
									aatransitionrate_events[aaiter.prevstates[1], aaiter.currstates[1]] += 1.0
									aatransitions += 1
								end
							end
						end
					end
				end
			end

			for h=1:modelparams.numhiddenstates
				modelparams.transitionrates[h,h] = 0.0			
				for h2=1:modelparams.numhiddenstates
					if h != h2
						modelparams.transitionrates[h,h2] = (transitionrate_counts[h,h2]+transitionrate_counts[h2,h])/(transitionrate_totals[h,h2]+transitionrate_totals[h2,h])
						modelparams.transitionrates[h,h] -= modelparams.transitionrates[h,h2]
					end
				end
			end

			println("HQ", modelparams.transitionrates)	
			aatransitionrates = zeros(Float64, modelparams.alphabet, modelparams.alphabet)
			for aa1=1:modelparams.alphabet
				aatransitionrates[aa1,aa1] = 0.0
				for aa2=1:modelparams.alphabet
					if aa1 != aa2
						aatransitionrates[aa1,aa2] = (aatransitionrate_events[aa1,aa2]+aatransitionrate_events[aa2,aa1])/(aatransitionrate_times[aa1,aa2]+aatransitionrate_times[aa2,aa1])
						aatransitionrates[aa1,aa1] -= aatransitionrates[aa1,aa2]
					end
				end
			end	
			println("AQ",aatransitionrates)		
			println("LG",[sum(LGexchangeability[aa,:]) for aa=1:20])
			if !uselgrates
				modelparams.aa_exchangeablities = aatransitionrates
			end

			#modelparams.rates = rate_cat_events./rate_cat_times
			#println("RATES ", modelparams.rates)
		end
		#optimize_gamma(rng, trainingexamples, modelparams)

		#println("LL ", rates_likelihood(trainingexamples, modelparams))
		#optimizerates3(trainingexamples,modelparams)
		
		#=
		if (iter) % 5 == 0
			optimizerates(trainingexamples, modelparams)
		end
		=#
		println("Hidden transitions: ", hiddentransitions)
		println("AA transitions: ", aatransitions)
		reset_matrix_cache(modelparams)
		for (proteins,nodelist,json_family,sequences) in trainingexamples
			for node in nodelist
				if !isroot(node)
					#node.branchlength /= modelparams.branchscalingfactor
				end
			end
		end
		#modelparams.aa_exchangeablities *= modelparams.branchscalingfactor
		scaleaarates = mean(modelparams.aarates)
		modelparams.aarates /= scaleaarates
		modelparams.aa_exchangeablities *= scaleaarates

		#modelparams.transitionrates *= modelparams.branchscalingfactor
		modelparams.branchscalingfactor = 1.0


		augmentedll_end,observationll_end, augmented_array, observation_array = calculateloglikelihood_array(modelparams, trainingexamples)
		totalll_end = augmentedll_end+observationll_end

		aaloglikelihood, sequence_array = aaintegrated_loglikelihood_array(rng, modelparams, trainingexamples)
		final_loglikelihood = aaloglikelihood+observationll_end

		
		print(familyoutfile, iter-1)
		print(familyoutfile,"\t",final_loglikelihood,"\t",aaloglikelihood,"\t",observationll_end,"\t",augmentedll_end)
		for (name,obsll,seqll) in zip(family_names, observation_array,sequence_array)
			print(familyoutfile, "\t",obsll+seqll,"\t",obsll,"\t",seqll)
		end
		println(familyoutfile)
		flush(familyoutfile)

		#println(logwriter, "iter\tll\taall\tstructurell\tpathll")
		#println(logwriter, iter-1,"\t",totalll_end,"\t",augmentedll_end,"\t",integrated_loglikelihood,"\t",aaloglikelihood,"\t",observationll_end)
		println(logwriter, iter-1,"\t",final_loglikelihood,"\t",aaloglikelihood,"\t",observationll_end,"\t",augmentedll_end)
		flush(logwriter)
		if noimprovement >= 3 || final_loglikelihood > maxloglikelihood
			if final_loglikelihood > maxloglikelihood		
				fout = open(modelfile, "w")
				reset_matrix_cache(modelparams)
				Serialization.serialize(fout, modelparams)
				close(fout)
			end 
			noimprovement = 0
			maxloglikelihood = final_loglikelihood			
		else
			noimprovement += 1 
			fin = open(modelfile, "r")	
			modelparams = Serialization.deserialize(fin)
			close(fin)
		end
	end
end

function parse_training_commandline()
    settings = ArgParseSettings()
    settings.prog = prog
    settings.version = version
    settings.add_version = true

    add_arg_group(settings, "training")
    @add_arg_table settings begin
        "numhiddenstates"
        	help = "train a model with this number of hidden states"
         	arg_type = Int
         	required = true
         "--sitethreshold"
         	arg_type = Int
         	help = "the minimum number of successful site resamplings required at each protein site at each iteration (maximum number of attempts is MAXSAMPLESPERITER)"
         	default = 3
         "--maxsamplesperiter"
         	arg_type = Int
         	default = 500
         "--trainingdirs"
         	help = "comma-seperated list of directories that contain training instance files to be used for training"
        	arg_type = String
    	"--traininginstances"
		 	help = "train using only the first N training files"
		 	arg_type = Int
    	"--samplebranchlengths"
        	help = "sample branch lengths during training"
      	    arg_type = Bool
      	    default = true
      	"--samplesiterates"
         	help = ""
          	action = :store_true
  	 	"--angles-cond-aa"
         	help = "angles depend on amino acids"
          	arg_type = Int
          	default = 0
        "--ratemode"
        	arg_type = Int
         	default = 1
     	"--independentsites"
         	help = "ignores all neighbouring dependencies, just learns a mixture model"
          	action = :store_true
      	"--sequenceonly"
         	help = "use only amino acid characters as observations, ignores all structure observations"
          	action = :store_true
      	"--loadfromcache"
         	help = ""
          	action = :store_true
      	"--uselgrates"
         	help = ""
          	action = :store_true
      	"--hiddenaascaling"
         	help = ""
      	    arg_type = Bool
      	    default = true
      	 "--precluster"
      	 	help = ""
          	action = :store_true

	        
    end
    return parse_args(settings)
end

parsed_args = parse_training_commandline()
train(parsed_args)