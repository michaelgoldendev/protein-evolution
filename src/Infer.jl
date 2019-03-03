include("Main.jl")

using TraitAssociation

function infer(parsed_args=Dict{String,Any}())
	rng = MersenneTwister(58498012421121)
	Random.seed!(12345)

	modelfile = parsed_args["inference_model_file"] != nothing ? parsed_args["inference_model_file"] : "model_h.5_learnrates.true.simultaneous.scaling.model"
	println("Using model: ", modelfile)

	fin = open(modelfile, "r")	
	modelparams = Serialization.deserialize(fin)
	close(fin)

	#samplebranchlengths = parsed_args["samplebranchlengths"]
	samplebranchlengths = false
	dosamplesiterates = true

	blindscores = Dict{String,Array{Float64,1}}()
	blindnodenames = parsed_args["blindproteins"] != nothing ? split(parsed_args["blindproteins"], r",|;") : String[]
	if length(blindnodenames) > 0
		println("Blinding protein sequences and structures: ", join(blindnodenames, ", "))
	end

	modelparams.scalingfactor = modelparams.branchscalingfactor
	println("Scaling factor: ", modelparams.scalingfactor)

	#=
	family_dir = "../data/families/"
	family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))
	family_file = joinpath(family_dir, family_files[11])
	
	#family_file = "../data/families/alpha-amylase_NC.ali_23.fam"

	full_path = abspath(family_file)
	json_family = JSON.parse(open(full_path, "r"))
	training_example = training_example_from_json_family(rng, modelparams, json_family, blindnodenames=blindnodenames, scalebranchlengths=samplebranchlengths)		
	println(json_family["newick_tree"])
	if length(training_example[2][1].children) == 1
		root = training_example[2][1].children[1]
		root.parent = Nullable{TreeNode}()
		training_example = (training_example[1], TreeNode[root], training_example[3], training_example[4])
	end
	proteins,nodelist,json_family,sequences = training_example=#
	
	#proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/alignments/HCV_REF_2014_ns5b_PRO_curated.fasta"), blindnodenames=blindnodenames)
	#proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/alignments/hiv-curated-sel.fasta"), blindnodenames=blindnodenames)

	#proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/influenza_a/HA/H1N1_selection2.fas"), newickfile=abspath("../data/influenza_a/HA/H1N1_selection2_rooted.fas.nwk"), blindnodenames=blindnodenames)
	#proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/influenza_a/HA/selection3.fasta"), newickfile=abspath("../data/influenza_a/HA/selection3.fasta.nwk"), blindnodenames=blindnodenames)
	#proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/influenza_a/HA/selection3.fasta"), newickfile=abspath("tree.mean.branch.consensus.nwk"), blindnodenames=blindnodenames)


	#proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/influenza_a/HA/selection4.fasta"), newickfile=abspath("tree.mean.branch.consensus.nwk"), blindnodenames=blindnodenames)

	#proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/hiv/curated6.fasta"), newickfile=abspath("../data/hiv/curated6.nwk"), blindnodenames=blindnodenames)
	
	
	fastafile = abspath("../data/influenza_a/HA/selection3.fasta")
	newickfile= abspath("../data/influenza_a/HA/selection3.fasta.nwk")	
	#newickfile= abspath("tree.mean.branch.consensus.nwk")		
	blindnodenames = String["6n41.pdb_1951"]
	#blindnodenames = String[]

	
	#fastafile = abspath("../data/hiv/curated6.fasta")
	#newickfile=abspath("../data/hiv/curated6.nwk")	
	#newickfile=abspath("tree.mean.branch.consensus.nwk")	
	#blindnodenames = String["B.US.1978.SF4.KJ704795"]
	#blindnodenames = String[]

	#fastafile = abspath("../data/test_data/hiv_pol_selection.fasta")
	#newickfile=abspath("../data/test_data/hiv_pol_selection.fasta.nwk")	
	#blindnodenames = String["CPZ.GA.1988.GAB1.X52154"]

	#fastafile = abspath("../data/test_data/westnile_dengue_selection.fasta")
	#newickfile=abspath("../data/test_data/westnile_dengue_selection.fasta.nwk")	
	#blindnodenames = String["AEN02430.1"]

	proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, fastafile, newickfile=newickfile, blindnodenames=blindnodenames)
	
	if length(blindnodenames) > 0
		LGreconstruction_score = TraitAssociation.pathreconstruction(fastafile,newickfile,blindnodenames,3)
		println("LGreconstruction_score ",LGreconstruction_score)
	end



	outputprefix = string("output/",basename(modelfile),".",basename(fastafile))

	println("Initialisation finished.")

	inputnodelist = deepcopy(nodelist)
	inputtreewriter = open("$(outputprefix).input.nwk", "w")
	println(inputtreewriter, getnewick(nodelist[1]))
	close(inputtreewriter)

	

	numcols = length(proteins[1])
	mcmcwriter = open("$(outputprefix).log", "w")
	treewriter = open("$(outputprefix).mcmc.log", "w")
	print(mcmcwriter, "iter\ttotalll\tpathll\tsequencell\tobservationll\tscalingfactor")
	for blindnodename in blindnodenames
		print(mcmcwriter,"\t",blindnodename)
	end
	for node in nodelist
		if !isroot(node)
			if node.name == ""
				print(mcmcwriter,"\tbranchlength_node$(node.nodeindex)")
			else 
				print(mcmcwriter,"\tbranchlength_$(node.name)")
			end
		end
	end
	for node in nodelist
		if !isroot(node)
			if node.name == ""
				print(mcmcwriter,"\thidden_events_node$(node.nodeindex)")
			else 
				print(mcmcwriter,"\thidden_events_$(node.name)")
			end
		end
	end
	for node in nodelist
		if !isroot(node)
			if node.name == ""
				print(mcmcwriter,"\taa_events_node$(node.nodeindex)")
			else 
				print(mcmcwriter,"\taa_events_$(node.name)")
			end
		end
	end
	for node in nodelist
		if !isroot(node)
			if node.name == ""
				print(mcmcwriter,"\tsubs_per_site_node$(node.nodeindex)")
			else 
				print(mcmcwriter,"\tsubs_per_site_$(node.name)")
			end
		end
	end
	println(mcmcwriter)
	for node in nodelist
		if node.name in blindnodenames
			node.data.fixbranchlength = true
			if !isnull(node.parent)
				get(node.parent).data.fixbranchlength = true
			end
		end
	end

	subspersite_cache = Dict{Int,Array{Float64,1}}()
	branchlength_cache = Dict{Int,Array{Float64,1}}()
	accepted_alpha = 0.0
	accepted_alpha_total = 0.0	
	reset_matrix_cache(modelparams)
	count_hidden_acceptance = zeros(Int, numcols)
	count_hidden_total = zeros(Int, numcols)
	for iter=1:10000

		if samplebranchlengths
			println("$(iter).1 Sampling branch lengths START")
			randindex = rand(2:length(nodelist))
			for node in nodelist
				if !isroot(node) && !node.data.fixbranchlength	
					t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
					#println(iter," branches ", node.branchlength,"\t",t)				
					#oldll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
					node.branchlength = t
					#newll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
					#println(iter," likeliho ", newll,"\t",oldll,"\t",propratio,"\t",newll-oldll+propratio)
				end
			end
			println("$(iter).1 Sampling branch lengths DONE")
		elseif iter > 10 && parsed_args["samplescalingfactor"]
			proposescalingfactor(rng, nodelist, Int[col for col=1:numcols], modelparams)
			reset_matrix_cache(modelparams)
		end
		
		println("$(iter).2 Sampling sites START")
		randcols = shuffle(rng, Int[i for i=1:numcols])
		for col in randcols
			
			#if iter % 2 == 0
			a1,a2,a3,a4, accepted_hidden, accepted_aa = samplepaths(rng, col, proteins,nodelist, modelparams, dosamplesiterates=dosamplesiterates)
			if accepted_hidden
				count_hidden_acceptance[col] += 1.0
			end			
			count_hidden_total[col] += 1.0
			#else
				#samplepaths(rng, col, proteins,nodelist, modelparams)
			#eend
			#samplepaths(rng, col, proteins,nodelist, modelparams)
			#samplepaths_seperate(rng, col, proteins,nodelist, modelparams)
			#felsensteinresample(rng, proteins, nodelist, col, cols, col, modelparams)
			#felsensteinresample_aa(rng, proteins, nodelist, col, cols,col, modelparams)
		end
		println(sort(count_hidden_acceptance./count_hidden_total))
		println(mean(count_hidden_acceptance./count_hidden_total))
		println("$(iter).2 Sampling sites DONE")

		println("$(iter).3 Sampling blind nodes START")
		@time begin
			for selnode in nodelist
				if selnode.name in blindnodenames					
					countmatches = 0.0
					counttotal = 0.0
					println(selnode.name)		
					alignedsequence = sequences[selnode.seqindex]
					sequence, phi_psi, omega, bond_angles, bond_lengths = protein_to_lists(sampletreenode(rng, selnode, modelparams, alignedsequence))
				    inputsequence = replace(alignedsequence, "-" => "")
				    for (aa1,aa2) in zip(sequence,inputsequence)
				    	if aa1 == aa2
				    		countmatches += 1.0
				    	end
				    	counttotal += 1.0
				    	#println(iter,"\t",aa1,"\t",aa2,"\t",countmatches/counttotal)
				    end
				    blindscoresarray = get(blindscores, selnode.name, Float64[])
				    push!(blindscoresarray, countmatches/counttotal)
				    blindscores[selnode.name] = blindscoresarray
				    println(blindscoresarray)
				    println(mean(blindscoresarray[max(1,div(length(blindscoresarray),2)):end]))
				end
			end
		end
		
		println("$(iter).3 Sampling blind nodes DONE")

		println("$(iter).4 Calculating likelihoods START")
		augmentedll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)	
		sequencell = 0.0
		for col=1:numcols
			sequencell += felsensteinresample_aa(rng, proteins, nodelist, Int[col], col, modelparams, false)
		end
		observationll = observationloglikelihood(proteins, nodelist, modelparams)
		print(mcmcwriter, iter-1,"\t",augmentedll+observationll,"\t",augmentedll,"\t",sequencell,"\t",observationll,"\t",modelparams.scalingfactor)
		for blindnodename in blindnodenames
			print(mcmcwriter,"\t",blindscores[blindnodename][end])
		end
		for node in nodelist
			if !isroot(node)				
				cached_branchlengths = get(branchlength_cache, node.nodeindex, Float64[])
				push!(cached_branchlengths, node.branchlength)
				branchlength_cache[node.nodeindex] = cached_branchlengths

				print(mcmcwriter,"\t$(node.branchlength)")
			end
		end
		println("$(iter).4 Calculating likelihoods DONE")	

		if dosamplesiterates
			currentll = augmentedll
			#=
			currentll = 0.0
			for col=1:numcols
				currentll += felsensteinresample_aa(rng, proteins, nodelist, Int[col], col, modelparams, false)
			end
			=#
			currentalpha = modelparams.rate_alpha
			modelparams.rate_alpha += randn(rng)*0.75
			if modelparams.rate_alpha > 0.0
				modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)
				reset_matrix_cache(modelparams)
				proposedll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
				#=
				proposedll = 0.0
				for col=1:numcols
					proposedll += felsensteinresample_aa(rng, proteins, nodelist, Int[col], col, modelparams, false)
				end=#

				if exp(proposedll-currentll) > rand(rng)
					augmentedll = proposedll
					accepted_alpha += 1.0
				else
					modelparams.rate_alpha = currentalpha
					modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)
					reset_matrix_cache(modelparams)
				end
			else
				modelparams.rate_alpha = currentalpha
				modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)
			end
			accepted_alpha_total += 1.0
			println(modelparams.rate_alpha,"\t",modelparams.rates,"\t", accepted_alpha/accepted_alpha_total)
		end
	
		for node in nodelist
			if !isroot(node)
				print(mcmcwriter,"\t$(sum([length(path)-1 for path in node.data.branchpath.paths]))")
			end
		end
		for node in nodelist
			if !isroot(node)
				#print(mcmcwriter,"\t$(sum([length(path)-1 for path in node.data.branchpath.paths]))")
				print(mcmcwriter,"\t",count_aminoacid_substitutions(rng,modelparams,node))
			end
		end

		outputnodelist = deepcopy(nodelist)
		for (node,outputnode) in zip(nodelist, outputnodelist)
			if !isroot(node)								
				subs_per_site = count_aminoacid_substitutions(rng,modelparams,node)/numcols
				outputnode.branchlength = subs_per_site				
				print(mcmcwriter,"\t$(subs_per_site)")
				cached_branchlengths = get(subspersite_cache, outputnode.nodeindex, Float64[])
				push!(cached_branchlengths, outputnode.branchlength)
				subspersite_cache[outputnode.nodeindex] = cached_branchlengths
			end
		end
		println(mcmcwriter)
		flush(mcmcwriter)

		println(treewriter, getnewick(outputnodelist[1]))
		flush(treewriter)

		consensustreewriter = open("$(outputprefix).mean.branch.consensus.nwk", "w")
		for outputnode in outputnodelist
			if !isroot(outputnode)
				outputnode.branchlength = mean(branchlength_cache[outputnode.nodeindex][max(1,div(length(branchlength_cache[outputnode.nodeindex]),2)):end])
			end
		end
		println(consensustreewriter, getnewick(outputnodelist[1]))
		close(consensustreewriter)

		consensustreewriter = open("$(outputprefix).median.branch.consensus.nwk", "w")
		for outputnode in outputnodelist
			if !isroot(outputnode)
				outputnode.branchlength = median(branchlength_cache[outputnode.nodeindex][max(1,div(length(branchlength_cache[outputnode.nodeindex]),2)):end])
			end
		end
		println(consensustreewriter, getnewick(outputnodelist[1]))
		close(consensustreewriter)


		consensustreewriter = open("$(outputprefix).mean.consensus.nwk", "w")
		for outputnode in outputnodelist
			if !isroot(outputnode)
				outputnode.branchlength = mean(subspersite_cache[outputnode.nodeindex][max(1,div(length(subspersite_cache[outputnode.nodeindex]),2)):end])
			end
		end
		println(consensustreewriter, getnewick(outputnodelist[1]))
		close(consensustreewriter)

		fout = open("$(outputprefix).branchlengths.csv", "w")
		println(fout, "node,input,branchlength,aa_subs_per_branchsite")
		for (inputnode, outputnode) in zip(inputnodelist, outputnodelist)
			if !isroot(inputnode)
				if inputnode.name == ""
					print(fout,"$(inputnode.nodeindex),")
				else
					print(fout,"$(inputnode.name),")
				end
				b1 = inputnode.branchlength
				b2 = mean(branchlength_cache[outputnode.nodeindex][max(1,div(length(branchlength_cache[outputnode.nodeindex]),2)):end])
				b3 = mean(subspersite_cache[outputnode.nodeindex][max(1,div(length(subspersite_cache[outputnode.nodeindex]),2)):end])
				println(fout, "$(b1),$(b2),$(b3)")
			end
		end
		close(fout)

		treewriter = open("$(outputprefix).scalings.nwk", "w")
		println(treewriter, getnewick(compare_branch_scalings(inputnodelist,outputnodelist)))
		close(treewriter)

		consensustreewriter = open("$(outputprefix).median.consensus.nwk", "w")
		for outputnode in outputnodelist
			if !isroot(outputnode)
				outputnode.branchlength = median(subspersite_cache[outputnode.nodeindex][max(1,div(length(subspersite_cache[outputnode.nodeindex]),2)):end])
			end
		end
		println(consensustreewriter, getnewick(outputnodelist[1]))
		close(consensustreewriter)
	end
	close(mcmcwriter)
	close(treewriter)
end

function parse_inference_commandline()
    settings = ArgParseSettings()
    settings.prog = prog
    settings.version = version
    settings.add_version = true

    add_arg_group(settings, "inference (prediction)")
    @add_arg_table settings begin
        "inference_model_file"
        	help = "specify model file to use for inference"
          	arg_type = String
          	required = true
        "--alignment"
        	help = "alignment in fasta format"
        	arg_type = String
        "--tree"
        	help = "tree in newick format"
        	arg_type = String        
        "--blindproteins"
        	help = "comma-seperated list of protein names, whose sequences and structures will be blinded"
        	arg_type = String        
        "--samplebranchlengths"
        	help = ""
      	    arg_type = Bool
      	    default = true
        "--samplescalingfactor"
        	help = "keep relative branch lengths fixed, but sample tree scaling factor"
          	action = :store_true
    end


    return parse_args(settings)
end

parsed_args = parse_inference_commandline()
infer(parsed_args)