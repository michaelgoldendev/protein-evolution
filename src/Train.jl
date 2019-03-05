include("Main.jl")
using Printf

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
	

	modelparams = ModelParams(LGmatrix,numhiddenstates,1.0)	
	maxloglikelihood = -Inf
	noimprovement = 0

	outputmodelname = string("_h.",modelparams.numhiddenstates)
	if independentsites
		outputmodelname = string(outputmodelname,".independentsites")
	end
	modelfile = string("models/model", outputmodelname, ".model")
	tracefile = string("models/trace", outputmodelname,".log")

	if parsed_args["loadfromcache"]
		modelparams = Serialization.deserialize(open(modelfile, "r"))
	end
	modelparams.usestructureobservations = usestructureobservations

	trainingexamples = Tuple[]

	#family_directories = ["../data/curated/curated_rna/", "../data/selected_families/"]
	family_directories = ["../data/selected_families/", "../data/curated/curated_rna/"]
	for family_dir in family_directories
		family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))
		for family_file in family_files
			full_path = abspath(joinpath(family_dir, family_file))
			json_family = JSON.parse(open(full_path, "r"))
			if 1 <= length(json_family["proteins"]) <= 1e10
				training_example = training_example_from_json_family(rng, modelparams, json_family)		
				println(json_family["newick_tree"])
				if length(training_example[2][1].children) == 1
					root = training_example[2][1].children[1]
					root.parent = Nullable{TreeNode}()
					training_example = (training_example[1], TreeNode[root], training_example[3], training_example[4])
				end
				push!(trainingexamples, training_example)
				if parsed_args["maxtraininginstances"] != nothing && length(trainingexamples) == parsed_args["maxtraininginstances"]
					break
				end
			end
		end
		if parsed_args["maxtraininginstances"] != nothing && length(trainingexamples) == parsed_args["maxtraininginstances"]
			break
		end
	end

	totalbranchlength_input = 0.0
	for training_example in trainingexamples
		for node in training_example[2]
			if !isroot(node)
				totalbranchlength_input += node.branchlength
			end
		end
	end

	for i=1:1
		for (proteins,nodelist,json_family,sequences) in trainingexamples
			numcols = length(proteins[1])		
		
			randcols = shuffle(rng, Int[i for i=1:numcols])
			for col in randcols
				samplepaths(rng,col,proteins,nodelist, modelparams)
			end
		end
	end	

	logwriter = open(tracefile, "w")
	println(logwriter, "iter\tll\taall\tstructurell\tpathll")
	for iter=1:10000
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

		totalbranchlength_output = 0.0
		accepted_hidden = 0.0
		accepted_hidden_total = 0.0
		accepted_aa = 0.0
		accepted_aa_total = 0.0
		for (proteins,nodelist,json_family,sequences) in trainingexamples
			numcols = length(proteins[1])
			
			accepted = zeros(Int, numcols)			
			for i=1:maxsamplesperiter
				randcols = shuffle(rng, Int[i for i=1:numcols])
				for col in randcols
					if accepted[col] < sitethreshold || i % 20 == 0 || (col > 1 && accepted[col-1] < sitethreshold) || (col < numcols && accepted[col+1] < sitethreshold) 
						a1,a2,a3,a4, hidden_accepted, aa_accepted = samplepaths(rng,col,proteins,nodelist, modelparams)
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
					println("min_accepted ", min_accepted," out of ", i, " mean is ", mean(accepted))
					break
				end


				if samplebranchlengths && (i <= sitethreshold || i % 10 == 0)
					for node in nodelist
						if !isroot(node)
							t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
							node.branchlength = t
							events = mean([length(p)-1.0 for p in node.data.branchpath.paths])
							aa_events = mean([length(p)-1.0 for p in node.data.aabranchpath.paths])
							#println(node.nodeindex,"\t",node.data.inputbranchlength,"\t",node.branchlength,"\t",events,"\t",aa_events)
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
							push!(modelparams.hiddennodes[h].phi_nodes[site.aa].data, site.phi)
							push!(modelparams.hiddennodes[h].omega_nodes[site.aa].data, site.omega)
							push!(modelparams.hiddennodes[h].psi_nodes[site.aa].data, site.psi)
							push!(modelparams.hiddennodes[h].bond_angle1_node.data, site.bond_angle1)
							push!(modelparams.hiddennodes[h].bond_angle2_node.data, site.bond_angle2)
							push!(modelparams.hiddennodes[h].bond_angle3_node.data, site.bond_angle3)
							add_point(modelparams.hiddennodes[h].bond_lengths_node, Float64[site.bond_length1, site.bond_length2, site.bond_length3])
						end
					end
				end
			end			
		end
		println("Acceptance:\t",accepted_hidden/accepted_hidden_total,"\t",accepted_aa/accepted_aa_total)

		for training_example in trainingexamples
			for node in training_example[2]
				if !isroot(node)
					totalbranchlength_output += node.branchlength
				end
			end
		end
		modelparams.branchscalingfactor = totalbranchlength_output/totalbranchlength_input

		if !independentsites
			estimate_hidden_transition_probs(modelparams)
		else
			estimate_hidden_mixture(modelparams)
		end

		for h=1:modelparams.numhiddenstates			
			estimate_categorical(modelparams.hiddennodes[h].aa_node, 1.0)
			for aa=1:modelparams.alphabet
				estimatevonmises(modelparams.hiddennodes[h].phi_nodes[aa])
				estimatevonmises(modelparams.hiddennodes[h].psi_nodes[aa])
				estimatevonmises(modelparams.hiddennodes[h].omega_nodes[aa])
				println(iter,"\t",h,"\t",aminoacids[aa],"\t",@sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]),"\t",modelparams.hiddennodes[h].phi_nodes[aa].mu,"\t",modelparams.hiddennodes[h].phi_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].phi_nodes[aa].N)		
				println(iter,"\t",h,"\t",aminoacids[aa],"\t",@sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]),"\t",modelparams.hiddennodes[h].psi_nodes[aa].mu,"\t",modelparams.hiddennodes[h].psi_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].psi_nodes[aa].N)
				println(iter,"\t",h,"\t",aminoacids[aa],"\t",@sprintf("%0.3f", modelparams.hiddennodes[h].aa_node.probs[aa]),"\t",modelparams.hiddennodes[h].omega_nodes[aa].mu,"\t",modelparams.hiddennodes[h].omega_nodes[aa].kappa,"\t",modelparams.hiddennodes[h].omega_nodes[aa].N)				
			end
			
			estimatevonmises(modelparams.hiddennodes[h].phi_node)
			estimatevonmises(modelparams.hiddennodes[h].psi_node)
			estimatevonmises(modelparams.hiddennodes[h].omega_node)			
			estimatevonmises(modelparams.hiddennodes[h].bond_angle1_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle2_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle3_node)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle1_node.mu,"\t",modelparams.hiddennodes[h].bond_angle1_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle1_node.N)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle2_node.mu,"\t",modelparams.hiddennodes[h].bond_angle2_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle2_node.N)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle3_node.mu,"\t",modelparams.hiddennodes[h].bond_angle3_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle3_node.N)
			println(iter,"\t",h,"\t", modelparams.hiddennodes[h].bond_lengths_node.mvn.μ)
			estimate_multivariate_node(modelparams.hiddennodes[h].bond_lengths_node)
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

		aatransitions = 0
		hiddentransitions = 0
		if learnrates
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
			modelparams.aa_exchangeablities = aatransitionrates

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
					node.branchlength /= modelparams.branchscalingfactor
				end
			end
		end
		modelparams.aa_exchangeablities *= modelparams.branchscalingfactor
		modelparams.transitionrates *= modelparams.branchscalingfactor
		modelparams.branchscalingfactor = 1.0

		augmentedll_end,observationll_end = calculateloglikelihood(modelparams, trainingexamples)
		totalll_end = augmentedll_end+observationll_end

		aaloglikelihood = aaintegrated_loglikelihood(rng, modelparams, trainingexamples)
		final_loglikelihood = aaloglikelihood+observationll_end

		#println(logwriter, "iter\tll\taall\tstructurell\tpathll")
		#println(logwriter, iter-1,"\t",totalll_end,"\t",augmentedll_end,"\t",integrated_loglikelihood,"\t",aaloglikelihood,"\t",observationll_end)
		println(logwriter, iter-1,"\t",final_loglikelihood,"\t",aaloglikelihood,"\t",observationll_end,"\t",augmentedll_end)
		flush(logwriter)
		if noimprovement >= Inf || final_loglikelihood > maxloglikelihood
			noimprovement = 0
			maxloglikelihood = final_loglikelihood			
			fout = open(modelfile, "w")
			reset_matrix_cache(modelparams)
			Serialization.serialize(fout, modelparams)
			close(fout)
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
    	"--maxtraininginstances"
		 	help = "train using only the first N training files"
		 	arg_type = Int
    	"--samplebranchlengths"
        	help = "sample branch lengths during training"
      	    arg_type = Bool
      	    default = true
     	"--independentsites"
         	help = "ignores all neighbouring dependencies, just learns a mixture model"
          	action = :store_true
      	"--sequenceonly"
         	help = "use only amino acid characters as observations, ignores all structure observations"
          	action = :store_true
      	"--loadfromcache"
         	help = ""
          	action = :store_true
	        
    end
    return parse_args(settings)
end

parsed_args = parse_training_commandline()
train(parsed_args)