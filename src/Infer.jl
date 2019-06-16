include("Main.jl")
using StructurePlots

using TraitAssociation
using Dates

function transitionmatrix(sequences::Array{String,1}, curri::Int, nexti::Int)
	transprobs = zeros(Float64,20,20)
	for seq in sequences
		aa1 = indexof(string(seq[curri]), aminoacids)
		aa2 = indexof(string(seq[nexti]), aminoacids)
		transprobs[aa1,aa2] += 1.0
	end 
	for aa=1:20
		transprobs[aa,:] = transprobs[aa,:]/(sum(transprobs[aa,:]) .+ 1e-20)
	end
	return transprobs
end

function getconsensus(seqs::Array{String,1})
	majorityseq = ""
	for col=1:length(seqs[1])
		aacount = zeros(Int,20)
		for seq in seqs
			aacount[CommonUtils.indexof(string(seq[col]),aminoacids)] += 1
		end
		majorityseq = string(majorityseq, aminoacids[findmax(aacount)[2]])
	end
	return majorityseq
end

function viterbi(sequences::Array{String,1})
	numcols = length(sequences[1])
	liks = zeros(Float64,numcols,20)
	for seq in sequences
		liks[1,indexof(string(seq[1]), aminoacids)] += 1.0
	end
	liks[1,:] = log.(liks[1,:]/sum(liks[1,:]))
	maxindices = zeros(Int,numcols,20)

	for nexti=2:numcols
		curri = nexti-1
		transprobs = log.(transitionmatrix(sequences, curri, nexti))
		for nextaa=1:20
			v = liks[curri,:] .+ transprobs[:,nextaa]
			liks[nexti,nextaa], maxindices[nexti,nextaa] = findmax(v)
		end
	end
	#println(maxindices)
	lastaa = findmax(liks[numcols,:])[2]	
	#println(lastaa)
	maxseq = string(aminoacids[lastaa])
	for col=numcols:-1:2
		lastaa = maxindices[col,lastaa]
		maxseq = string(aminoacids[lastaa], maxseq)
		#println(col,"\t",lastaa)
	end

	#println(liks)
	#rintln(transitionmatrix(sequences, 1, 2))
	return maxseq
end

function calculatematch(realseq::String, predseq::String)
	countmatches = 0.0
	counttotal = 0.0
    for (aa1,aa2) in zip(predseq,realseq)
    	if uppercase(aa2) in aminoacids
	    	if aa1 == aa2
	    		countmatches += 1.0
	    	end
	    	counttotal += 1.0
	    end
    end
    println("COUNTTOTAL ", counttotal)
    return countmatches/counttotal
end

function infer(parsed_args=Dict{String,Any}())
	rng = MersenneTwister(58498012421121)
	Random.seed!(12345)

	modelfile = parsed_args["inference_model_file"] != nothing ? parsed_args["inference_model_file"] : "model_h.5_learnrates.true.simultaneous.scaling.model"
	println("Using model: ", modelfile)

	fin = open(modelfile, "r")	
	modelparams = Serialization.deserialize(fin)
	close(fin)

	samplebranchlengths = parsed_args["samplebranchlengths"]
	dosamplesiterates = parsed_args["samplesiterates"]

	sequencescores = Dict{String,Array{Float64,1}}()
	structurescores = Dict{String,Array{Float64,1}}()
	sequencescores_current = Dict{String,Float64}()
	blindstructures = parsed_args["blindstructures"] != nothing ? convert(Array{String,1}, split(parsed_args["blindstructures"], r",|;")) : String[]
	blindproteins = parsed_args["blindproteins"] != nothing ? convert(Array{String,1}, split(parsed_args["blindproteins"], r",|;")) : String[]
	println(blindstructures)
	println(blindproteins)
	if length(blindproteins) > 0
		println("Blinding protein sequences and structures: ", join(blindproteins, ", "))
	end

	modelparams.scalingfactor = modelparams.branchscalingfactor
	println("Scaling factor: ", modelparams.scalingfactor)

	datafile = parsed_args["dataset"]
	newickfile = parsed_args["tree"]

	modelparams.numrates = 1
	if dosamplesiterates
		modelparams.numrates = 20
	end 
	mu = 1.0
	alpha = 1.0
	stationaryprobs = (modelparams.transitionprobs^50)[1,:]
	aafreqs = zeros(Float64,20)
	for h=1:modelparams.numhiddenstates
		aafreqs = aafreqs .+ modelparams.hiddennodes[h].aa_node.probs*stationaryprobs[h]
	end
	R = zeros(20,20)
	for i=1:20
		for j=1:20
		    if i != j
		        R[i,j] = modelparams.aa_exchangeablities[i,j]*aafreqs[j]
		        R[i,i] -= R[i,j]
		    end
		end
	end


	nodenames = AbstractString[]
	if endswith(datafile, ".fam")
		json_family = JSON.parse(open(datafile, "r"))
		newickstring = json_family["newick_tree"]
		root = gettreefromnewick(newickstring)
		nodelist = getnodelist(root)
		for node in nodelist
			if node.name != ""
				push!(nodenames,node.name)
			end
		end
	elseif newickfile != nothing
		root = gettreefromnewick(readlines(open(newickfile,"r"))[1])
		nodelist = getnodelist(root)
		for node in nodelist
			if node.name != ""
				push!(nodenames,node.name)
			end
		end
	end

	if parsed_args["random"]
		for nodename in nodenames
			if !(nodename in blindproteins)
				push!(blindproteins, nodename)
			end
		end
	elseif parsed_args["sequencesonly"]
		for nodename in nodenames
			if !(nodename in blindstructures)
				push!(blindstructures, nodename)
			end
		end
	elseif parsed_args["unblindproteins"] != nothing
		blindexcept = convert(Array{String,1}, split(parsed_args["unblindproteins"], r",|;"))
		for nodename in nodenames
			if !(nodename in blindproteins) && !(nodename in blindexcept)
				push!(blindstructures, nodename)
			end
		end	
	end
	println("blindproteins: ", blindproteins)
	println("blindstructures: ", blindstructures)

	proteins = nothing
	nodelist = nothing
	json_family = nothing
	sequences = nothing
	fastafile = datafile
	if endswith(datafile, ".fam")
		json_family = JSON.parse(open(datafile, "r"))
		proteins, nodelist,json_family, sequences = training_example_from_json_family(rng, modelparams, json_family, blindproteins=blindproteins, blindstructures=blindstructures)
		newickfile = string(tempname(),".nwk")
		fout = open(newickfile, "w")
		println(fout, json_family["newick_tree"])
		close(fout)

		fastafile = string(tempname(),".fas")
		println(fastafile)
		fout = open(fastafile, "w")
		for p in json_family["proteins"]
			println(fout, ">", p["name"])
			println(fout, p["aligned_sequence"])
		end
		close(fout)
	else
		proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, datafile, newickfile=newickfile, blindproteins=blindproteins)
	end

	mu = 1.0
	alpha = 1.0
	if dosamplesiterates
		mu,alpha = TraitAssociation.find_mu_and_alpha(fastafile,newickfile,blindproteins,modelparams.numrates,R,aafreqs)
	end
	println("mu: ", mu)
	println("alpha: ", alpha)
	println("END")

	if length(proteins) == 1
		println(length(nodelist))
		println(json_family["newick_tree"])
		marginalprobs = forwardbackward(rng, nodelist[1], modelparams)
		println(marginalprobs)
		exit()
	end

	LGreconstruction_score = 0.0
	if length(blindproteins) > 0		
		#LGreconstruction_score,mu,alpha = TraitAssociation.pathreconstruction(datafile,newickfile,blindproteins,modelparams.numrates,R,aafreqs,mu,alpha)
		LGreconstruction_score,mu,alpha = TraitAssociation.pathreconstruction(fastafile,newickfile,blindproteins,modelparams.numrates)
		println("LGreconstruction_score ",LGreconstruction_score)
	end

	
	if dosamplesiterates
		modelparams.scalingfactor = mu
		modelparams.rate_alpha = alpha
		modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)
		modelparams.rate_freqs = ones(Float64,modelparams.numrates)/modelparams.numrates		
		reset_matrix_cache(modelparams)
	end

	structuresamples = Dict{String,Any}()
	for name in blindstructures
		structuresamples[name] = Sample(name, modelparams)
	end
	structuresamples["metadata:blindproteins"] = blindproteins
	structuresamples["metadata:blindstructures"] = blindstructures


	
	outdir = string("output/", Dates.format(Dates.now(), "yyyy-mm-dd.HH\\hMM\\mSS\\s.s"), ".", string(hash(parsed_args), base=36), "/")
	println(outdir)
	if !isdir(outdir)
		mkpath(outdir)
	end
	outputprefix = string(outdir,"output")
	fout = open(joinpath(outdir,"arguments.json"),"w")
	JSON.print(fout, parsed_args)
	close(fout)

	println("Initialisation finished.")

	inputnodelist = deepcopy(nodelist)
	inputtreewriter = open("$(outputprefix).input.nwk", "w")
	println(inputtreewriter, getnewick(nodelist[1]))
	close(inputtreewriter)

	

	numcols = length(proteins[1])
	mcmcwriter = open("$(outputprefix).log", "w")
	treewriter = open("$(outputprefix).mcmc.log", "w")
	print(mcmcwriter, "iter\ttotalll\tpathll\tsequencell\tobservationll\tscalingfactor\talpha\tsecondsperiter\tacceptancehidden\tacceptanceaa")

	blindseq_writers = Dict{String,Any}()
	for blindnodename in blindproteins
		print(mcmcwriter,"\t",blindnodename)
		try
			blindseq_writers[blindnodename] = open(string("$(outputprefix).", blindnodename,".fasta"), "w")
		catch 

		end
	end
	for blindnodename in blindproteins
		print(mcmcwriter,"\t",string("LGmean:", blindnodename))
		print(mcmcwriter,"\t",string("majorityseq:", blindnodename))
		print(mcmcwriter,"\t",string("viterbi:", blindnodename))
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
	for name in blindstructures
		print(mcmcwriter,"\t",string("angularrmsd:", name))
	end
	for name in blindstructures
		print(mcmcwriter,"\t",string("angularrmsd25:", name))
	end
	for name in blindstructures
		print(mcmcwriter,"\t",string("angularrmsd50:", name))
	end
	for name in blindstructures
		print(mcmcwriter,"\t",string("angularrmsd75:", name))
	end
	for name in blindstructures
		print(mcmcwriter,"\t",string("angularrmsd90:", name))
	end
	println(mcmcwriter)
	
	majority = Dict{String, Array{String,1}}()
	for node in nodelist
		if node.name in blindproteins
			majority[node.name] = String[]
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
	count_aa_acceptance = zeros(Int, numcols)

	maxaugmentedll = -Inf
	sequencescoresatmax = Float64[]
	ratetotals = zeros(Float64, numcols)
	maxiters = parsed_args["maxiters"]
	cache = OUDiffusionCache()
	for iter=1:maxiters
		starttime = time()
		if iter > 10 && samplebranchlengths
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
		end

		if iter > 10 && parsed_args["samplescalingfactor"]
			#augmentedll1 = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)	
			propratio = proposescalingfactor(rng, nodelist, Int[col for col=1:numcols], modelparams)
			#maxscalingfactor(rng, nodelist, Int[col for col=1:numcols], modelparams)
			reset_matrix_cache(modelparams)
			#augmentedll2 = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)	
			#println((augmentedll2-augmentedll1+propratio),"\t",augmentedll1,"\t", augmentedll2,"\t",propratio)
		end
		
		println("$(iter).2 Sampling sites START")
		randcols = shuffle(rng, Int[i for i=1:numcols])
		for col in randcols
			a1,a2,a3,a4, accepted_hidden, accepted_aa = samplepaths_seperate_new(rng, col,nodelist, modelparams, cache, dosamplesiterates=dosamplesiterates, accept_everything=(iter<=3))
			if accepted_hidden
				count_hidden_acceptance[col] += 1.0
			end			
			if accepted_aa
				count_aa_acceptance[col] += 1.0
			end
			count_hidden_total[col] += 1.0
		end
		println("min acceptance: ", minimum(count_hidden_acceptance./count_hidden_total))
		println("mean acceptance: ", mean(count_hidden_acceptance./count_hidden_total))
		println("$(iter).2 Sampling sites DONE")

		println("$(iter).3 Sampling blind nodes START")
		@time begin
			for selnode in nodelist
				if selnode.name in blindproteins					
					alignedsequence = sequences[selnode.seqindex]								
					sampledseq = ""
					for col=1:numcols
						if  alignedsequence[col] != '-'
							sampledseq = string(sampledseq, aminoacids[selnode.data.aabranchpath.paths[col][end]])
						end
					end
					println(selnode.name)		
					try
						println(blindseq_writers[selnode.name], ">sample$(iter)")
						println(blindseq_writers[selnode.name], sampledseq)
						flush(blindseq_writers[selnode.name])
					catch 

					end
					majorityseq = ""
					viterbiseq = ""		
					push!(majority[selnode.name], sampledseq)
					startindex = max(1, div(length(majority[selnode.name]),3))
					seqs = majority[selnode.name][startindex:end]
					
					for col=1:length(sampledseq)
						aacount = zeros(Int,20)
						for seq in seqs
							aacount[CommonUtils.indexof(string(seq[col]),aminoacids)] += 1
						end
						majorityseq = string(majorityseq, aminoacids[findmax(aacount)[2]])
					end
					viterbiseq = viterbi(seqs)

				    inputsequence = replace(alignedsequence, "-" => "")
				    sequencescoresarray = get(sequencescores, selnode.name, Float64[])
				    push!(sequencescoresarray, calculatematch(inputsequence,sampledseq))
				    sequencescores[selnode.name] = sequencescoresarray
				    sequencescores_current[string("majorityseq:",selnode.name)] = calculatematch(inputsequence,majorityseq)
				    sequencescores_current[string("viterbi:",selnode.name)] = calculatematch(inputsequence,viterbiseq)
				    println(sequencescoresarray)
				    println("mean: ", mean(sequencescoresarray[max(1,div(length(sequencescoresarray),3)):end]))
				    println("std: ", std(sequencescoresarray[max(1,div(length(sequencescoresarray),3)):end]))
				end
			end
		end
		
		println("$(iter).3 Sampling blind nodes DONE")

		println("$(iter).4 Calculating likelihoods START")
		augmentedll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)	
		sequencell = 0.0
		for col=1:numcols
			sequencell += felsensteinresample_aa(rng, nodelist, Int[col], col, modelparams, false)
		end
		observationll = observationloglikelihood(proteins, nodelist, modelparams)
		elapsedtime = time()-starttime 
		acceptance_rate_hidden =  mean(count_hidden_acceptance./count_hidden_total)
		acceptance_rate_aa =  mean(count_aa_acceptance./count_hidden_total)
		print(mcmcwriter, iter-1,"\t",augmentedll+observationll,"\t",augmentedll,"\t",sequencell,"\t",observationll,"\t",modelparams.scalingfactor, "\t", modelparams.rate_alpha, "\t",elapsedtime,"\t", acceptance_rate_hidden, "\t", acceptance_rate_aa)
		for blindnodename in blindproteins
			print(mcmcwriter,"\t",sequencescores[blindnodename][end])
		end
		for blindnodename in blindproteins
			print(mcmcwriter,"\t",LGreconstruction_score)
			print(mcmcwriter,"\t",sequencescores_current[string("majorityseq:",blindnodename)])
			print(mcmcwriter,"\t",sequencescores_current[string("viterbi:",blindnodename)])
		end
		if augmentedll > maxaugmentedll
			maxaugmentedll = augmentedll
			sequencescoresatmax = Float64[sequencescores[blindnodename][end] for blindnodename in blindproteins]
		end
		#println("maxll","\t",sequencescoresatmax,"\t", maxaugmentedll)

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
			#=
			for z=1:4
				currentll = augmentedll
				#=
				currentll = 0.0
				for col=1:numcols
					currentll += felsensteinresample_aa(rng, proteins, nodelist, Int[col], col, modelparams, false)
				end
				=#
				currentalpha = modelparams.rate_alpha
				modelparams.rate_alpha += randn(rng)*0.20
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
				println(modelparams.rate_alpha,"\t",modelparams.rates,"\tacceptance: ", accepted_alpha/accepted_alpha_total)

				for col=1:numcols
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
					samplesiterates(rng, cols, col, nodelist, modelparams)
				end
				augmentedll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
			end
			=#
		end

		
		for node in nodelist
			if !isleafnode(node) && node.name != ""			
				ancestorwriter = open("reconstruction_$(node.name).fas", iter == 1 ? "w" : "a")
				sampledseq = ""
				for col=1:numcols
					sampledseq = string(sampledseq, aminoacids[node.data.aabranchpath.paths[col][end]])
				end
				seqs = get(majority, node.name, String[])
				push!(seqs, sampledseq)
				majority[node.name] = seqs
				println(ancestorwriter, ">iter$(iter)")
				println(ancestorwriter, sampledseq)
				close(ancestorwriter)

				startindex = max(1, div(length(seqs),3))
				seqs = majority[node.name][startindex:end]

				consensusseqwriter = open("reconstruction_$(node.name).fas", "w")
				println(consensusseqwriter, ">consensus")
				println(consensusseqwriter, getconsensus(seqs))
				close(consensusseqwriter)
			end
		end
	
		for node in nodelist
			if !isroot(node)
				print(mcmcwriter,"\t$(sum([length(path)-1 for path in node.data.branchpath.paths]))")
			end
		end

		for node in nodelist
			if !isroot(node)
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

		if dosamplesiterates
			for col=1:numcols
				ratecat = nodelist[1].data.ratesbranchpath.paths[col][end]
				ratetotals[col] += modelparams.rates[ratecat]
			end
		end

		if iter % 1 == 0

			othernames = AbstractString[]
			for selnode in nodelist
				if selnode.seqindex > 0
					if json_family != nothing && haskey(json_family, "proteins") && haskey(json_family["proteins"][selnode.seqindex], "aligned_phi_psi")
						sequence, phi_psi, omega, bond_angles, bond_lengths = protein_to_lists(sampletreenode(rng,selnode,modelparams,sequences[selnode.seqindex]))
						phipsi_real = Tuple{Float64,Float64}[(p[1],p[2]) for p in json_family["proteins"][selnode.seqindex]["aligned_phi_psi"]]
						if count(x -> x[1] > -100.0 || x[2] > -100.0, phipsi_real) > 0 && !(selnode.name in blindstructures) && !(selnode.name in blindproteins)
							push!(othernames,selnode.name)
						end
						#angulardist = angular_rmsd(phi_psi, phipsi_real)
						#println("LENGTHS ", length(phi_psi), "\t", length(phipsi_real))
						angulardist = angular_rmsd_percentile(phi_psi, phipsi_real)
						
						#println("angular_rmsd ($(selnode.name)): ", angulardist)

						key = selnode.name
					    structurescoresarray = get(structurescores, key, Float64[])
					    push!(structurescoresarray,  angular_rmsd(phi_psi, phipsi_real))
						structurescores[key] = structurescoresarray

						key = string(selnode.name,"_perc25")
					    structurescoresarray = get(structurescores, key, Float64[])
					    push!(structurescoresarray, angular_rmsd_percentile(phi_psi, phipsi_real, 25.0))
						structurescores[key] = structurescoresarray

						key = string(selnode.name,"_perc50")
					    structurescoresarray = get(structurescores, key, Float64[])
					    push!(structurescoresarray, angular_rmsd_percentile(phi_psi, phipsi_real, 50.0))
						structurescores[key] = structurescoresarray

						key = string(selnode.name,"_perc75")
					    structurescoresarray = get(structurescores, key, Float64[])
					    push!(structurescoresarray, angular_rmsd_percentile(phi_psi, phipsi_real, 75.0))
						structurescores[key] = structurescoresarray

						key = string(selnode.name,"_perc90")
					    structurescoresarray = get(structurescores, key, Float64[])
					    push!(structurescoresarray, angular_rmsd_percentile(phi_psi, phipsi_real, 90.0))
						structurescores[key] = structurescoresarray

						if !haskey(structuresamples, selnode.name)
							structuresamples[selnode.name] = Sample(string(selnode.name), modelparams)
						end

						parent = get(selnode.parent)
						nodedata = Tuple[(copy(selnode.data.branchpath.paths[col]), copy(selnode.data.aabranchpath.paths[col]), parent.data.protein.sites[col].phi, parent.data.protein.sites[col].psi, selnode.data.protein.sites[col].phi, selnode.data.protein.sites[col].psi, selnode.branchlength) for col=1:numcols]
						#println(nodedata)
						push!(structuresamples[selnode.name].nodedata, nodedata)
						push!(structuresamples[selnode.name].phipsisamples, phi_psi)
						push!(structuresamples[selnode.name].aasamples, Int[selnode.data.aabranchpath.paths[col][end] for col=1:numcols])
						push!(structuresamples[selnode.name].hiddensamples, Int[selnode.data.branchpath.paths[col][end] for col=1:numcols])
						if json_family != nothing
							structuresamples[selnode.name].json_family = json_family
						end

						reset_matrix_cache(modelparams, false)
					end
				end
			end

			if iter % 25 == 1
				samplesfile = string("$(outputprefix).samples")

				fout = open(samplesfile, "w")
				Serialization.serialize(fout, structuresamples)
				close(fout)
				benchmarktype = ""
				if parsed_args["random"]
					benchmarktype = "Random"
				elseif parsed_args["sequencesonly"]
					benchmarktype = "Sequence only"
				elseif parsed_args["unblindproteins"] != nothing
					benchmarktype = "Homologous structure"
				end

				StructurePlots.plotaccuracy(samplesfile, othernames, benchmarktype)

				if iter > 1 && iter % 200 == 1
					println(othernames)
					StructurePlots.plotstructuresamples(samplesfile, othernames, benchmarktype)
				end
			end

		end
		for name in blindstructures
			print(mcmcwriter,"\t",structurescores[name][end])
		end	
		for name in blindstructures
			print(mcmcwriter,"\t",structurescores[string(name,"_perc25")][end])
		end
		for name in blindstructures
			print(mcmcwriter,"\t",structurescores[string(name,"_perc50")][end])
		end	
		for name in blindstructures
			print(mcmcwriter,"\t",structurescores[string(name,"_perc75")][end])
		end	
		for name in blindstructures
			print(mcmcwriter,"\t",structurescores[string(name,"_perc90")][end])
		end		
		println(mcmcwriter)
		flush(mcmcwriter)
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
        "--dataset"
        	help = "alignment in fasta format or .fam file"
        	arg_type = String
        "--tree"
        	help = "tree in newick format"
        	arg_type = String        
        "--blindproteins"
        	help = "comma-seperated list of protein names, whose sequences and structures will be blinded"
        	arg_type = String
    	"--blindstructures"
        	help = "comma-seperated list of protein names, whose structures only  will be blinded"
        	arg_type = String           
        "--samplebranchlengths"
        	help = ""
      	    action = :store_true
        "--samplescalingfactor"
        	help = "keep relative branch lengths fixed, but sample tree scaling factor"
          	action = :store_true
        "--samplesiterates"
        	help = ""
          	action = :store_true
        "--maxiters"
        	help = ""
        	arg_type = Int
        	default = typemax(Int)
        "--random"
        	help = ""
          	action = :store_true
        "--sequencesonly"
        	help = ""
        	action = :store_true
        "--unblindproteins"
        	help = ""
        	arg_type = String
    end


    return parse_args(settings)
end

parsed_args = parse_inference_commandline()
infer(parsed_args)