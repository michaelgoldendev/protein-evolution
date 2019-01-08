using Distributions
using LinearAlgebra
using Serialization
using JSON

push!(LOAD_PATH,string(@__DIR__,"/../../../dev/MolecularEvolution/src/"))
using MolecularEvolution

push!(LOAD_PATH,@__DIR__)
using EMNodes
using PDBBuilder

using Binaries
using BranchPaths
using CommonUtils
using LG
using Random
using CTMCs
using Nullables
using FastaIO

secondarystructure = "HBEGITSC"
aminoacids = "ACDEFGHIKLMNPQRSTVWY"

function binarize!(tree::TreeNode)
    nodes = getnodelist(tree)
    counter = 0
    for n in nodes
        while length(n.children) > 2
            c1 = pop!(n.children)
            c2 = pop!(n.children)
            counter +=1 
            push!(n.children, TreeNode(0.0, "binarized_$counter"))
            n.children[end].children = [c1,c2]
            n.children[end].parent = n
        end
    end
end

function pimod(angle::Float64)
  theta = mod2pi(angle)
  if theta > pi
    return theta -2.0*pi
  else
    return theta
  end
end



mutable struct ModelParams
    alphabet::Int
    aminoacidQ::Array{Float64,2}
    numhiddenstates::Int
    initialprobs::Array{Float64,1}
    transitionprobs::Array{Float64,2}
    transitionrates::Array{Float64,2}
	hiddennodes::Array{HiddenNode,1}
	mu::Float64
	matrixcache::Dict{Tuple{Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}

    function ModelParams(aminoacidQ::Array{Float64,2}, numhiddenstates::Int,mu::Float64)
        initialprobs = ones(Float64, numhiddenstates)./numhiddenstates
        transitionprobs = ones(Float64, numhiddenstates, numhiddenstates)./numhiddenstates
        transitionrates = ones(Float64, numhiddenstates, numhiddenstates)
        for h1=1:numhiddenstates
            transitionrates[h1,h1] = 0.0
            for h2=1:numhiddenstates
                if h1 != h2
                    transitionrates[h1,h1] -= transitionrates[h1,h2]
                end
            end
        end
		hiddennodes = HiddenNode[]
		for h=1:numhiddenstates
			push!(hiddennodes, HiddenNode())
		end
		matrixcache = Dict{Tuple{Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}()
        new(numhiddenstates,aminoacidQ,numhiddenstates,initialprobs,transitionprobs,transitionrates,hiddennodes,mu,matrixcache)
    end
end


function absmat(M::Array{Float64,2})
  dim1 = size(M,1)
  dim2 = size(M,2)
  for i=1:dim1
    for j=1:dim2
      if M[i,j] <= 0.0
        M[i,j] = 1e-50
      end
    end
  end
  return M
end

countcachemisses = 0
countcachehits = 0
function reset_matrix_cache(modelparams::ModelParams)
	global countcachemisses
	global countcachehits
	println("RESET: CACHE HITS $(countcachehits) / $(countcachehits+countcachemisses) ($(countcachehits / (countcachehits+countcachemisses)))")
	countcachemisses = 0
	countcachehits = 0
	modelparams.matrixcache = Dict{Tuple{Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}()
end

function getPt(modelparams::ModelParams, prevh::Int, nexth::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = (prevh, nexth)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1

		prevprobs = ones(Float64, modelparams.numhiddenstates)
		if prevh > 0
			prevprobs = modelparams.transitionprobs[prevh,:]
		end
		nextprobs = ones(Float64, modelparams.numhiddenstates)
		if nexth > 0
			nextprobs = modelparams.transitionprobs[:,nexth]
		end
		Q = constructJointMatrix(modelparams, prevprobs, nextprobs)
		
		decomposition = eigen(Q)
		D, V = decomposition.values, decomposition.vectors
		Vi = inv(V)
		modelparams.matrixcache[key] = (Q,V,D,Vi)
		#return Q,exp(Q*t)

		println("CACHE HITS $(countcachehits) / $(countcachehits+countcachemisses) ($(countcachehits / (countcachehits+countcachemisses)))")
	else		
		countcachehits += 1
	end
	
	Q,V,D,Vi = modelparams.matrixcache[key]

	return Q, absmat(real(V*Diagonal(exp.(D*t))*Vi))
end

mutable struct AugmentedNodeData <: NodeData
	branchpath::BranchPath
	dummy::Int
	Rmatrices::Array{Array{Float64,2},1}
	Pmatrices::Array{Array{Float64,2},1}
	vs::Array{Array{Float64,1},1}
	time::Array{Float64,1}
	protein::Protein

	AugmentedNodeData(branchpath::BranchPath, dummy::Int) = new(branchpath, dummy, Array{Float64,2}[], Array{Float64,2}[], Array{Float64,1}[], Float64[], Protein())
end

function gettransprobs(node::TreeNode, selcolin::Int, cols::Array{Int,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	branchiterator = BranchPathIterator(node.data.branchpath,cols)
	P = Matrix{Float64}(I, modelparams.alphabet, modelparams.alphabet)
	Pmatrices = Array{Float64,2}[]
	dummypath = Int[]
	dummytime = Float64[]
	for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
		#=
		prevprobs = ones(Float64,modelparams.numhiddenstates)
		if changecol == 2
			prevh = currstates[1]
			prevprobs = modelparams.transitionprobs[prevh,:]
		end
		succprobs = ones(Float64,modelparams.numhiddenstates)
		if length(currstates) == 3
			succh = currstates[3]
			succprobs = modelparams.transitionprobs[:,succh]
		end
		R = constructJointMatrix(modelparams, prevprobs, succprobs)
		dt = (currtime-prevtime)*node.branchlength
		Pi = exp(R*dt)=#

		
		dt = (currtime-prevtime)*node.branchlength
		prevh = 0
		succh = 0
		if changecol == 2
			prevh = currstates[1]
		end
		if length(currstates) == 3
			succh = currstates[3]
		end
		R,Pi = getPt(modelparams, prevh, succh, dt)

		P *= Pi
		push!(dummypath,0)
		push!(dummytime,prevtime)
	end
	return P
end

function felsensteinhelper(node::TreeNode, selcolin::Int, cols::Array{Int,1}, v::Array{Float64,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	branchiterator = BranchPathIterator(node.data.branchpath,cols)
	P = Matrix{Float64}(I, modelparams.alphabet, modelparams.alphabet)
	Rmatrices = Array{Float64,2}[]
	Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
		#=
		prevprobs = ones(Float64,modelparams.numhiddenstates)
		if changecol == 2
			prevh = currstates[1]
			prevprobs = modelparams.transitionprobs[prevh,:]
		end
		succprobs = ones(Float64,modelparams.numhiddenstates)
		if length(currstates) == 3
			succh = currstates[3]
			succprobs = modelparams.transitionprobs[:,succh]
		end
		R = constructJointMatrix(modelparams, prevprobs, succprobs)
		dt = (currtime-prevtime)*node.branchlength
    	Pi = exp(R*dt)=#
    	
    	dt = (currtime-prevtime)*node.branchlength
    	prevh = 0
		succh = 0
		if changecol == 2
			prevh = currstates[1]
		end
		if length(currstates) == 3
			succh = currstates[3]
		end
		R,Pi = getPt(modelparams, prevh, succh, dt)
    	push!(Rmatrices, R*dt)
    	push!(Pmatrices,Pi)
    	push!(dummytime,prevtime)
    end
    push!(dummytime,1.0)

    tempv = copy(v)
    pushfirst!(vs,v)
    for P in reverse(Pmatrices)
        tempv = P*tempv
    	pushfirst!(vs,copy(tempv))
    end
    popfirst!(vs)

    node.data.Rmatrices = Rmatrices
    node.data.Pmatrices = Pmatrices
    node.data.vs = vs
    node.data.time = dummytime
    return Pmatrices,vs
end


function felsensteinresample(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, selcolin::Int, cols::Array{Int,1}, modelparams::ModelParams)
	#selcol = findfirst(x -> x == selcolin, cols)
	selcol = selcolin
	likelihoods = ones(Float64, length(nodelist), modelparams.alphabet)*-Inf
	logm = zeros(Float64,length(nodelist))

	stack = Int[1]
	while length(stack) > 0
		nodeindex = stack[end]
		node = nodelist[nodeindex]
		if isleafnode(node)
			v = observationlikelihood(node.data.protein, selcolin, modelparams)
			for a=1:length(v)
				likelihoods[nodeindex,a] = v[a]
			end
			pop!(stack)
		else
			leftchildindex = node.children[1].nodeindex
			rightchildindex = node.children[2].nodeindex
			cont = true
			if likelihoods[leftchildindex,1] == -Inf
				push!(stack, leftchildindex)
				cont = false
			end
			if likelihoods[rightchildindex,1] == -Inf
				push!(stack, rightchildindex)
				cont = false
			end

			if cont
                lefttransprobs = gettransprobs(nodelist[leftchildindex], selcol, cols, modelparams)
				righttransprobs = gettransprobs(nodelist[rightchildindex], selcol, cols, modelparams)

        		#likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])
        		likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:]).*observationlikelihood(node.data.protein, selcolin, modelparams)

				Pmatrices_left,vs_left = felsensteinhelper(nodelist[leftchildindex], selcol, cols, likelihoods[leftchildindex,:], modelparams)
				Pmatrices_right,vs_right = felsensteinhelper(nodelist[rightchildindex], selcol, cols, likelihoods[rightchildindex,:], modelparams)

				m = maximum(likelihoods[nodeindex,:])
				likelihoods[nodeindex,:] = likelihoods[nodeindex,:] ./ m
				logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
        		pop!(stack)
        	end
        end
    end

	rootnode = nodelist[1]
	len = length(rootnode.data.branchpath.paths)
	rootseq = Int[rootnode.data.branchpath.paths[i][1] for i=1:len]
	logfreqs = zeros(Float64,modelparams.alphabet)
	prevprobs = ones(Float64, modelparams.numhiddenstates)
	nextprobs = ones(Float64, modelparams.numhiddenstates)
	if selcol > 1
		prevh = rootseq[selcol-1]
		prevprobs = modelparams.transitionprobs[prevh,:]
	end
	if selcol < len
		nexth = rootseq[selcol+1]
		nextprobs = modelparams.transitionprobs[:,nexth]
	end
	for a=1:modelparams.alphabet
		logfreqs[a] = log(prevprobs[a]*nextprobs[a]) # TODO: include initial probs
	end
	freqs = exp.(logfreqs.-maximum(logfreqs))
	freqs /= sum(freqs)
	rootliks = freqs.*likelihoods[1,:]
	rootstate = CommonUtils.sample(rng,rootliks)
	rootnode.data.branchpath.paths[selcol] = Int[rootstate]
	rootnode.data.branchpath.times[selcol] = Float64[0.0]
	print = false
	backwardsampling(rng,nodelist[1], rootstate, selcol,likelihoods,print,modelparams)
end

function backwardsampling(rng::AbstractRNG,node::TreeNode, state::Int, selcol::Int,likelihoods,print::Bool,modelparams::ModelParams)
	for child in node
		path = Int[state]
		for (Pi,v) in zip(child.data.Pmatrices, child.data.vs)
			liks = Pi[path[end],:].*v
			samplestate = CommonUtils.sample(rng,liks)
			push!(path,samplestate)
		end

		newpath = Int[]
		newtime = Float64[]
		for z=1:length(path)-1
			dt = child.data.time[z+1]-child.data.time[z]
			samplepath, sampletimes = modifiedrejectionsampling(rng, child.data.Rmatrices[z], path[z], path[z+1],(modelparams))
			append!(newpath,samplepath)
			append!(newtime,(sampletimes*dt) .+ child.data.time[z])
		end

		newpath, newtime = removevirtualjumps(newpath, newtime)

		if print && isleafnode(child) && child.data.branchpath.paths[selcol][end] != newpath[end]
			println(isleafnode(child),"\t",selcol)
			println(child.data.Pmatrices)
			println(child.data.vs)
			println("J", newpath,newtime)
			println("K", child.data.branchpath.paths[selcol], child.data.branchpath.times[selcol])
		end
		child.data.branchpath.paths[selcol] = newpath
		child.data.branchpath.times[selcol] = newtime
		backwardsampling(rng,child, path[end],selcol, likelihoods,print,modelparams)
	end
end


function augmentedloglikelihood(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, modelparams::ModelParams)
	loglikelihood = 0.0
	for node in nodelist
		if !isroot(node)
			branchiterator = BranchPathIterator(node.data.branchpath,cols)

			for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
				Qii = 0.0
				len = length(cols)
				for selcol=1:len
					prevstate = prevstates[selcol]
					prevprobs = ones(Float64, modelparams.numhiddenstates)
					if selcol > 1
						prevh = prevstates[selcol-1]
						prevprobs = modelparams.transitionprobs[prevh,:]
					end
					nextprobs = ones(Float64, modelparams.numhiddenstates)
					if selcol < len
						nexth = prevstates[selcol+1]
						nextprobs = modelparams.transitionprobs[:,nexth]
					end
					#Qii += getratematrixrow(prevstates, selcol, params, modelspecification,componentindex)[prevstate]
					Qii += constructJointMatrix(modelparams, prevprobs, nextprobs)[prevstate,prevstate]
				end
				dt = (currtime-prevtime)*node.branchlength
				if changecol > 0
					prevstate = prevstates[changecol]
					currstate = currstates[changecol]

					prevprobs = ones(Float64, modelparams.numhiddenstates)
					if changecol > 1
						prevh = currstates[changecol-1]
						prevprobs = modelparams.transitionprobs[prevh,:]
					end
					nextprobs = ones(Float64, modelparams.numhiddenstates)
					if changecol < len
						nexth = currstates[changecol+1]
						nextprobs = modelparams.transitionprobs[:,nexth]
					end
					#Qhi = getratematrixrow(prevstates, changecol, params, modelspecification,componentindex)[currstate]
					Qhi = constructJointMatrix(modelparams, prevprobs, nextprobs)[prevstate,currstate]
					loglikelihood += log(Qhi)
				end
				loglikelihood +=  Qii*dt
			end
		end
	end
	return loglikelihood
end

function constructJointMatrix(modelparams::ModelParams, prevprobs::Array{Float64,1}, succprobs::Array{Float64,1})
	rate = modelparams.mu
	Q = zeros(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)
	for h1=1:modelparams.numhiddenstates
        for h2=1:modelparams.numhiddenstates
			if h1 != h2
				Q[h1,h2] = prevprobs[h2]*succprobs[h2]*rate*modelparams.transitionrates[h1,h2]
				Q[h1,h1] -= Q[h1,h2]
			end
		end
	end
    return Q
end

function siteloglikelihood(site::SiteObservation, h::Int, modelparams::ModelParams)
	ll = 0.0
	if site.aa > 0
		ll += log(modelparams.hiddennodes[h].aa_node.probs[site.aa])
	end

	if site.phi > -100.0
		ll += logpdf(modelparams.hiddennodes[h].phi_node.dist, site.phi)
	end
	if site.omega > -100.0
		ll += logpdf(modelparams.hiddennodes[h].omega_node.dist, site.omega)
	end
	if site.psi > -100.0
		ll += logpdf(modelparams.hiddennodes[h].psi_node.dist, site.psi)
	end
	#=
	if site.bond_angle1 > -100.0
		ll += logpdf(modelparams.hiddennodes[h].bond_angle1_node.dist, site.bond_angle1)
	end
	if site.bond_angle2 > -100.0
		ll += logpdf(modelparams.hiddennodes[h].bond_angle2_node.dist, site.bond_angle2)
	end
	if site.bond_angle3 > -100.0
		ll += logpdf(modelparams.hiddennodes[h].bond_angle3_node.dist, site.bond_angle3)
	end
	if site.bond_length1 > -100.0 && site.bond_length2 > -100.0 && site.bond_length3 > -100.0
		ll += logpdf(modelparams.hiddennodes[h].bond_lengths_node.mvn, Float64[site.bond_length1, site.bond_length2, site.bond_length3])
	end=#
	return ll
end

function observationlikelihood(protein::Protein, col::Int, modelparams::ModelParams)
	if 1 <= col <= length(protein.sites)
	    v = zeros(Float64, modelparams.alphabet)
		for h=1:modelparams.numhiddenstates
			v[h] = siteloglikelihood(protein.sites[col], h, modelparams)
		end
		return exp.(v .- maximum(v))
	else
		return ones(Float64, modelparams.numhiddenstates)
	end
end


function observationloglikelihood(proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, modelparams::ModelParams)
	ll = 0.0
	for node in nodelist
		for col=1:min(length(node.data.branchpath.paths), length(node.data.protein))
			h = node.data.branchpath.paths[col][end]				 
			ll += siteloglikelihood(proteins[node.seqindex].sites[col], h, modelparams)
		end
	end
	return ll
end

function sampletreenode(rng::AbstractRNG, node::TreeNode, modelparams::ModelParams, aligned_sequence::String)
	numcols = length(node.data.branchpath.paths)
	sampled = Protein()
	sampled.name = node.name
	for col=1:numcols
		aa = aligned_sequence[col]
		if aa != '-'
			site = SiteObservation()
			site.aa = indexof(string(aa), aminoacids)
			site.h =  node.data.branchpath.paths[col][end]
			h = site.h
			hiddennode = modelparams.hiddennodes[h]
			site.phi = pimod(vonmisesrand(rng, hiddennode.phi_node.dist))
			site.omega = pimod(vonmisesrand(rng, hiddennode.omega_node.dist))
			site.psi = pimod(vonmisesrand(rng, hiddennode.psi_node.dist))
			bond_lengths = rand(rng, hiddennode.bond_lengths_node.mvn)
			site.bond_length1 = bond_lengths[1]
			site.bond_length2 = bond_lengths[2]
			site.bond_length3 = bond_lengths[3]
			site.bond_angle1 = pimod(vonmisesrand(rng, hiddennode.bond_angle1_node.dist))
			site.bond_angle2 = pimod(vonmisesrand(rng, hiddennode.bond_angle2_node.dist))
			site.bond_angle3 = pimod(vonmisesrand(rng, hiddennode.bond_angle3_node.dist))
			push!(sampled.sites, site)
		end
	end
	return sampled
end

function protein_to_lists(protein::Protein)
	sequence = ""
    phi_psi = Tuple{Float64,Float64}[] 
    omega = Float64[]
    bond_angles = Tuple{Float64,Float64,Float64}[]
    bond_lengths = Tuple{Float64,Float64,Float64}[]
    numcols = length(protein)
    for col=1:numcols
    	sequence = string(sequence, get(aminoacids,protein.sites[col].aa,'-'))
    	push!(phi_psi, (protein.sites[col].phi, protein.sites[col].psi))
    	push!(omega, protein.sites[col].omega)
    	push!(bond_angles, (protein.sites[col].bond_angle1, protein.sites[col].bond_angle2, protein.sites[col].bond_angle3))
    	push!(bond_lengths, (protein.sites[col].bond_length1, protein.sites[col].bond_length2, protein.sites[col].bond_length3))
    end
	return sequence, phi_psi, omega, bond_angles, bond_lengths
end

function compare_branch_scalings(nodelist1,nodelist2)
	branch_scaling = Float64[]
	outputlist = deepcopy(nodelist1)
	branch_scalings = Float64[0.0]
	totalbranchlength1 = 0.0
	totalbranchlength2 = 0.0
	index = 1
	for (n1, n2, o1) in zip(nodelist1, nodelist2, outputlist)
		n1.nodeindex = index
		n2.nodeindex = index
		o1.nodeindex = index
		if !isroot(n1)
			totalbranchlength1 += n1.branchlength
			totalbranchlength2 += n2.branchlength
		end
		index += 1
	end
	for (n1, n2, o1) in zip(nodelist1, nodelist2, outputlist)
		if !isroot(n1)
			b1 = n1.branchlength / totalbranchlength1
			b2 = n2.branchlength / totalbranchlength2
			push!(branch_scalings, b2/b1)
		end
	end

	for o1 in outputlist
		o1.name = string(o1.name, "_$(branch_scalings[o1.nodeindex])")
	end
	return outputlist[1]
end

function training_example_from_sequence_alignment(rng::AbstractRNG, modelparams::ModelParams, fastafile::String)
	println(fastafile)
	sequences = AbstractString[]
    names = AbstractString[]

	FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            push!(sequences, seq)
        end
    end

    newickstring, cachefile = Binaries.fasttreeaa(fastafile)
    root = gettreefromnewick(newickstring)
	binarize!(root)
	nodelist = getnodelist(root)
	for (index,node) in enumerate(nodelist)
		node.nodeindex = index
	end

	numcols = length(sequences[1])
	V = constructJointMatrix(modelparams, ones(Float64,modelparams.numhiddenstates), ones(Float64,modelparams.numhiddenstates))
	paths = Array{Int,1}[]
	times = Array{Float64,1}[]
	for col=1:numcols
		state = rand(rng,1:modelparams.numhiddenstates)
		push!(paths,Int[state,state])
		push!(times,Float64[0.0, 1.0])
	end	
	root.data = AugmentedNodeData(BranchPath(paths,times), 1)
	for node in nodelist
		if !isroot(node)
			paths = Array{Int,1}[]
			times = Array{Float64,1}[]
			parentnode = get(node.parent)
			for col=1:numcols
				parentstate = parentnode.data.branchpath.paths[col][end]
				nodestate = rand(rng,1:modelparams.numhiddenstates)
				path,time = modifiedrejectionsampling(rng, V*0.5, parentstate, nodestate, nothing)
				push!(paths,path)
				push!(times,time)
			end
			node.data = AugmentedNodeData(BranchPath(paths,times), 1)
		end
	end
    
    proteins = Protein[]
	name_protein_dict = Dict{String,Tuple{Int,Protein}}()
	for (p,(sequence,name)) in enumerate(zip(sequences,names))
	    protein = Protein(name)
	    for aa in sequence
		    site = SiteObservation()
		    site.aa = CommonUtils.indexof(string(aa),aminoacids)
		    push!(protein.sites,site)
	    end
	    name_protein_dict[name] = (p, protein)
	    push!(proteins, protein)
	end
	for node in nodelist
		if haskey(name_protein_dict, node.name)
			proteinindex, protein = name_protein_dict[node.name]
			node.seqindex = proteinindex
			node.data.protein = protein
		end
	end

	return (proteins, nodelist)
end

function training_example_from_json_family(rng::AbstractRNG, modelparams::ModelParams, json_family)
	root = gettreefromnewick(json_family["newick_tree"])
	binarize!(root)
	nodelist = getnodelist(root)
	for (index,node) in enumerate(nodelist)
		node.nodeindex = index
	end

	numcols = length(json_family["proteins"][1]["aligned_sequence"])

	V = constructJointMatrix(modelparams, ones(Float64,modelparams.numhiddenstates), ones(Float64,modelparams.numhiddenstates))
	paths = Array{Int,1}[]
	times = Array{Float64,1}[]
	for col=1:numcols
		state = rand(rng,1:modelparams.numhiddenstates)
		push!(paths,Int[state,state])
		push!(times,Float64[0.0, 1.0])
	end	
	root.data = AugmentedNodeData(BranchPath(paths,times), 1)
	for node in nodelist
		if !isroot(node)
			paths = Array{Int,1}[]
			times = Array{Float64,1}[]
			parentnode = get(node.parent)
			for col=1:numcols
				parentstate = parentnode.data.branchpath.paths[col][end]
				nodestate = rand(rng,1:modelparams.numhiddenstates)
				path,time = modifiedrejectionsampling(rng, V*0.5, parentstate, nodestate, nothing)
				push!(paths,path)
				push!(times,time)
			end
			node.data = AugmentedNodeData(BranchPath(paths,times), 1)
		end
	end

	proteins = Protein[]
	name_protein_dict = Dict{String,Tuple{Int,Protein}}()
	for p=1:length(json_family["proteins"])
	    protein = Protein()
	    json_protein = json_family["proteins"][p]
	    for (aa,phi_psi,omega,bond_lengths,bond_angles) in zip(json_protein["aligned_sequence"], json_protein["aligned_phi_psi"], json_protein["aligned_omega"], json_protein["aligned_bond_lengths"], json_protein["aligned_bond_angles"])
	    	phi = phi_psi[1]
	    	psi = phi_psi[2]
	        push!(protein.sites, SiteObservation(0,indexof(string(aa), aminoacids),phi,omega,psi,bond_lengths[1],bond_lengths[2],bond_lengths[3],bond_angles[1],bond_angles[2],bond_angles[3]))
	    end
	    name_protein_dict[json_family["proteins"][p]["name"]] = (p, protein)
	    push!(proteins, protein)
	end
	for node in nodelist
		if haskey(name_protein_dict, node.name)
			proteinindex, protein = name_protein_dict[node.name]
			node.seqindex = proteinindex			
			println("node ", proteinindex)
			node.data.protein = protein
		end
	end
	return (proteins, nodelist,json_family)
end

function getexitrate(node::TreeNode, cols::Array{Int,1}, modelparams::ModelParams)
	len = length(cols)
	exitrate = 0.0
	branchiterator = BranchPathIterator(node.data.branchpath,cols)
	for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
		Qii = 0.0

		for selcol=1:len
			prevstate = prevstates[selcol]
			prevprobs = ones(Float64, modelparams.numhiddenstates)
			if selcol > 1
				prevh = prevstates[selcol-1]
				prevprobs = modelparams.transitionprobs[prevh,:]
			end
			nextprobs = ones(Float64, modelparams.numhiddenstates)
			if selcol < len
				nexth = prevstates[selcol+1]
				nextprobs = modelparams.transitionprobs[:,nexth]
			end
			Qii += constructJointMatrix(modelparams, prevprobs, nextprobs)[prevstate,prevstate]
		end
		dt = (currtime-prevtime)
		exitrate +=  Qii*dt
	end
	return exitrate
end

function proposebranchlength(rng::AbstractRNG, node::TreeNode, cols::Array{Int,1}, modelparams::ModelParams)
	totalexitrate = getexitrate(node, cols, modelparams)
	t =  log(1.0 - rand(rng))/totalexitrate
	newll = log(-totalexitrate) + totalexitrate*t
	oldll = log(-totalexitrate) + totalexitrate*node.branchlength
	propratio = oldll-newll
	return t,propratio
end

function train(numhiddenstates::Int=5)
	rng = MersenneTwister(10498012421321)
	Random.seed!(1234)

	family_dir = "../data/families/"
	family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))

	modelparams = ModelParams(LGmatrix,numhiddenstates,0.2)

	trainingexamples = Tuple[]
	for family_file in family_files[1:30]
		full_path = abspath(joinpath(family_dir, family_file))
		json_family = JSON.parse(open(full_path, "r"))
		if 1 <= length(json_family["proteins"]) <= 1e10
			training_example = training_example_from_json_family(rng, modelparams, json_family)		
			println(json_family["newick_tree"])
			if length(training_example[2][1].children) == 1
				root = training_example[2][1].children[1]
				root.parent = Nullable{TreeNode}()
				training_example = (training_example[1], TreeNode[root], training_example[3])
			end
			push!(trainingexamples, training_example)
		end
	end

	logwriter = open(string("trace", modelparams.numhiddenstates,".log"), "w")
	println(logwriter, "iter\ttotalll\tpathll\tobservationll")
	for iter=1:1000
		reset_matrix_cache(modelparams)
		aacounts = ones(Float64, modelparams.numhiddenstates, 20)*0.1
		transitioncounts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
		transitionratecounts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
		transitionratetotals = zeros(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)
		for h=1:modelparams.numhiddenstates
			transitionratecounts[h,h] = 0.0
		end

		for (proteins,nodelist,json_family) in trainingexamples
			numcols = length(proteins[1])

			for node in nodelist
				if !isroot(node)
					t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
					node.branchlength = t
				end
			end
			
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
				felsensteinresample(rng, proteins, nodelist, col, cols, modelparams)
			end

			for col=1:numcols-1
				for node in nodelist
					branchiterator = BranchPathIterator(node.data.branchpath, Int[col,col+1])
					for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
						transitioncounts[prevstates[1],prevstates[2]] += 1.0
					end
				end
			end

			for col=1:numcols
				for node in nodelist
					if !isroot(node)
						branchiterator = BranchPathIterator(node.data.branchpath, Int[col])
						for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
							if changecol > 0
								transitionratecounts[prevstates[changecol],currstates[changecol]] += 1.0
								dt = currtime-prevtime
								transitionratetotals[prevstates[changecol],currstates[changecol]] += dt*node.branchlength
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
							modelparams.hiddennodes[h].aa_node.counts[site.aa] += 1.0
							push!(modelparams.hiddennodes[h].phi_node.data, site.phi)
							push!(modelparams.hiddennodes[h].omega_node.data, site.omega)
							push!(modelparams.hiddennodes[h].psi_node.data, site.psi)
							push!(modelparams.hiddennodes[h].bond_angle1_node.data, site.bond_angle1)
							push!(modelparams.hiddennodes[h].bond_angle2_node.data, site.bond_angle2)
							push!(modelparams.hiddennodes[h].bond_angle3_node.data, site.bond_angle3)
							add_point(modelparams.hiddennodes[h].bond_lengths_node, Float64[site.bond_length1, site.bond_length2, site.bond_length3])
						end
					end
				end
			end			
		end

		observationll = 0.0
		augmentedll = 0.0
		for (proteins,nodelist) in trainingexamples
			numcols = length(proteins[1])
			augmentedll += augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
			observationll += observationloglikelihood(proteins, nodelist, modelparams)
		end
		println(logwriter, iter-1,"\t",augmentedll+observationll,"\t",augmentedll,"\t",observationll)
		flush(logwriter)

		for h=1:modelparams.numhiddenstates
			for h2=1:modelparams.numhiddenstates
				modelparams.transitionprobs[h,h2] = transitioncounts[h,h2]
			end
			modelparams.transitionprobs[h,:] = modelparams.transitionprobs[h,:] ./ sum(transitioncounts[h,:])

			#transitionratecounts[h,:] = transitionratecounts[h,:] ./ sum(transitionratecounts[h,:])
			#=
			transitionratetotals = transitionratetotals ./ length(trainingexamples)
			modelparams.transitionrates[h,h] = 0.0
			for h2=1:modelparams.numhiddenstates
				if h != h2
					modelparams.transitionrates[h,h2] = (transitionratetotals[h,h2]+transitionratetotals[h2,h])/2.0
					modelparams.transitionrates[h2,h] = modelparams.transitionrates[h,h2]
					modelparams.transitionrates[h,h] -= modelparams.transitionrates[h,h2]
				end
			end=#

			estimate_categorical(modelparams.hiddennodes[h].aa_node)
			estimatevonmises(modelparams.hiddennodes[h].phi_node)
			estimatevonmises(modelparams.hiddennodes[h].psi_node)
			estimatevonmises(modelparams.hiddennodes[h].omega_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle1_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle2_node)
			estimatevonmises(modelparams.hiddennodes[h].bond_angle3_node)
			estimate_multivariate_node(modelparams.hiddennodes[h].bond_lengths_node)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].phi_node.mu,"\t",modelparams.hiddennodes[h].phi_node.kappa,"\t",modelparams.hiddennodes[h].phi_node.N)		
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].psi_node.mu,"\t",modelparams.hiddennodes[h].psi_node.kappa,"\t",modelparams.hiddennodes[h].psi_node.N)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].omega_node.mu,"\t",modelparams.hiddennodes[h].omega_node.kappa,"\t",modelparams.hiddennodes[h].omega_node.N)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle1_node.mu,"\t",modelparams.hiddennodes[h].bond_angle1_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle1_node.N)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle2_node.mu,"\t",modelparams.hiddennodes[h].bond_angle2_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle2_node.N)
			println(iter,"\t",h,"\t",modelparams.hiddennodes[h].bond_angle3_node.mu,"\t",modelparams.hiddennodes[h].bond_angle3_node.kappa,"\t",modelparams.hiddennodes[h].bond_angle3_node.N)
			println(iter,"\t",h,"\t", modelparams.hiddennodes[h].bond_lengths_node.mvn.μ)
			for aa=1:20
				println(iter,"\t",h,"\t",aminoacids[aa],"\t",modelparams.hiddennodes[h].aa_node.probs[aa])
			end			
		end
		println("probs ", modelparams.transitionprobs)
		println("rates ", transitionratecounts)
		println("rates ", transitionratetotals)


		proteins,nodelist,json_family = trainingexamples[1]
	    sequence, phi_psi, omega, bond_angles, bond_lengths = protein_to_lists(sampletreenode(rng, nodelist[2], modelparams, json_family["proteins"][1]["aligned_sequence"]))

	    chain = build_structure_from_angles(sequence, phi_psi, omega, bond_angles, bond_lengths; use_input_bond_angles=true, use_input_bond_lengths=true)
	    fout = open("sample.pdb", "w")
		writepdb(fout, chain)
		close(fout)
	    #=
		for col=1:length(phipsisample)
			node = nodelist[2]
			h = node.data.branchpath.paths[col][end]
			println(pimod(phipsisample[col][3]),"\t", proteins[1][col].omega, "\t", pimod(phipsisample[col][1]),"\t", proteins[1][col].phi,"\t",modelparams.hiddennodes[h].phi_node.dist,"\t",pimod(phipsisample[col][2]),"\t", proteins[1][col].psi,"\t",modelparams.hiddennodes[h].psi_node.dist)
		end=#

		fout = open(string("model_h_",modelparams.numhiddenstates, ".model"), "w")
		Serialization.serialize(fout, modelparams)
		close(fout)
	end
end

function count_aminoacid_substitutions(rng::AbstractRNG, modelparams::ModelParams, node::TreeNode)
	aajumps = 0.0
	for path in node.data.branchpath.paths
		prevaa = CommonUtils.sample(rng, modelparams.hiddennodes[path[1]].aa_node.probs)
		for p=2:length(path)
			nextaa = CommonUtils.sample(rng, modelparams.hiddennodes[path[p]].aa_node.probs)
			if prevaa != nextaa
				aajumps += 1.0
			end
			prevaa = nextaa
		end
	end
	return aajumps
end	

function infer()
	rng = MersenneTwister(10498012421321)
	Random.seed!(1234)

	fin = open("model_h_25.model", "r")
	modelparams = Serialization.deserialize(fin)
	close(fin)

	
	family_dir = "../data/families/"
	family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))
	family_file = joinpath(family_dir, family_files[11])
	
	#family_file = "../data/families/alpha-amylase_NC.ali_23.fam"

	full_path = abspath(family_file)
	json_family = JSON.parse(open(full_path, "r"))
	training_example = training_example_from_json_family(rng, modelparams, json_family)		
	println(json_family["newick_tree"])
	if length(training_example[2][1].children) == 1
		root = training_example[2][1].children[1]
		root.parent = Nullable{TreeNode}()
		training_example = (training_example[1], TreeNode[root], training_example[3])
	end

	proteins,nodelist,json_family = training_example

	#=
	proteins,nodelist = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/alignments/hiv-curated-sel.fasta"))=#

	inputnodelist = deepcopy(nodelist)

	inputtreewriter = open("tree.input.nwk", "w")
	println(inputtreewriter, getnewick(nodelist[1]))
	close(inputtreewriter)

	numcols = length(proteins[1])
	mcmcwriter = open("mcmc.log", "w")
	treewriter = open("tree.mcmc.log", "w")
	print(mcmcwriter, "iter\ttotalll\tpathll\tobservationll")
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
				print(mcmcwriter,"\tevents_node$(node.nodeindex)")
			else 
				print(mcmcwriter,"\tevents_$(node.name)")
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
	branchlength_cache = Dict{Int,Array{Float64,1}}()
	for iter=1:10000
		for node in nodelist
			if !isroot(node)
				t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
				node.branchlength = t
			end
		end
		
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
			felsensteinresample(rng, proteins, nodelist, col, cols, modelparams)
		end

		augmentedll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)	
		observationll = observationloglikelihood(proteins, nodelist, modelparams)
		print(mcmcwriter, iter-1,"\t",augmentedll+observationll,"\t",augmentedll,"\t",observationll)
		for node in nodelist
			if !isroot(node)
				print(mcmcwriter,"\t$(node.branchlength)")
			end
		end		
		for node in nodelist
			if !isroot(node)
				print(mcmcwriter,"\t$(sum([length(path)-1 for path in node.data.branchpath.paths]))")
			end
		end

		outputnodelist = deepcopy(nodelist)
		for (node,outputnode) in zip(nodelist, outputnodelist)
			if !isroot(node)				
				subs_per_site = count_aminoacid_substitutions(rng,modelparams,node)/numcols
				outputnode.branchlength = subs_per_site				
				print(mcmcwriter,"\t$(subs_per_site)")
				cached_branchlengths = get(branchlength_cache, outputnode.nodeindex, Float64[])
				push!(cached_branchlengths, outputnode.branchlength)
				branchlength_cache[outputnode.nodeindex] = cached_branchlengths
			end
		end
		println(mcmcwriter)
		flush(mcmcwriter)

		println(treewriter, getnewick(outputnodelist[1]))
		flush(treewriter)

		consensustreewriter = open("tree.mean.consensus.nwk", "w")
		for outputnode in outputnodelist
			if !isroot(outputnode)
				outputnode.branchlength = mean(branchlength_cache[outputnode.nodeindex][max(1,div(length(branchlength_cache[outputnode.nodeindex]),2)):end])
			end
		end
		println(consensustreewriter, getnewick(outputnodelist[1]))
		close(consensustreewriter)

		treewriter = open("tree.scalings.nwk", "w")
		println(treewriter, getnewick(compare_branch_scalings(inputnodelist,outputnodelist)))
		close(treewriter)

		consensustreewriter = open("tree.median.consensus.nwk", "w")
		for outputnode in outputnodelist
			if !isroot(outputnode)
				outputnode.branchlength = median(branchlength_cache[outputnode.nodeindex][max(1,div(length(branchlength_cache[outputnode.nodeindex]),2)):end])
			end
		end
		println(consensustreewriter, getnewick(outputnodelist[1]))
		close(consensustreewriter)
	end
	close(mcmcwriter)
	close(treewriter)
end


using PyPlot
function plot_nodes()
	fin = open("model_h_25.model", "r")
	modelparams = Serialization.deserialize(fin)
	close(fin)

	nrows = 2
	ncols = 3


	for h=1:modelparams.numhiddenstates
		fig = plt[:figure](figsize=(10, 6))
		plt[:clf]
		plt[:rc]("text", usetex=true)
		plt[:rc]("font", family="serif")

		plt[:suptitle]("Hidden state $(h)", fontsize=16)

		N = 500
		mat = zeros(Float64, N, N)
		for (i, x) in enumerate(range(-pi,stop=pi,length=N))
			for (j, y) in enumerate(range(-pi,stop=pi,length=N))
				mat[N-j+1,i] = pdf(modelparams.hiddennodes[h].phi_node.dist, x)
				mat[N-j+1,i] *= pdf(modelparams.hiddennodes[h].psi_node.dist, y)
			end
		end
		ax = plt[:subplot](nrows, ncols, 1, aspect="equal")
		ax[:imshow](mat)
		angle_tick_positions = [0, div(N-1,4), div(N-1,2), div((N-1)*3,4), N-1]
		angle_labels = ["\$-\\pi\$","\$-\\pi/2\$", "0", "\$\\pi/2\$", "\$\\pi\$"]
		ax[:set_xticks](angle_tick_positions)
		ax[:set_yticks](angle_tick_positions)
		ax[:set_xticklabels](angle_labels)
		ax[:set_yticklabels](reverse(angle_labels))
		plt[:xlabel]("Phi (\$\\phi\$)")
		plt[:ylabel]("Psi (\$\\psi\$)")
		
		ax = plt[:subplot](nrows, ncols, 2, aspect="auto")
		x = range(-pi,stop=pi,length=N)
		y = pdf.(modelparams.hiddennodes[h].omega_node.dist, x)
		ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
		ax[:set_xlim](-pi, pi)
		ax[:set_ylim](0.0, 10.0)
		ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
		ax[:set_xticklabels](angle_labels)
		plt[:xlabel]("Omega (\$\\omega\$)")

		aminoacidstext = ["Ala","Cys","Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr"]
	    barcolors = ["#777775", "#fedd00", "#ef3340", "#ef3340", "#000000", "#fedd00", "#0087c7", "#333334", "#0087c7", "#333334", "#333334", "#65428a", "#fedd00", "#65428a", "#0087c7", "#0087c7", "#333334", "#777775", "#000000", "#000000"]    
	    barwidth = 0.5
	    ax = plt[:subplot](nrows, ncols, 3, aspect="auto")
		ax[:bar](1:20 + barwidth, modelparams.hiddennodes[h].aa_node.probs, barwidth, color=barcolors, label="Men")
		ax[:set_ylim](0.0, 1.0)
		ax[:set_xticks](1:20 + barwidth + barwidth/2.0)    
	    ax[:set_xticklabels](aminoacidstext, rotation="vertical", fontsize=7)

	    ax = plt[:subplot](nrows, ncols, 4, aspect="auto")
		x = range(-pi,stop=pi,length=N)
		y = pdf.(modelparams.hiddennodes[h].bond_angle1_node.dist, x)
		ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
		ax[:set_xlim](-pi, pi)
		ax[:set_ylim](0.0, 15.0)
		ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
		ax[:set_xticklabels](angle_labels)
		plt[:xlabel]("\${\\angle}C_{i-1}-N_{i}-CA_{i}\$")

		ax = plt[:subplot](nrows, ncols, 5, aspect="auto")
		x = range(-pi,stop=pi,length=N)
		y = pdf.(modelparams.hiddennodes[h].bond_angle2_node.dist, x)
		ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
		ax[:set_xlim](-pi, pi)
		ax[:set_ylim](0.0, 15.0)
		ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
		ax[:set_xticklabels](angle_labels)
		plt[:xlabel]("\${\\angle}N_{i}-CA_{i}-C_{i}\$")

		ax = plt[:subplot](nrows, ncols, 6, aspect="auto")
		x = range(-pi,stop=pi,length=N)
		y = pdf.(modelparams.hiddennodes[h].bond_angle3_node.dist, x)
		ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
		ax[:set_xlim](-pi, pi)
		ax[:set_ylim](0.0, 15.0)
		ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
		ax[:set_xticklabels](angle_labels)
		plt[:xlabel]("\${\\angle}CA_{i}-C_{i}-N_{i+1}\$")

		plt[:subplots_adjust](hspace=0.4)
		plt[:savefig]("plot$(h).png")
		plt[:close]
	end
end

train()
#plot_nodes()
#infer()
