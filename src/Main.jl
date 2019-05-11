prog = "Protein sequence and structure evolution"
version = "0.1"

using Distributions
using LinearAlgebra
using Serialization
using JSON
using Statistics

push!(LOAD_PATH,string(@__DIR__,"/../../MolecularEvolution/src/"))
using MolecularEvolution

push!(LOAD_PATH,@__DIR__)
using EMNodes
using BivariateVonMises
using Backbone
using AngleUtils

using Binaries
using BranchPaths
using CommonUtils
using LG
using Random
using CTMCs
using Nullables
using FastaIO
using Formatting
using Printf
using NLopt
using SpecialFunctions
using ArgParse



secondarystructure = "HBEGITSC"
aminoacids = "ACDEFGHIKLMNPQRSTVWY"

mutable struct ModelParams
    alphabet::Int
    aminoacidQ::Array{Float64,2}
    aa_exchangeablities::Array{Float64,2}
    rate_alpha::Float64
    numrates::Int
    rates::Array{Float64,1}
    rate_freqs::Array{Float64,1}
    numhiddenstates::Int
    initialprobs::Array{Float64,1}
    transitioncounts::Array{Float64,2}
    transitionprobs::Array{Float64,2}
    transitionrates::Array{Float64,2}
	hiddennodes::Array{HiddenNode,1}
	mu::Float64
	hiddenmu::Float64
	matrixcache::Dict{Tuple{Int,Int,Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}
	branchscalingfactor::Float64
	scalingfactor::Float64
	usestructureobservations::Bool
	hidden_conditional_on_aa::Bool
	ratemode::Int
	aarates::Array{Float64,1}
	use_bivariate_von_mises::Bool

    function ModelParams(aminoacidQ::Array{Float64,2}, numhiddenstates::Int,mu::Float64=1.0,hiddenmu::Float64=1.0)
        initialprobs = ones(Float64, numhiddenstates)./numhiddenstates
        transitioncounts = ones(Float64, numhiddenstates, numhiddenstates)*0.1
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
			push!(hiddennodes, HiddenNode(LGfreqs))
		end
		println("AAFreqs",hiddennodes[1].aa_node.probs)
		matrixcache = Dict{Tuple{Int,Int,Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}()
		numrates = 10
		rate_alpha = 0.5
		#rates = discretizegamma(rate_alpha,1.0/rate_alpha,numrates)
		rates = ones(Float64,numrates)
		rate_freqs = ones(Float64,numrates)/numrates
		#
        new(20,aminoacidQ,LGexchangeability, rate_alpha, numrates, rates,rate_freqs,numhiddenstates,initialprobs,transitioncounts,transitionprobs,transitionrates,hiddennodes,mu,hiddenmu,matrixcache,1.0,1.0,true,true,1, ones(Float64,numhiddenstates), true)
    end
end

function  modifiedrejectionsampling(rng::AbstractRNG, Q::Array{Float64,2}, x0::Int, xt::Int, modelparams::ModelParams, key::Tuple{Int,Int,Int,Int}=(-1000,-1000,-1000,-1000))
	maxiters1 = 1000
	maxiters2 = 2000
	path = Int[]
	times = Float64[]
	S = copy(Q)
	for i=1:size(S,1)
		S[i,i] = 0.0
		S[i,:] ./= sum(S[i,:])
	end

	time1 = time()
	count = 0
	success = true
	while success
		path = Int[x0]
		times = Float64[0.0]
		totalt = 0.0
		if x0 != xt
			r = rand(rng)
			totalt += log(1.0-r*(1.0-exp(Q[path[end],path[end]])))/Q[path[end],path[end]]
			push!(path, CommonUtils.sample(rng, S[path[end],:]))
			push!(times, totalt)
		end
		if count > maxiters1
			success = false
			break
		end
		while true
			if count > maxiters2
				success = false
				break
			end

			r = rand(rng)
			samplet = log(1.0-r)/Q[path[end],path[end]]
			totalt += samplet
			if totalt > 1.0
				break
			end
			push!(path, CommonUtils.sample(rng, S[path[end],:]))
			push!(times, totalt)
			count += 1
		end
		if path[end] == xt
			break
		end

	end
	time2 = time()

	if !success
		if key[1] != -1000 && haskey(modelparams.matrixcache,key)
			dummy,V,D,Vi = modelparams.matrixcache[key]
			#println("cached")
			path,times = CTMCs.recursivesampling(rng, Q, x0, xt, V, D, Vi)
			time3 = time()
			if rand(rng) < 0.001
				println(time2-time1,"\t",time3-time2)	
			end
		else
			#println("not cached")
			path,times = CTMCs.recursivesampling(rng, Q, x0, xt)
		end
		#path, times = approximatesampling(rng, Q, 1.0, x0, xt)
		return path,times,true
	else
		return path, times, success
	end
end

mutable struct Sample 
	name::String
	phipsisamples::Array{Array{Tuple{Float64,Float64},1},1}
	aasamples::Array{Array{Int,1},1}
	hiddensamples::Array{Array{Int,1},1}
	json_family::Dict{String,Any}
	modelparams::ModelParams

	function Sample(name::String, modelparams::ModelParams)
		new(name, Array{Tuple{Float64,Float64},1}[], Int[], Int[], Dict{String,Any}(), modelparams)
	end
end



function binarize!(tree::TreeNode)
    nodes = getnodelist(tree)
    counter = 0
    for n in nodes
    	binarized = false
        while length(n.children) > 2
            c1 = pop!(n.children)
            c2 = pop!(n.children)
            counter +=1 
            push!(n.children, TreeNode(0.0, "binarized_$counter"))
            n.children[end].children = [c1,c2]
            c1.parent = Nullable{TreeNode}(n.children[end])
            c2.parent = Nullable{TreeNode}(n.children[end])
            n.children[end].parent = Nullable{TreeNode}(n)
            binarized = true
        end
        if binarized
        	totalbranchlength = sum(Float64[n.branchlength for n in n.children])
        	avglen = totalbranchlength/length(n.children)
        	for c in n.children
        		c.branchlength = avglen
        	end
        end
    end


end

function nodedistances(nodelist::Array{TreeNode,1})
	for (index,node) in enumerate(nodelist)
		node.nodeindex = index
	end

	distmatrix = ones(Float64, length(nodelist), length(nodelist))*Inf
	for node in nodelist
		distmatrix[node.nodeindex,node.nodeindex] = 0.0
		for child in node
			distmatrix[node.nodeindex,child.nodeindex] = child.branchlength
			distmatrix[child.nodeindex,node.nodeindex] = child.branchlength
		end
	end

	mindistmatrix = copy(distmatrix)
	for n1 in nodelist
		for v1 in nodelist
			for v2 in nodelist
				mindistmatrix[n1.nodeindex,v2.nodeindex] = min(mindistmatrix[n1.nodeindex,v2.nodeindex], mindistmatrix[n1.nodeindex,v1.nodeindex]+mindistmatrix[v1.nodeindex,v2.nodeindex])
			end
		end
	end

	return mindistmatrix
end

function reorient(node, new_parent, new_branch_length)
    newchildren = TreeNode[]
    for c in node
    	if c != new_parent
    		push!(newchildren,c)
    	end
    end
    node.children = newchildren
    # If this node has a parent, reorient the parent.
    if !isnull(node.parent)
        parent = node.parent.value
        reorient(parent, node, node.branchlength)
        # Add the parent as a child
        push!(node.children, parent)
    end
    # Set then new parent as the parent, with the new branch length
    node.parent = Nullable{TreeNode}(new_parent)
    node.branchlength = new_branch_length
end

function reroot(child_node, dist_above_child=(child_node.branchlength/2))
    if dist_above_child > child_node.branchlength
        print("This isn't going to work")
    end
        
    # Remembering stuff
    dist_below_parent = child_node.branchlength - dist_above_child
    old_parent = child_node.parent.value
        
    new_root = TreeNode(0.0,"root")
    child_node.branchlength = dist_above_child
    reorient(old_parent, child_node, dist_below_parent)
    new_root.children = [child_node, old_parent]
    return new_root
end


function midpoint_root(nodelistin::Array{TreeNode,1})
	nodelist = deepcopy(nodelistin)
	for node in nodelist
		if startswith(node.name,":") || node.name == ""
			node.name = string("node",node.nodeindex)
		end
	end

	mindistmatrix = nodedistances(nodelist)
	minnodetotipdistance = Inf
	minnodeindex = 0
	for nonleaf in nodelist
		if !isleafnode(nonleaf)			
			nodetotipdistance = 0.0
			vals = Float64[]
			for leaf in nodelist
				if isleafnode(leaf)
					nodetotipdistance += mindistmatrix[nonleaf.nodeindex, leaf.nodeindex]
					push!(vals, mindistmatrix[nonleaf.nodeindex, leaf.nodeindex])
				end
			end
			stdeviation = std(vals)
			if stdeviation < minnodetotipdistance
				minnodetotipdistance = stdeviation
				minnodeindex = nonleaf.nodeindex
			end
		end
	end

	newroot = reroot(nodelist[minnodeindex], nodelist[minnodeindex].branchlength/2.0)
	nodelist = getnodelist(newroot)
	for node in nodelist
		if length(node.children) == 1
			node.children = node.children[1].children
			for c in node.children
				c.parent = Nullable{TreeNode}(node)
			end
		end
	end

	return newroot
end

function testprint(node::TreeNode, seen::Array{Int,1}=Int[])
	println(node.nodeindex,"\t", [c.nodeindex for c in node])
	if node.nodeindex in seen
		exit()
	end
	push!(seen, node.nodeindex)	

	for c in node
		testprint(c, seen)
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

function discretizegamma(shape::Float64, scale::Float64, numcategories::Int)
  if numcategories == 1
    return Float64[1.0]
  else
    catwidth = 1.0 / numcategories
    vs = Float64[catwidth/2.0 + (i-1.0)*catwidth for i=1:numcategories]
    gammadist = Gamma(shape, scale)
    return Float64[quantile(gammadist, v) for v in vs]
  end
end




function estimate_hidden_transition_probs(modelparams::ModelParams)
	freqs = zeros(Float64,modelparams.numhiddenstates)
	for h=1:modelparams.numhiddenstates
		freqs[h] = sum(modelparams.transitioncounts[:,h])
	end
	freqs /= sum(freqs)
	
	#=
	transitioncounts = copy(modelparams.transitioncounts)
	#=
	for h1=1:modelparams.numhiddenstates
		for h2=1:modelparams.numhiddenstates
			transitioncounts[h1,h2] /= freqs[h2]
		end
	end=#

	for h1=1:modelparams.numhiddenstates
		for h2=1:modelparams.numhiddenstates
			modelparams.transitionprobs[h1,h2] = (transitioncounts[h1,h2]+transitioncounts[h2,h1])
		end
	end
	for h1=1:modelparams.numhiddenstates
		modelparams.transitionprobs[h1,:] = modelparams.transitionprobs[h1,:] / sum(modelparams.transitionprobs[h1,:])
	end=#
	
	for h=1:modelparams.numhiddenstates
		modelparams.transitionprobs[h,:] = modelparams.transitioncounts[h,:] ./ sum(modelparams.transitioncounts[h,:])
	end
	println(modelparams.transitionprobs)
	println(modelparams.transitionprobs^40)
	println(freqs)
	modelparams.transitioncounts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
end

function estimate_hidden_mixture(modelparams::ModelParams)
	freqs = zeros(Float64,modelparams.numhiddenstates)
	for h=1:modelparams.numhiddenstates
		freqs[h] = sum(modelparams.transitioncounts[:,h])
	end
	freqs /= sum(freqs)
	for h=1:modelparams.numhiddenstates
		modelparams.transitionprobs[h,:] = freqs
	end
	println(modelparams.transitionprobs)
	println(modelparams.transitionprobs^40)
	println(freqs)
	modelparams.transitioncounts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
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
function reset_matrix_cache(modelparams::ModelParams, printflag::Bool=true)
	global countcachemisses
	global countcachehits
	if printflag
		println("RESET: CACHE HITS $(countcachehits) / $(countcachehits+countcachemisses) ($(countcachehits / (countcachehits+countcachemisses)))")
	end
	countcachemisses = 0
	countcachehits = 0
	modelparams.matrixcache = Dict{Tuple{Int,Int,Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}()
end

function getAAkey(h::Int)
	return (-1,-1,h,0)
end

function getAAandPt(modelparams::ModelParams, prev_hmm::Int, next_hmm::Int, h::Int, rate::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = getAAkey(h)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1
		Q = constructAAMatrix(modelparams, prev_hmm, next_hmm, h)
		
		
		decomposition = eigen(Q)
		D, V = decomposition.values, decomposition.vectors
		Vi = inv(V)
		modelparams.matrixcache[key] = (Q,V,D,Vi)
		#return Q,exp(Q*t)

		#println("CACHE HITS $(countcachehits) / $(countcachehits+countcachemisses) ($(countcachehits / (countcachehits+countcachemisses)))")
	else		
		countcachehits += 1
	end
	
	Q,V,D,Vi = modelparams.matrixcache[key]

	sitescale = modelparams.rates[rate]
	Q *= sitescale

	if t < 0.0
		return key, Q, Q
	else
		return key, Q, absmat(real(V*Diagonal(exp.(D*t*sitescale))*Vi))
	end
end

function getQkey(prevh::Int, nexth::Int, aa::Int)
	return (-1,prevh, nexth, aa)
end

function getQandPt(modelparams::ModelParams, prevh::Int, nexth::Int, aa::Int, rate::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = getQkey(prevh,nexth,aa)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1

		Q = constructHiddenMatrix(modelparams, prevh, nexth, aa)
		
		decomposition = eigen(Q)
		D, V = decomposition.values, decomposition.vectors
		Vi = inv(V)
		modelparams.matrixcache[key] = (Q,V,D,Vi)

		#println("CACHE HITS $(countcachehits) / $(countcachehits+countcachemisses) ($(countcachehits / (countcachehits+countcachemisses)))")
	else		
		countcachehits += 1
	end
	
	
	Q,V,D,Vi = modelparams.matrixcache[key]

	sitescale = modelparams.rates[rate]
	Q *= sitescale

	if t < 0.0
		return key, Q, Q
	else
		return key, Q, absmat(real(V*Diagonal(exp.(D*t*sitescale))*Vi))
	end
end

function getJointPt(modelparams::ModelParams, prev_hmm::Int, next_hmm::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = (-2,-2, prev_hmm,next_hmm)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1
		Q = constructJointMatrix(modelparams, prev_hmm, next_hmm)
		
		decomposition = eigen(Q)
		D, V = decomposition.values, decomposition.vectors
		Vi = inv(V)
		modelparams.matrixcache[key] = (Q,V,D,Vi)
		#return Q,exp(Q*t)

		#println("CACHE HITS $(countcachehits) / $(countcachehits+countcachemisses) ($(countcachehits / (countcachehits+countcachemisses)))")
	else		
		countcachehits += 1
	end
	
	Q,V,D,Vi = modelparams.matrixcache[key]

	if t < 0.0
		return Q, Q
	else
		return Q, absmat(real(V*Diagonal(exp.(D*t))*Vi))
	end
end

mutable struct AugmentedNodeData <: NodeData
	branchpath::BranchPath
	aabranchpath::BranchPath
	ratesbranchpath::BranchPath
	jointbranchpath::BranchPath
	dummy::Int
	protein::Protein
	inputbranchlength::Float64
	fixbranchlength::Bool
	proposallikelihood::Array{Float64,1}

	AugmentedNodeData(col::Int) = new(BranchPath(col),BranchPath(col),BranchPath(col),BranchPath(col),1, Protein(),1.0, false, zeros(Float64,col)) 
	AugmentedNodeData(branchpath::BranchPath, aabranchpath::BranchPath, ratesbranchpath::BranchPath, jointbranchpath::BranchPath, dummy::Int) = new(branchpath, aabranchpath, ratesbranchpath, jointbranchpath, dummy, Protein(),1.0,false, zeros(Float64,1))
end

function felsensteinhelper_joint(node::TreeNode, selcolin::Int, cols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol])])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	Pret = Matrix{Float64}(I, modelparams.alphabet*modelparams.numhiddenstates, modelparams.alphabet*modelparams.numhiddenstates)
	Rmatrices = Array{Float64,2}[]
	Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	
    for (index, it) in enumerate(multi_iter)
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		R,Pi = getJointPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol+1, 0), dt)
		Pret *= Pi
		push!(Rmatrices, R*dt)
    	push!(Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)

    	if index == 1
    		node.data.jointbranchpath.R, node.data.jointbranchpath.P2 =  getJointPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol+1, 0), node.branchlength)
    	end
	end
	push!(dummytime,1.0)

    tempv = copy(v)
    pushfirst!(vs,v)
    for P in reverse(Pmatrices)
        tempv = P*tempv
    	pushfirst!(vs,copy(tempv))
    end
    popfirst!(vs)

    node.data.jointbranchpath.Rmatrices = Rmatrices
    node.data.jointbranchpath.Pmatrices = Pmatrices
    node.data.jointbranchpath.vs = vs
    node.data.jointbranchpath.time = dummytime
    node.data.jointbranchpath.P = Pret
    return Pret,node.data.branchpath.Pmatrices,vs
end


function felsensteinresample_joint(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, selcolin::Int, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams, sample::Bool=true)
	#selcol = findfirst(x -> x == selcolin, cols)
	selcol = selcolin
	likelihoods = ones(Float64, length(nodelist), modelparams.numhiddenstates*modelparams.alphabet)*-Inf
	logm = zeros(Float64,length(nodelist))

	stack = Int[1]
	while length(stack) > 0
		nodeindex = stack[end]
		node = nodelist[nodeindex]
		if isleafnode(node)
			v = observationlikelihood(node.data.protein, selcolin, modelparams)
			if 1 <= selcolin <= length(node.data.protein.sites) && node.data.protein.sites[selcolin].aa > 0
				numstates = modelparams.numhiddenstates*modelparams.alphabet
				for a=1:numstates
					likelihoods[nodeindex,a] = 0.0
				end
				for h=1:modelparams.numhiddenstates
					aa = node.data.protein.sites[selcolin].aa
					index = (h-1)*modelparams.alphabet + aa
					likelihoods[nodeindex,index] = v[h]
				end
			else
				for h=1:modelparams.numhiddenstates
					for aa=1:modelparams.alphabet
						index = (h-1)*modelparams.alphabet + aa
						likelihoods[nodeindex,index] = v[h]
					end
				end
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
				lefttransprobs, Pmatrices_left, vs_left = felsensteinhelper_joint(nodelist[leftchildindex], selcol, cols, aacol, likelihoods[leftchildindex,:], modelparams)
				righttransprobs, Pmatrices_right, vs_right = felsensteinhelper_joint(nodelist[rightchildindex], selcol, cols, aacol, likelihoods[rightchildindex,:], modelparams)

				likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])

				m = maximum(likelihoods[nodeindex,:])
				likelihoods[nodeindex,:] = likelihoods[nodeindex,:] ./ m
				logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
        		pop!(stack)
        	end
        end
    end

	rootnode = nodelist[1]

	len = length(rootnode.data.jointbranchpath.paths)
	prevh = 0
	prevprobs = ones(Float64, modelparams.numhiddenstates)
	if selcol > 1
		prevh = rootnode.data.branchpath.paths[selcol-1][1]
		prevprobs = modelparams.transitionprobs[prevh,:]
	end
	nexth = 0
	nextprobs = ones(Float64, modelparams.numhiddenstates)
	if selcol < len
		nexth = rootnode.data.branchpath.paths[selcol+1][1]
		nextprobs = modelparams.transitionprobs[:,nexth]
	end
	v = Float64[]
	for h=1:modelparams.numhiddenstates
		for aa=1:modelparams.alphabet
			push!(v, prevprobs[h]*nextprobs[h]*modelparams.hiddennodes[h].aa_node.probs[aa])
		end
	end
	rootliks = v.*likelihoods[1,:]
	if sample
		rootstate = CommonUtils.sample(rng,rootliks)
		rootnode.data.jointbranchpath.paths[selcol] = Int[rootstate]
		rootnode.data.jointbranchpath.times[selcol] = Float64[0.0]
		print = false
		cont = backwardsampling_joint(rng, nodelist[1], rootstate, selcol, modelparams)
		#cont = backwardsampling_joint2(rng, nodelist[1], rootstate, selcol, modelparams)

		for node in nodelist
			joint_to_aa_h!(modelparams, node, selcol)
		end
		return cont
	else
		return log(sum(rootliks)) + logm[1]
	end
end

function joint_to_aa_h!(modelparams::ModelParams, node::TreeNode, selcol::Int)
	newpath = node.data.jointbranchpath.paths[selcol]
	newtime = node.data.jointbranchpath.times[selcol]

	node.data.aabranchpath.paths[selcol] = Int[(a-1)%modelparams.alphabet+1 for a in newpath]
	node.data.aabranchpath.times[selcol] = copy(newtime)
	node.data.branchpath.paths[selcol] = Int[div(a-1,modelparams.alphabet)+1 for a in newpath]
	node.data.branchpath.times[selcol] = copy(newtime)

	node.data.aabranchpath.paths[selcol], node.data.aabranchpath.times[selcol] = removevirtualjumps(node.data.aabranchpath.paths[selcol], node.data.aabranchpath.times[selcol])
	node.data.branchpath.paths[selcol], node.data.branchpath.times[selcol] = removevirtualjumps(node.data.branchpath.paths[selcol], node.data.branchpath.times[selcol])
end

function aa_h_to_joint!(modelparams::ModelParams, node::TreeNode, selcol::Int)
	i = 1
	j = 1

	aapath = node.data.aabranchpath.paths[selcol]
	aatime = node.data.aabranchpath.times[selcol] 
	hpath = node.data.branchpath.paths[selcol]
	htime = node.data.branchpath.times[selcol]

	mintime = 0.0
	newpath = Int[(hpath[j]-1)*20 + aapath[i]]
	newtime = Float64[0.0]
	while true
		if i+1 <= length(aapath) && j+1 <= length(hpath)
			if aatime[i+1] < htime[j+1]
				i += 1
				push!(newpath, (hpath[j]-1)*20 + aapath[i])
				push!(newtime, aatime[i])
			else
				j += 1
				push!(newpath, (hpath[j]-1)*20 + aapath[i])
				push!(newtime, htime[j])
			end
		elseif i+1 <= length(aapath)
			i += 1
			push!(newpath, (hpath[j]-1)*20 + aapath[i])
			push!(newtime, aatime[i])
		elseif j+1 <= length(hpath)
			j += 1
			push!(newpath, (hpath[j]-1)*20 + aapath[i])
			push!(newtime, htime[j])
		else 
			break
		end
	end

	node.data.jointbranchpath.paths[selcol] = copy(newpath)
	node.data.jointbranchpath.times[selcol] = copy(newtime)
end

function backwardsampling_joint(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		proposallikelihood = 0.0
		P = child.data.jointbranchpath.P
		R = child.data.jointbranchpath.R*child.branchlength
		liks = P[state,:].*child.data.jointbranchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime, success = modifiedrejectionsampling(rng, R, state, b, modelparams)
		if !success
			return false
		end
		child.data.jointbranchpath.paths[aacol] = newpath
		child.data.jointbranchpath.times[aacol] = newtime
		for z=1:length(newpath)-1
			a = newpath[z]
			b = newpath[z+1]
			dt = newtime[z+1]-newtime[z]
			proposallikelihood += log(R[a,b])
			proposallikelihood += R[a,a]*dt
		end
		a = newpath[end]
		proposallikelihood += R[a,a]*(1.0 - newtime[end])
		child.data.proposallikelihood[aacol] = proposallikelihood

		success = backwardsampling_joint(rng,child, b,aacol, modelparams)
		if !success
			return false
		end
	end
	return true
end

function backwardsampling_joint2(rng::AbstractRNG,node::TreeNode, state::Int, selcol::Int, modelparams::ModelParams)
	if isroot(node)
		aacol = selcol
		rootnode = node
		len = length(rootnode.data.branchpath.paths)
		prevh = 0
		if aacol > 1
			prevh = rootnode.data.branchpath.paths[aacol-1][end]
		end
		nexth = 0
		if aacol < len
			nexth = rootnode.data.branchpath.paths[aacol+1][end]
		end
		jstate = rootnode.data.jointbranchpath.paths[aacol][end]
		aa =  ((jstate-1) % modelparams.alphabet)+1
		freqs = gethiddeninitialprobs(modelparams, prevh, nexth, aa)
		h = div(jstate-1, modelparams.alphabet)+1 
		node.data.proposallikelihood[selcol] = log(freqs[h])
	end

	for child in node
		proposallikelihood = 0.0
		path = Int[state]
		for (Pi,v) in zip(child.data.jointbranchpath.Pmatrices, child.data.jointbranchpath.vs)
			liks = Pi[path[end],:].*v
			samplestate = CommonUtils.sample(rng,liks)
			push!(path,samplestate)
		end

		newpath = Int[]
		newtime = Float64[]
		for z=1:length(path)-1
			dt = child.data.jointbranchpath.time[z+1]-child.data.jointbranchpath.time[z]
			R = child.data.jointbranchpath.Rmatrices[z]
			samplepath, sampletimes, success = modifiedrejectionsampling(rng, R, path[z], path[z+1], modelparams)
			if !success
				return false
			end
			append!(newpath,samplepath)
			append!(newtime,(sampletimes*dt) .+ child.data.jointbranchpath.time[z])

			for z=1:length(samplepath)-1
				dt2 = sampletimes[z+1]-sampletimes[z]
				a = samplepath[z]
				b = samplepath[z+1]
				proposallikelihood += log(R[a,b])
				proposallikelihood += R[a,a]*dt2
			end
			a = samplepath[end]
			proposallikelihood += R[a,a]*(1.0-sampletimes[end])
		end
		child.data.proposallikelihood[selcol] = proposallikelihood

		

		newpath, newtime = removevirtualjumps(newpath, newtime)

		child.data.jointbranchpath.paths[selcol] = newpath
		child.data.jointbranchpath.times[selcol] = newtime

		joint_to_aa_h!(modelparams, child, selcol)

		success = backwardsampling_joint2(rng,child, path[end],selcol,modelparams)
		if !success
			return false
		end
	end
	return true
end

function felsensteinhelper_likelihood(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
	selcol = findfirst(x -> x == aacol, cols)
	loglikelihood = log(modelparams.hiddennodes[nodelist[1].data.branchpath.paths[aacol][end]].aa_node.probs[nodelist[1].data.aabranchpath.paths[aacol][end]])
	for node in nodelist
		if !isroot(node)
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols),BranchPathIterator(node.data.aabranchpath,Int[aacol])])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]

		    for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength	

				loglikelihood += gethiddenentry(modelparams, get(hiddeniter.prevstates,selcol-1,0), get(hiddeniter.prevstates,selcol+1,0), hiddeniter.prevstates[selcol],  hiddeniter.prevstates[selcol], aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end])*dt
				changecol = 0
				if multi_iter.branchpathindex == 1
					changecol = hiddeniter.mincol
				else
					changecol = aaiter.mincol
				end
				if multi_iter.branchpathindex == 1 && changecol == selcol
					loglikelihood += log(gethiddenentry(modelparams, get(hiddeniter.prevstates,selcol-1,0), get(hiddeniter.prevstates,selcol+1,0), hiddeniter.prevstates[selcol], hiddeniter.currstates[selcol], aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end])*node.branchlength)
				end
			end
		end
	end
    return loglikelihood
end

function felsensteinhelper_hidden(node::TreeNode, selcolin::Int, incols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, incols)
	cols = copy(incols)
	deleteat!(cols,selcol)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol])])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	Pret = Matrix{Float64}(I, modelparams.numhiddenstates, modelparams.numhiddenstates)
	node.data.branchpath.Rmatrices = Array{Float64,2}[]
	node.data.branchpath.RmatricesX = Array{Float64,2}[]
	node.data.branchpath.Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	index = 1
	fastsample = true
	if fastsample
		dt = node.branchlength
		key,R,Pi = getQandPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol, 0), aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], dt)
		node.data.branchpath.key = key
		node.data.branchpath.R = R
		node.data.branchpath.P2 = Pi
		Pret = Pi
		push!(node.data.branchpath.Rmatrices, R*dt)
		push!(node.data.branchpath.RmatricesX, R)
    	push!(node.data.branchpath.Pmatrices,Pi)
    	push!(dummytime,0.0)
	else
	    for it in multi_iter
			dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
			key,R,Pi = getQandPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol, 0), aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], dt)
			#println("R ",R)
			#println("Pi ",Pi)
			#println(get(hiddeniter.prevstates, selcol-1, 0),"\t", get(hiddeniter.prevstates, selcol, 0),"\t", aaiter.prevstates[1], "\t", node.data.ratesbranchpath.paths[aacol][end], "\t", dt)
			Pret *= Pi
			push!(node.data.branchpath.Rmatrices, R*dt)
			push!(node.data.branchpath.RmatricesX, R)
	    	push!(node.data.branchpath.Pmatrices,Pi)
	    	push!(dummytime,multi_iter.prevtime)
	    	if index == 1
	    		node.data.branchpath.key , node.data.branchpath.R, node.data.branchpath.P2 =  getQandPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol, 0), aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], node.branchlength)
	    	end
	    	index += 1
		end
		push!(dummytime,1.0)
	end

    tempv = copy(v)
    pushfirst!(vs,v)
    for P in reverse(node.data.branchpath.Pmatrices)
        tempv = P*tempv
    	pushfirst!(vs,copy(tempv))
    end
    popfirst!(vs)

    node.data.branchpath.vs = vs
    node.data.branchpath.time = dummytime
    node.data.branchpath.P = Pret
    return Pret,node.data.branchpath.Pmatrices,vs
end




function felsensteinresample_hidden(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, selcolin::Int, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams, sample::Bool=true, useoldsampling::Bool=false)
	#selcol = findfirst(x -> x == selcolin, cols)
	selcol = selcolin
	likelihoods = ones(Float64, length(nodelist), modelparams.numhiddenstates)*-Inf
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
				lefttransprobs,Pmatrices_left,vs_left = felsensteinhelper_hidden(nodelist[leftchildindex], selcol, cols, aacol, likelihoods[leftchildindex,:], modelparams)
				righttransprobs,Pmatrices_right,vs_right = felsensteinhelper_hidden(nodelist[rightchildindex], selcol, cols, aacol, likelihoods[rightchildindex,:], modelparams)

        		#likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])
        		likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:]).*observationlikelihood(node.data.protein, selcolin, modelparams)

				m = maximum(likelihoods[nodeindex,:])
				likelihoods[nodeindex,:] = likelihoods[nodeindex,:] ./ m
				logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
        		pop!(stack)

        	end
        end
    end

	rootnode = nodelist[1]
	len = length(rootnode.data.branchpath.paths)
	prevh = 0
	if selcol > 1
		prevh = rootnode.data.branchpath.paths[selcol-1][end]
	end
	nexth = 0
	if selcol < len
		nexth = rootnode.data.branchpath.paths[selcol+1][end]
	end
	freqs = gethiddeninitialprobs(modelparams, prevh, nexth, rootnode.data.aabranchpath.paths[selcol][end])
	rootliks = freqs.*likelihoods[1,:]
	Z = log(sum(rootliks)) + logm[1]
	if sample
		#=
		rootstate = CommonUtils.sample(rng,rootliks)
		rootnode.data.branchpath.paths[selcol] = Int[rootstate]
		rootnode.data.branchpath.times[selcol] = Float64[0.0]
		if useoldsampling
			backwardsampling_old(rng,nodelist[1], rootstate, selcol, modelparams)
		else
			return backwardsampling(rng,nodelist[1], rootstate, selcol, modelparams)
			#return backwardsampling_stack(rng, nodelist, rootstate, selcol, modelparams)
		end
		=#
		return rootliks
	else
		return Z
	end
end

function joint_proposal_likelihood(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams, paths::Array{Array{Int,1},1}=Array{Int,1}[], times::Array{Array{Float64,1},1}=Array{Float64,1}[])
	rootnode = nodelist[1]
	len = length(rootnode.data.branchpath.paths)
	prevh = 0
	if aacol > 1
		prevh = rootnode.data.branchpath.paths[aacol-1][end]
	end
	nexth = 0
	if aacol < len
		nexth = rootnode.data.branchpath.paths[aacol+1][end]
	end
	aa =  ((rootnode.data.jointbranchpath.paths[aacol][end]-1) % modelparams.alphabet)+1 
	if length(paths) > 0
		aa = ((paths[1][end]-1) % modelparams.alphabet)+1 	
	end
	freqs = gethiddeninitialprobs(modelparams, prevh, nexth, aa)
	loglikelihood = 0.0
	for (nodeindex, node) in enumerate(nodelist)
		path = node.data.jointbranchpath.paths[aacol]
		time = node.data.jointbranchpath.times[aacol]
		if length(paths) > 0
			path = paths[nodeindex]
			time = times[nodeindex]
		end
		if isroot(node)
			h = div(path[end]-1, modelparams.alphabet)+1 
			loglikelihood += log(freqs[h])
		else	
			loglikelihood += log(node.data.jointbranchpath.P[path[1],path[end]])
			
			R = node.data.jointbranchpath.R
			for i=1:length(path)-1
				dt = (time[i+1]-time[i])*node.branchlength
				loglikelihood += R[path[i],path[i]]*dt
				loglikelihood += log(R[path[i],path[i+1]]*node.branchlength)
			end
			dt = (1.0 - time[end])*node.branchlength
			loglikelihood += R[path[end],path[end]]*dt
			loglikelihood -= log(node.data.jointbranchpath.P2[path[1],path[end]])
		end
	end
	return loglikelihood
end

function hidden_proposal_likelihood(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams, paths::Array{Array{Int,1},1}=Array{Int,1}[], times::Array{Array{Float64,1},1}=Array{Float64,1}[])
	rootnode = nodelist[1]
	len = length(rootnode.data.branchpath.paths)
	prevh = 0
	if aacol > 1
		prevh = rootnode.data.branchpath.paths[aacol-1][end]
	end
	nexth = 0
	if aacol < len
		nexth = rootnode.data.branchpath.paths[aacol+1][end]
	end
	freqs = gethiddeninitialprobs(modelparams, prevh, nexth, rootnode.data.aabranchpath.paths[aacol][end])	
	loglikelihood = 0.0
	for (nodeindex, node) in enumerate(nodelist)
		path = node.data.branchpath.paths[aacol]
		time = node.data.branchpath.times[aacol]
		if length(paths) > 0
			path = paths[nodeindex]
			time = times[nodeindex]
		end
		if isroot(node)
			h = path[end]
			loglikelihood += log(freqs[h])
		else	
			loglikelihood += log(node.data.branchpath.P[path[1],path[end]])
			
			R = node.data.branchpath.R
			for i=1:length(path)-1
				dt = (time[i+1]-time[i])*node.branchlength
				loglikelihood += R[path[i],path[i]]*dt
				loglikelihood += log(R[path[i],path[i+1]]*node.branchlength)
			end
			dt = (1.0 - time[end])*node.branchlength
			loglikelihood += R[path[end],path[end]]*dt
			loglikelihood -= log(node.data.branchpath.P2[path[1],path[end]])
		end
	end
	return loglikelihood
end

function backwardsampling(rng::AbstractRNG, nodelist::Array{TreeNode,1}, rootstate::Int, aacol::Int,modelparams::ModelParams)
	nodelist[1].data.branchpath.paths[aacol] = Int[rootstate]
	nodelist[1].data.branchpath.times[aacol] = Float64[0.0]
	states = zeros(Int, length(nodelist))
	states[1] = rootstate
	for node in nodelist
		state = states[node.nodeindex]
		for child in node
			if states[child.nodeindex] == 0
				P = child.data.branchpath.P
				R = child.data.branchpath.R
				liks = P[state,:].*child.data.branchpath.vs[end]
				b = CommonUtils.sample(rng,liks)
				states[child.nodeindex] = b
				newpath, newtime, success = modifiedrejectionsampling(rng, R*child.branchlength, state, b, modelparams,child.data.branchpath.key)
				if !success
					return false
				end
				child.data.branchpath.paths[aacol] = newpath
				child.data.branchpath.times[aacol] = newtime
			end
		end
	end
	return true
end


#=
function backwardsampling(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		P = child.data.branchpath.P
		R = child.data.branchpath.R
		liks = P[state,:].*child.data.branchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime, success = modifiedrejectionsampling(rng, R*child.branchlength, state, b, modelparams,child.data.branchpath.key)
		#println(state,"\t",P[state,:],"\t",child.data.branchpath.vs[end],"\t",liks,"\t",b, "\t", newpath,"\t",liks[b])
		if !success
			return false
		end
		child.data.branchpath.paths[aacol] = newpath
		child.data.branchpath.times[aacol] = newtime
		success = backwardsampling(rng,child, b,aacol, modelparams)
		if !success
			return false
		end
	end
	return true
end=#

function backwardsampling_loglikelihood(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams)
	rootnode = nodelist[1]
	len = length(rootnode.data.branchpath.paths)
	prevh = 0
	if aacol > 1
		prevh = rootnode.data.branchpath.paths[aacol-1][end]
	end
	nexth = 0
	if aacol < len
		nexth = rootnode.data.branchpath.paths[aacol+1][end]
	end
	freqs = gethiddeninitialprobs(modelparams, prevh, nexth, rootnode.data.aabranchpath.paths[aacol][end])
	h = rootnode.data.branchpath.paths[aacol][end]
	ll = log(freqs[h]*likelihoods[1,h])
	for node in nodelist
		if !isroot(node)
			a = node.data.branchpath.paths[aacol][1]
			b = node.data.branchpath.paths[aacol][end]
			ll += log(node.data.branchpath.P[a,b]*child.data.branchpath.vs[b])
		end
	end
	return ll
end

function proposallikelihood_stack(nodelist::Array{TreeNode,1}, selcol::Int, modelparams::ModelParams, inpaths::Array{Array{Int,1},1}=Array{Int,1}[], intimes::Array{Array{Float64,1},1}=Array{Float64,1}[])	
	rootnode = nodelist[1]
	len = length(rootnode.data.branchpath.paths)
	prevh = 0
	if selcol > 1
		prevh = rootnode.data.branchpath.paths[selcol-1][end]
	end
	nexth = 0
	if selcol < len
		nexth = rootnode.data.branchpath.paths[selcol+1][end]
	end
	freqs = gethiddeninitialprobs(modelparams, prevh, nexth, rootnode.data.aabranchpath.paths[selcol][end])	
	loglikelihood = 0.0
	for node in nodelist
		time = node.data.branchpath.times[selcol]
		path = node.data.branchpath.paths[selcol]
		if length(intimes) > 0
			time = intimes[node.nodeindex]
			path = inpaths[node.nodeindex]
		end

		if isroot(node)
			h = path[end]
			loglikelihood += log(freqs[h])
		else
			combinedpaths = Int[path[1]]
			combinedtimes = Float64[0.0]
			Rindices = Int[1]
			i1 = 2
			i2 = 2			
			while true
				if i1 <= length(time) && i2 <= length(node.data.branchpath.time)-1
					t1 = time[i1]
					t2 = node.data.branchpath.time[i2]
					if t1 < t2
						push!(combinedpaths, path[i1])
						t1 = time[i1]		
						push!(combinedtimes, t1)	
						push!(Rindices, Rindices[end])
						i1 += 1
					elseif t2 < t1
						push!(combinedpaths, combinedpaths[end])
						t2 = node.data.branchpath.time[i2]
						push!(combinedtimes, t2)
						push!(Rindices,i2)
						i2 += 1
					else
						push!(combinedpaths, path[i1])
						t1 = time[i1]		
						push!(combinedtimes, t1)
						push!(Rindices, Rindices[end])			
						i1 += 1
						i2 += 1
					end
				elseif i1 <= length(time)	
					push!(combinedpaths, path[i1])				
					t1 = time[i1]		
					push!(combinedtimes, t1)	
					push!(Rindices, Rindices[end])			
					i1 += 1
				elseif i2 <= length(node.data.branchpath.time)-1
					push!(combinedpaths, combinedpaths[end])
					t2 = node.data.branchpath.time[i2]
					push!(combinedtimes, t2)
					push!(Rindices,i2)
					i2 += 1
				else
					break
				end
			end
			#=
			if length(times) > 2 && length(node.data.branchpath.time)-1 > 2
				println("A ", times,"\t", paths)
				println("B ", node.data.branchpath.time)
				println("C ", combinedtimes, "\t", combinedpaths)
				println("D ", Rindices)
				println("E ", length(intimes))
			end=#

			for z=1:length(Rindices)-1
				Rindex = Rindices[z]
				dt = combinedtimes[z+1] - combinedtimes[z]
				a = combinedpaths[z]
				b = combinedpaths[z+1]
				if a != b
					loglikelihood += log(node.data.branchpath.RmatricesX[Rindex][a,b]*node.branchlength)
				end
				loglikelihood += node.data.branchpath.RmatricesX[Rindex][a,a]*dt*node.branchlength
			end
			Rindex = Rindices[end]
			a = combinedpaths[end]
			dt = 1.0 - combinedtimes[end]
			loglikelihood += node.data.branchpath.RmatricesX[Rindex][a,a]*dt*node.branchlength
			#=
			if length(times[node.nodeindex]) > 2 && length(node.data.branchpath.time)-1 > 2
				
				println("ll= ", loglikelihood)
			end=#
		end
	end
	return loglikelihood
end

function backwardsampling_stack(rng::AbstractRNG, nodelist::Array{TreeNode,1}, rootstate::Int, selcol::Int,modelparams::ModelParams)
	states = zeros(Int, length(nodelist))
	states[1] = rootstate
	for node in nodelist
		state = states[node.nodeindex]
		for child in node
			if states[child.nodeindex] == 0
				path = Int[state]
				for (Pi,v) in zip(child.data.branchpath.Pmatrices, child.data.branchpath.vs)
					liks = Pi[path[end],:].*v			
					samplestate = CommonUtils.sample(rng,liks)
					push!(path,samplestate)			
				end
				states[child.nodeindex] = path[end]

				newpath = Int[]
				newtime = Float64[]
				for z=1:length(path)-1
					dt = child.data.branchpath.time[z+1]-child.data.branchpath.time[z]
					samplepath, sampletimes, success = modifiedrejectionsampling(rng, child.data.branchpath.Rmatrices[z], path[z], path[z+1], modelparams)
					if !success
						return false
					end
					append!(newpath,samplepath)
					append!(newtime,(sampletimes*dt) .+ child.data.branchpath.time[z])
				end

				newpath, newtime = removevirtualjumps(newpath, newtime)
				child.data.branchpath.paths[selcol] = newpath
				child.data.branchpath.times[selcol] = newtime
			end
		end
	end
	#println("states ",states)
	return true
end

function felsensteinhelper_aa_likelihood(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
	cols = Int[aacol]
	loglikelihood = log(modelparams.hiddennodes[nodelist[1].data.branchpath.paths[aacol][end]].aa_node.probs[nodelist[1].data.aabranchpath.paths[aacol][end]])
	for node in nodelist
		if !isroot(node)
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols),BranchPathIterator(node.data.aabranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]

		    for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength	

				loglikelihood += getaaentry(modelparams, 0, 0, hiddeniter.prevstates[1],  aaiter.prevstates[1], aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end])*dt

				changecol = 0
				if multi_iter.branchpathindex == 1
					changecol = hiddeniter.mincol
				else
					changecol = aaiter.mincol
				end
				if multi_iter.branchpathindex == 2
					loglikelihood += log(getaaentry(modelparams, 0, 0, hiddeniter.prevstates[1], aaiter.prevstates[1], aaiter.currstates[1], node.data.ratesbranchpath.paths[aacol][end])*node.branchlength)
				end
			end
		end
	end
    return loglikelihood
end

function felsensteinhelper_aa(node::TreeNode, incols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	cols = Int[aacol]
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols)])
	hiddeniter = multi_iter.branchpathiterators[1]
	node.data.aabranchpath.Rmatrices = Array{Float64,2}[]
	node.data.aabranchpath.RmatricesX = Array{Float64,2}[]
	node.data.aabranchpath.Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	Pret = Matrix{Float64}(I, modelparams.alphabet, modelparams.alphabet)
	node.data.aabranchpath.R = Pret
	index = 1
	fastsample = true
	if fastsample
		dt = node.branchlength
		key,R,Pi = getAAandPt(modelparams, 0, 0, hiddeniter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], dt)				
		node.data.aabranchpath.key = key
		node.data.aabranchpath.R = R
		node.data.aabranchpath.P2 = Pi
		Pret = Pi
		push!(node.data.aabranchpath.Rmatrices, R*dt)
		push!(node.data.aabranchpath.RmatricesX, R)
    	push!(node.data.aabranchpath.Pmatrices,Pi)
    	push!(dummytime,0.0)
	else
	    for it in multi_iter
			dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
			key,R,Pi = getAAandPt(modelparams, 0, 0, hiddeniter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], dt)				
	    	Pret *= Pi
			push!(node.data.aabranchpath.Rmatrices, R*dt)
			push!(node.data.aabranchpath.RmatricesX, R)
	    	push!(node.data.aabranchpath.Pmatrices,Pi)
	    	push!(dummytime,multi_iter.prevtime)
	    	if index == 1
	    		node.data.aabranchpath.key,node.data.aabranchpath.R, node.data.aabranchpath.P2 = getAAandPt(modelparams, 0, 0, hiddeniter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], node.branchlength)
	    	end
	    	index += 1
		end
	end
	push!(dummytime,1.0)

    tempv = copy(v)
    pushfirst!(vs,v)
    for P in reverse(node.data.aabranchpath.Pmatrices)    	
        tempv = P*tempv
    	pushfirst!(vs,copy(tempv))
    end
    popfirst!(vs)

    node.data.aabranchpath.vs = vs
    node.data.aabranchpath.time = dummytime
    node.data.aabranchpath.P = Pret
    return Pret,node.data.aabranchpath.Pmatrices,vs
end

function felsensteinresample_aa(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams, sample::Bool=true, useoldsampling::Bool=false)
	likelihoods = ones(Float64, length(nodelist), modelparams.alphabet)*-Inf
	logm = zeros(Float64,length(nodelist))

	stack = Int[1]
	while length(stack) > 0
		nodeindex = stack[end]
		node = nodelist[nodeindex]
		if isleafnode(node)
			if 1 <= aacol <= length(node.data.protein.sites) && node.data.protein.sites[aacol].aa > 0
				for a=1:modelparams.alphabet
					likelihoods[nodeindex,a] = 0.0
				end
				likelihoods[nodeindex,node.data.protein.sites[aacol].aa] = 1.0
			else
				if modelparams.hidden_conditional_on_aa					
					likelihoods[nodeindex,:] = observationlikelihood_aa(node.data.protein, aacol, node.data.branchpath.paths[aacol][end], modelparams)
					#=
					for a=1:modelparams.alphabet
						if likelihoods[nodeindex,a] != 1.0
							println(likelihoods)							
							#exit()
						end
					end=#
				else
					for a=1:modelparams.alphabet
						likelihoods[nodeindex,a] = 1.0
					end
				end
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
				lefttransprobs,Pmatrices_left,vs_left = felsensteinhelper_aa(nodelist[leftchildindex], cols, aacol, likelihoods[leftchildindex,:], modelparams)
				righttransprobs,Pmatrices_right,vs_right = felsensteinhelper_aa(nodelist[rightchildindex], cols, aacol, likelihoods[rightchildindex,:], modelparams)
        		likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])
        		if modelparams.hidden_conditional_on_aa
        			# TODO add something here for internal nodes!!!!
        		end

        		m = maximum(likelihoods[nodeindex,:])
				likelihoods[nodeindex,:] = likelihoods[nodeindex,:] ./ m
				logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
        		pop!(stack)
        	end
        end
    end

	rootnode = nodelist[1]
	len = length(rootnode.data.aabranchpath.paths)
	h = rootnode.data.branchpath.paths[aacol][end]
	rootliks = modelparams.hiddennodes[h].aa_node.probs.*likelihoods[1,:]
	if sample
		#=
		rootstate = CommonUtils.sample(rng,rootliks)
		rootnode.data.aabranchpath.paths[aacol] = Int[rootstate]
		rootnode.data.aabranchpath.times[aacol] = Float64[0.0]
		if useoldsampling
			return backwardsampling_aa_old(rng,nodelist[1], rootstate, aacol,modelparams)
		else
			return backwardsampling_aa(rng,nodelist[1], rootstate, aacol,modelparams)
		end
		=#
		return rootliks
	else
		return log(sum(rootliks)) + logm[1]
	end
end


function aa_proposal_likelihood(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams, paths::Array{Array{Int,1},1}=Array{Int,1}[], times::Array{Array{Float64,1},1}=Array{Float64,1}[])	
	loglikelihood = 0.0
	for (nodeindex, node) in enumerate(nodelist)
		path = node.data.aabranchpath.paths[aacol]
		time = node.data.aabranchpath.times[aacol]
		if length(paths) > 0
			path = paths[nodeindex]
			time = times[nodeindex]
		end
		if isroot(node)
			h = node.data.branchpath.paths[aacol][end]
			loglikelihood += log(modelparams.hiddennodes[h].aa_node.probs[path[end]])
		else
			P = node.data.aabranchpath.P
			R = node.data.aabranchpath.R			
			loglikelihood += log(P[path[1],path[end]])
			
			for i=1:length(path)-1
				loglikelihood += R[path[i],path[i]]*(time[i+1]-time[i])*node.branchlength
				loglikelihood += log(R[path[i],path[i+1]]*node.branchlength)
			end
			dt = (1.0 - time[end])*node.branchlength
			loglikelihood += R[path[end],path[end]]*dt
			loglikelihood -= log(node.data.aabranchpath.P2[path[1],path[end]])
		end
	end
	return loglikelihood
end

function backwardsampling_aa_loglikelihood(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams)
	rootnode = nodelist[1]
	len = length(rootnode.data.aabranchpath.paths)
	h = rootnode.data.branchpath.paths[aacol][end]
	aa  = rootnode.data.aabranchpath.paths[aacol][end]
	ll = log(modelparams.hiddennodes[h].aa_node.probs[aa]*likelihoods[1,aa])
	for node in nodelist
		if !isroot(node)
			a = node.data.aabranchpath.paths[aacol][1]
			b = node.data.aabranchpath.paths[aacol][end]
			ll += log(node.data.aabranchpath.P[a,b]*child.data.aabranchpath.vs[b])
		end
	end
	return ll
end

function backwardsampling_aa(rng::AbstractRNG, nodelist::Array{TreeNode,1}, rootstate::Int, aacol::Int,modelparams::ModelParams)
	nodelist[1].data.aabranchpath.paths[aacol] = Int[rootstate]
	nodelist[1].data.aabranchpath.times[aacol] = Float64[0.0]
	states = zeros(Int, length(nodelist))
	states[1] = rootstate
	for node in nodelist
		state = states[node.nodeindex]
		for child in node
			if states[child.nodeindex] == 0
				P = child.data.aabranchpath.P
				R = child.data.aabranchpath.R
				liks = P[state,:].*child.data.aabranchpath.vs[end]
				b = CommonUtils.sample(rng,liks)
				states[child.nodeindex] = b
				newpath, newtime, success = modifiedrejectionsampling(rng, R*child.branchlength, state, b, modelparams,child.data.aabranchpath.key)
				if !success
					return false
				end
				child.data.aabranchpath.paths[aacol] = newpath
				child.data.aabranchpath.times[aacol] = newtime
			end
		end
	end
	return true
end
#=
function backwardsampling_aa(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		P = child.data.aabranchpath.P
		R = child.data.aabranchpath.R
		liks = P[state,:].*child.data.aabranchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime, success = modifiedrejectionsampling(rng, R*child.branchlength, state, b, modelparams, child.data.aabranchpath.key)
		if !success
			return false
		end
		child.data.aabranchpath.paths[aacol] = newpath
		child.data.aabranchpath.times[aacol] = newtime
		success = backwardsampling_aa(rng,child, b,aacol, modelparams)
		if !success
			return false
		end
	end
	return true
end=#

function augmentedloglikelihood_site(nodelist::Array{TreeNode,1}, col::Int, modelparams::ModelParams)
	numcols = length(nodelist[1].data.branchpath.paths)
	loglikelihood = 0.0
	prevh = 0
	if col > 1
		prevh = nodelist[1].data.branchpath.paths[col-1][end]
	end
	#=
	nexth = 0
	if col < numcols
		nexth = nodelist[1].data.branchpath.paths[col+1][end]
	end=#
	h = nodelist[1].data.branchpath.paths[col][end]		
	if prevh > 0
		loglikelihood += log(modelparams.transitionprobs[prevh,h])			
	end
	loglikelihood += log(modelparams.hiddennodes[h].aa_node.probs[nodelist[1].data.aabranchpath.paths[col][end]])
	for node in nodelist
		if !isroot(node)			
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

			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]	

			for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)
				changecol = selcol
				Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol], node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
				loglikelihood += Qii*dt*node.branchlength

				changecol = 0
				if multi_iter.branchpathindex == 1
					changecol = hiddeniter.mincol
				else
					changecol = aaiter.mincol
				end

				if changecol == selcol
					prevstatesh = hiddeniter.prevstates[changecol]
					currstatesh = prevstatesh
					prevstatesaa = aaiter.prevstates[changecol]
					currstatesaa = prevstatesaa
					if multi_iter.branchpathindex == 1
						currstatesh = hiddeniter.currstates[changecol]
					else
						currstatesaa = aaiter.currstates[changecol]
					end

					Qhi = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), prevstatesh, currstatesh, prevstatesaa, currstatesaa, node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
					loglikelihood += log(Qhi*node.branchlength)
				end
			end
		end
	end

	#println("loglikelihood ", loglikelihood)
	return loglikelihood
end

function augmentedloglikelihood(nodelist::Array{TreeNode,1}, inputcols::Array{Int,1}, modelparams::ModelParams)
	numcols = length(inputcols)
	loglikelihood = 0.0

	numcols = length(nodelist[1].data.branchpath.paths)
	for col in inputcols
		prevh = 0
		if col > 1
			prevh = nodelist[1].data.branchpath.paths[col-1][end]
		end
		#=
		nexth = 0		
		if col < numcols
			nexth = nodelist[1].data.branchpath.paths[col+1][end]
		end=#
		h = nodelist[1].data.branchpath.paths[col][end]
		if prevh > 0
			loglikelihood += log(modelparams.transitionprobs[prevh,h])			
		end
		loglikelihood += log(modelparams.hiddennodes[h].aa_node.probs[nodelist[1].data.aabranchpath.paths[col][end]])
	end
	for node in nodelist
		if !isroot(node)
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

				multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
				hiddeniter = multi_iter.branchpathiterators[1]
				aaiter = multi_iter.branchpathiterators[2]	

				for it in multi_iter
					dt = (multi_iter.currtime-multi_iter.prevtime)
					changecol = selcol
					Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol], node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
					loglikelihood += Qii*dt*node.branchlength

					changecol = 0
					if multi_iter.branchpathindex == 1
						changecol = hiddeniter.mincol
					else
						changecol = aaiter.mincol
					end

					if changecol == selcol
						prevstatesh = hiddeniter.prevstates[changecol]
						currstatesh = prevstatesh
						prevstatesaa = aaiter.prevstates[changecol]
						currstatesaa = prevstatesaa
						if multi_iter.branchpathindex == 1
							currstatesh = hiddeniter.currstates[changecol]
						else
							currstatesaa = aaiter.currstates[changecol]
						end

						Qhi = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), prevstatesh, currstatesh, prevstatesaa, currstatesaa, node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
						loglikelihood += log(Qhi*node.branchlength)
					end
				end
			end
		end
	end

	#println("loglikelihood ", loglikelihood)
	return loglikelihood
end



function augmentedloglikelihood_slow(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, modelparams::ModelParams)
	numcols = length(cols)
	loglikelihood = 0.0

	numcols = length(nodelist[1].data.branchpath.paths)
	for col in cols
		prevh = 0
		if col > 1
			prevh = nodelist[1].data.branchpath.paths[col-1][end]
		end
		nexth = 0
		if col < numcols
			nexth = nodelist[1].data.branchpath.paths[col+1][end]
		end		
		if prevh > 0 && nexth > 0
			loglikelihood += log(modelparams.transitionprobs[prevh,nexth])			
		end

		h = nodelist[1].data.branchpath.paths[col][end]
		loglikelihood += log(modelparams.hiddennodes[h].aa_node.probs[nodelist[1].data.aabranchpath.paths[col][end]])
	end
	for node in nodelist
		Qiitotal = 0.0
		Qhitotal = 0.0
		N = 0.0
		if !isroot(node)
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]	

			for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)
				for col=1:numcols
					changecol = col
					Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol])
					loglikelihood += Qii*dt*node.branchlength
				end

				changecol = 0
				if multi_iter.branchpathindex == 1
					changecol = hiddeniter.mincol
				else
					changecol = aaiter.mincol
				end

				if changecol > 0
					prevstatesh = hiddeniter.prevstates[changecol]
					currstatesh = prevstatesh
					prevstatesaa = aaiter.prevstates[changecol]
					currstatesaa = prevstatesaa
					if multi_iter.branchpathindex == 1
						currstatesh = hiddeniter.currstates[changecol]
					else
						currstatesaa = aaiter.currstates[changecol]
					end

					Qhi = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), prevstatesh, currstatesh, prevstatesaa, currstatesaa)
					loglikelihood += log(Qhi*node.branchlength)
				end
			end
		end
	end

	println("HERE", augmentedloglikelihood_fast(nodelist,cols,modelparams))
	println("HERE", loglikelihood)
	return loglikelihood
end

function gethiddeninitialprobs(modelparams::ModelParams, prevh::Int, nexth::Int, aa::Int)
	prevprobs = ones(Float64, modelparams.numhiddenstates)
	if prevh > 0
		prevprobs = modelparams.transitionprobs[prevh,:]
	end	
	nextprobs = ones(Float64, modelparams.numhiddenstates)
	if nexth > 0
		nextprobs = modelparams.transitionprobs[:,nexth]
	end
	logfreqs = zeros(Float64,modelparams.numhiddenstates)
	for h=1:modelparams.numhiddenstates
		logfreqs[h] = log(prevprobs[h]) + log(nextprobs[h]) + log(modelparams.hiddennodes[h].aa_node.probs[aa])
	end
	freqs = exp.(logfreqs.-maximum(logfreqs))
	freqs /= sum(freqs)
end

function getaaentry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, h::Int, prevaa::Int, curraa::Int, sitescale::Float64)
	if prevaa != curraa
		if modelparams.ratemode == 1
			return modelparams.aarates[h]*sitescale*modelparams.scalingfactor*modelparams.mu*modelparams.aa_exchangeablities[prevaa,curraa]*modelparams.hiddennodes[h].aa_node.probs[curraa]
		elseif modelparams.ratemode == 2
			return 0.05*modelparams.aarates[h]*sitescale*modelparams.scalingfactor*modelparams.mu*modelparams.aa_exchangeablities[prevaa,curraa]*sqrt(modelparams.hiddennodes[h].aa_node.probs[curraa]/modelparams.hiddennodes[h].aa_node.probs[prevaa])
		end
	else
		q = 0.0
		for aa=1:modelparams.alphabet
			if aa != prevaa
				q -= getaaentry(modelparams,prevh_hmm,nexth_hmm, h, prevaa, aa, sitescale)
			end
		end
		return q
	end
end

function getaaentry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, h::Int, prevaa::Int, curraa::Int, rate::Int)
	return getaaentry(modelparams, prevh_hmm, nexth_hmm, h, prevaa, curraa, modelparams.rates[rate])
end

function gethiddenentry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, prevh::Int, currh::Int, aa::Int, sitescale::Float64)
	if prevh != currh
		prevprob = 1.0
		prevprob2 = 1.0
		if prevh_hmm != 0
			prevprob = modelparams.transitionprobs[prevh_hmm,currh]
			prevprob2 = modelparams.transitionprobs[prevh_hmm,prevh]
		end
		nextprob = 1.0
		nextprob2 = 1.0
		if nexth_hmm != 0
			nextprob = modelparams.transitionprobs[currh,nexth_hmm]
			nextprob2 = modelparams.transitionprobs[prevh,nexth_hmm]
		end
		if modelparams.ratemode == 1
			return sitescale*modelparams.scalingfactor*sqrt((prevprob*nextprob)/(prevprob2*nextprob2))*modelparams.hiddenmu*modelparams.transitionrates[prevh,currh]*modelparams.hiddennodes[currh].aa_node.probs[aa]
			#return modelparams.scalingfactor*sqrt((prevprob*nextprob)/(prevprob2*nextprob2))*modelparams.hiddenmu*modelparams.transitionrates[prevh,currh]*modelparams.hiddennodes[currh].aa_node.probs[aa]
		elseif modelparams.ratemode == 2
			return 0.05*sitescale*modelparams.scalingfactor*sqrt((prevprob*nextprob)/(prevprob2*nextprob2))*modelparams.hiddenmu*modelparams.transitionrates[prevh,currh]*sqrt(modelparams.hiddennodes[currh].aa_node.probs[aa]/modelparams.hiddennodes[prevh].aa_node.probs[aa])
		end
	else
		q = 0.0
		for h=1:modelparams.numhiddenstates
			if h != prevh
				q -= gethiddenentry(modelparams, prevh_hmm, nexth_hmm, prevh, h, aa, sitescale)
			end
		end
		return q
	end
end 

function gethiddenentry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, prevh::Int, currh::Int, aa::Int, rate::Int)
	return gethiddenentry(modelparams, prevh_hmm, nexth_hmm, prevh, currh, aa, modelparams.rates[rate])
end

function entry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, prevh::Int, currh::Int, prevaa::Int, curraa::Int, prevrate::Int, currrate::Int)
	if prevh == currh && prevaa == curraa
		q = 0.0
		for h=1:modelparams.numhiddenstates
			for aa=1:modelparams.alphabet
				if h != prevh || aa != prevaa
					q -= entry(modelparams, prevh_hmm, nexth_hmm, prevh, h, prevaa, aa, prevrate, currrate)
				end
			end
		end
		return q
	elseif prevh == currh
		return getaaentry(modelparams, prevh_hmm, nexth_hmm, prevh, prevaa, curraa, prevrate)
	elseif prevaa == curraa
		return gethiddenentry(modelparams, prevh_hmm, nexth_hmm, prevh, currh, prevaa, prevrate)
	else
		return 0.0
	end
end



function constructJointMatrix(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int)
	J = zeros(modelparams.numhiddenstates*modelparams.alphabet,modelparams.numhiddenstates*modelparams.alphabet)
	for h1=1:modelparams.numhiddenstates
		for h2=1:modelparams.numhiddenstates
			for aa1=1:modelparams.alphabet
				for aa2=1:modelparams.alphabet
					if (h1 == h2 && aa1 != aa2) || (h1 != h2 && aa1 == aa2)
						i1 = (h1-1)*modelparams.alphabet + aa1
						i2 = (h2-1)*modelparams.alphabet + aa2
						J[i1,i2] = entry(modelparams, prevh_hmm, nexth_hmm, h1, h2, aa1, aa2, 1, 1)
						J[i1,i1] -= J[i1,i2]
					end
				end
			end
		end
	end
	return J
end


function constructAAMatrix(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, h::Int, sitescale::Float64=1.0)
	Q = zeros(Float64, modelparams.alphabet, modelparams.alphabet)
	for aa1=1:modelparams.alphabet
        for aa2=1:modelparams.alphabet
			if aa1 != aa2
				Q[aa1,aa2] = getaaentry(modelparams, prevh_hmm, nexth_hmm, h, aa1, aa2, sitescale)
				Q[aa1,aa1] -= Q[aa1,aa2]
			end
		end
	end

    return Q
end

function constructHiddenMatrix(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, aa::Int, sitescale::Float64=1.0)
	Q = zeros(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)
	for h1=1:modelparams.numhiddenstates
        for h2=1:modelparams.numhiddenstates
			if h1 != h2
				Q[h1,h2] = gethiddenentry(modelparams, prevh_hmm, nexth_hmm, h1, h2, aa, sitescale)
				Q[h1,h1] -= Q[h1,h2]
			end
		end
	end

    return Q
end


function backwardsamplesingle(rng::AbstractRNG, node::TreeNode, modelparams::ModelParams)
	""" use backwardsampling to efficiently sample a tree consisting of a single node """
	numcols = length(node.data.protein.sites)
	liks = zeros(Float64, numcols, modelparams.numhiddenstates)
	siteliks = zeros(Float64, modelparams.numhiddenstates)
	loglikelihood = 0.0
	for col=1:numcols
		m = -Inf
		for h=1:modelparams.numhiddenstates
			aa = node.data.protein.sites[col].aa
			siteliks[h] = siteloglikelihood(node.data.protein.sites[col], h, aa, modelparams)
			if aa > 0
				 siteliks[h] += log(modelparams.hiddennodes[h].aa_node.probs[aa])
			end

			m = max(m, siteliks[h])
		end
		loglikelihood += m
		siteliks = exp.(siteliks .- m)

		if col == 1
			for h=1:modelparams.numhiddenstates
				liks[col,h] = siteliks[h]
			end
		else
			res = (((liks[col-1,:]')*modelparams.transitionprobs)').*siteliks
			m = -Inf
			for h=1:modelparams.numhiddenstates 
				liks[col,h] = res[h]
				m = max(m, res[h])
			end
			loglikelihood += log(m)
			liks[col,:] = liks[col,:] ./ m
		end
	end
	loglikelihood += log(sum(liks[numcols,:]))

	col = numcols
	sampledstates = zeros(Int, numcols)
	sampledstates[col] = CommonUtils.sample(rng, liks[col,:])
	while col > 1
		col -= 1
		sampledstates[col] = CommonUtils.sample(rng, liks[col,:].*modelparams.transitionprobs[:,sampledstates[col+1]])
	end

	for col=1:numcols
		node.data.branchpath.paths[col][end] = sampledstates[col]
	end

	return loglikelihood
end

function forwardbackward(rng::AbstractRNG, node::TreeNode, modelparams::ModelParams)
	numcols = length(node.data.protein.sites)
	forwardliks = ones(Float64, numcols, modelparams.numhiddenstates)*-Inf	
	for h=1:modelparams.numhiddenstates
		forwardliks[1,h] = siteloglikelihood(node.data.protein.sites[1], h, node.data.protein.sites[1].aa, modelparams)
	end
	for col=2:numcols
		for prevh=1:modelparams.numhiddenstates
			for h=1:modelparams.numhiddenstates
				aa = node.data.protein.sites[col].aa
				ll = siteloglikelihood(node.data.protein.sites[col], h, aa, modelparams)
				if aa > 0
					 ll += log(modelparams.hiddennodes[h].aa_node.probs[aa])
				end				
				forwardliks[col,h] = logsumexp(forwardliks[col,h], forwardliks[col-1,prevh] + log(modelparams.transitionprobs[prevh,h]) + ll)
			end
		end
	end

	backwardliks = ones(Float64, numcols, modelparams.numhiddenstates)*-Inf
	for h=1:modelparams.numhiddenstates
		backwardliks[numcols,h] = 0.0
	end
	for col=numcols-1:-1:1
		for h=1:modelparams.numhiddenstates
			for nexth=1:modelparams.numhiddenstates			
				aa = node.data.protein.sites[col+1].aa
				ll = siteloglikelihood(node.data.protein.sites[col+1], nexth, aa, modelparams)
				if aa > 0
					 ll += log(modelparams.hiddennodes[nexth].aa_node.probs[aa])
				end
				backwardliks[col,h] = logsumexp(backwardliks[col,h], backwardliks[col+1,nexth] + log(modelparams.transitionprobs[h,nexth]) + ll)				
			end
		end
	end

	marginalprobabilities = zeros(Float64, numcols, modelparams.numhiddenstates)
	for col=1:numcols
		sumll = -Inf
		for h=1:modelparams.numhiddenstates
			marginalprobabilities[col,h] = forwardliks[col,h]+backwardliks[col,h]
			sumll = logsumexp(sumll, marginalprobabilities[col,h])
		end
		for h=1:modelparams.numhiddenstates
			marginalprobabilities[col,h] = exp(marginalprobabilities[col,h]-sumll)
		end
	end

	return marginalprobabilities
end


function siteloglikelihood(site::SiteObservation, h::Int, aa::Int, modelparams::ModelParams)
	ll = 0.0
	if modelparams.usestructureobservations
		if modelparams.hidden_conditional_on_aa
			if aa > 0
				if modelparams.use_bivariate_von_mises && site.phi > -100.0 && site.psi > -100.0
					ll += BivariateVonMises.logpdf(modelparams.hiddennodes[h].phipsi_nodes[aa], Float64[site.phi, site.psi])
				else
					if site.phi > -100.0
						ll += logpdf(modelparams.hiddennodes[h].phi_nodes[aa].dist, site.phi)
					end
					if site.psi > -100.0
						ll += logpdf(modelparams.hiddennodes[h].psi_nodes[aa].dist, site.psi)
					end
				end
				if site.omega > -100.0
					ll += logpdf(modelparams.hiddennodes[h].omega_nodes[aa].dist, site.omega)
				end
			else
				sumll = -Inf
				for aa=1:modelparams.alphabet
					temp = 0.0
					if modelparams.use_bivariate_von_mises && site.phi > -100.0 && site.psi > -100.0
						temp += BivariateVonMises.logpdf(modelparams.hiddennodes[h].phipsi_nodes[aa], Float64[site.phi, site.psi])
					else
						if site.phi > -100.0
							temp += logpdf(modelparams.hiddennodes[h].phi_nodes[aa].dist, site.phi)
						end
						if site.psi > -100.0
							temp += logpdf(modelparams.hiddennodes[h].psi_nodes[aa].dist, site.psi)
						end
					end
					if site.omega > -100.0
						temp += logpdf(modelparams.hiddennodes[h].omega_nodes[aa].dist, site.omega)
					end
					
					sumll = logsumexp(sumll, log(modelparams.hiddennodes[h].aa_node.probs[aa])+temp)
				end
				ll += sumll
			end
		else
			if modelparams.use_bivariate_von_mises && site.phi > -100.0 && site.psi > -100.0
				ll += BivariateVonMises.logpdf(modelparams.hiddennodes[h].phipsi_node, Float64[site.phi, site.psi])
			else
				if site.phi > -100.0
					ll += logpdf(modelparams.hiddennodes[h].phi_node.dist, site.phi)
				end				
				if site.psi > -100.0
					ll += logpdf(modelparams.hiddennodes[h].psi_node.dist, site.psi)
				end
			end			
			if site.omega > -100.0
				ll += logpdf(modelparams.hiddennodes[h].omega_node.dist, site.omega)
			end
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
	end
	return ll
end

function observationlikelihood(protein::Protein, col::Int, modelparams::ModelParams)
	if 1 <= col <= length(protein.sites)
	    v = zeros(Float64, modelparams.numhiddenstates)
	    Z = -Inf
		for h=1:modelparams.numhiddenstates
			v[h] = siteloglikelihood(protein.sites[col], h, protein.sites[col].aa, modelparams)
			Z = CommonUtils.logsumexp(Z, v[h])
		end
		return exp.(v .- Z)
	else
		return ones(Float64, modelparams.numhiddenstates)
	end
end

function observationlikelihood_aa(protein::Protein, col::Int, h::Int, modelparams::ModelParams)
	if 1 <= col <= length(protein.sites)
	    v = zeros(Float64, modelparams.alphabet)
	    #=
	    Z = -Inf
		for aa=1:modelparams.alphabet
			v[aa] = siteloglikelihood(protein.sites[col], h, aa, modelparams)
			Z = CommonUtils.logsumexp(Z, v[aa])
		end
		return exp.(v .- Z)=#
		for aa=1:modelparams.alphabet
			v[aa] = siteloglikelihood(protein.sites[col], h, aa, modelparams)
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
			ll += siteloglikelihood(proteins[node.seqindex].sites[col], h, proteins[node.seqindex].sites[col].aa, modelparams)
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
		site = SiteObservation()
		site.phi = -1000.0
		site.psi = -1000.0
		if aa != '-'
			#site.aa = indexof(string(aa), aminoacids)
			site.h =  node.data.branchpath.paths[col][end]
			h = site.h
			hiddennode = modelparams.hiddennodes[h]
			#site.aa = CommonUtils.sample(rng, hiddennode.aa_node.probs)
			site.aa = node.data.aabranchpath.paths[col][end]
			if modelparams.hidden_conditional_on_aa
				if modelparams.use_bivariate_von_mises
					site.phi, site.psi = BivariateVonMises.sample(rng, hiddennode.phipsi_nodes[site.aa])
				else
					site.phi = pimod(vonmisesrand(rng, hiddennode.phi_nodes[site.aa].dist))
					site.psi = pimod(vonmisesrand(rng, hiddennode.psi_nodes[site.aa].dist))
				end
			else
				if modelparams.use_bivariate_von_mises
					site.phi, site.psi = BivariateVonMises.sample(rng, hiddennode.phipsi_node)
				else
					site.phi = pimod(vonmisesrand(rng, hiddennode.phi_node.dist))
					site.psi = pimod(vonmisesrand(rng, hiddennode.psi_node.dist))
				end
			end
			site.omega = pimod(vonmisesrand(rng, hiddennode.omega_nodes[site.aa].dist))
			bond_lengths = rand(rng, hiddennode.bond_lengths_node.mvn)
			site.bond_length1 = bond_lengths[1]
			site.bond_length2 = bond_lengths[2]
			site.bond_length3 = bond_lengths[3]
			#=
			site.bond_angle1 = pimod(vonmisesrand(rng, hiddennode.bond_angle1_node.dist))
			site.bond_angle2 = pimod(vonmisesrand(rng, hiddennode.bond_angle2_node.dist))
			site.bond_angle3 = pimod(vonmisesrand(rng, hiddennode.bond_angle3_node.dist))
			println("J", col)=#
		end
		push!(sampled.sites, site)
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

using Viridis
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

	gradient = viridisgradient(_viridis_data)

	for o1 in outputlist
		if 0.0 < branch_scalings[o1.nodeindex] < Inf
			val = (log(branch_scalings[o1.nodeindex])+2.0)/4.0
			col = getcolor(gradient, val)
			hex = string(string(col[1], base=16), string(col[2], base=16), string(col[3], base=16))
			if isleafnode(o1)
				o1.name = "'$(o1.name)'[&label=\"$(format(branch_scalings[o1.nodeindex],precision=1))\",!color=#$(hex)]"
			else
				o1.name = "[&label=\"$(format(branch_scalings[o1.nodeindex],precision=1))\",!color=#$(hex)]"			
			end
		end
	end
	newick = getnewick(outputlist[1])

	nexus = """
	#NEXUS
	begin trees;
		tree tree_1 = [&R] $(newick)
	end;

	begin figtree;
		set appearance.backgroundColorAttribute="Default";
		set appearance.backgroundColour=#ffffff;
		set appearance.branchColorAttribute="User selection";
		set appearance.branchColorGradient=false;
		set appearance.branchLineWidth=2.0;
		set appearance.branchMinLineWidth=0.0;
		set appearance.branchWidthAttribute="Fixed";
		set appearance.foregroundColour=#000000;
		set appearance.hilightingGradient=false;
		set appearance.selectionColour=#2d3680;
		set branchLabels.colorAttribute="User selection";
		set branchLabels.displayAttribute="Branch times";
		set branchLabels.fontName="sansserif";
		set branchLabels.fontSize=8;
		set branchLabels.fontStyle=0;
		set branchLabels.isShown=false;
		set branchLabels.significantDigits=4;
		set layout.expansion=0;
		set layout.layoutType="RECTILINEAR";
		set layout.zoom=0;
		set legend.attribute="label";
		set legend.fontSize=10.0;
		set legend.isShown=false;
		set legend.significantDigits=4;
		set nodeBars.barWidth=4.0;
		set nodeBars.displayAttribute=null;
		set nodeBars.isShown=false;
		set nodeLabels.colorAttribute="User selection";
		set nodeLabels.displayAttribute="label";
		set nodeLabels.fontName="sansserif";
		set nodeLabels.fontSize=8;
		set nodeLabels.fontStyle=0;
		set nodeLabels.isShown=true;
		set nodeLabels.significantDigits=4;
		set nodeShapeExternal.colourAttribute="User selection";
		set nodeShapeExternal.isShown=false;
		set nodeShapeExternal.minSize=10.0;
		set nodeShapeExternal.scaleType=Width;
		set nodeShapeExternal.shapeType=Circle;
		set nodeShapeExternal.size=4.0;
		set nodeShapeExternal.sizeAttribute="Fixed";
		set nodeShapeInternal.colourAttribute="User selection";
		set nodeShapeInternal.isShown=false;
		set nodeShapeInternal.minSize=10.0;
		set nodeShapeInternal.scaleType=Width;
		set nodeShapeInternal.shapeType=Circle;
		set nodeShapeInternal.size=4.0;
		set nodeShapeInternal.sizeAttribute="Fixed";
		set polarLayout.alignTipLabels=false;
		set polarLayout.angularRange=0;
		set polarLayout.rootAngle=0;
		set polarLayout.rootLength=100;
		set polarLayout.showRoot=true;
		set radialLayout.spread=0.0;
		set rectilinearLayout.alignTipLabels=false;
		set rectilinearLayout.curvature=0;
		set rectilinearLayout.rootLength=100;
		set scale.offsetAge=0.0;
		set scale.rootAge=1.0;
		set scale.scaleFactor=1.0;
		set scale.scaleRoot=false;
		set scaleAxis.automaticScale=true;
		set scaleAxis.fontSize=8.0;
		set scaleAxis.isShown=false;
		set scaleAxis.lineWidth=1.0;
		set scaleAxis.majorTicks=1.0;
		set scaleAxis.minorTicks=0.5;
		set scaleAxis.origin=0.0;
		set scaleAxis.reverseAxis=false;
		set scaleAxis.showGrid=true;
		set scaleBar.automaticScale=true;
		set scaleBar.fontSize=10.0;
		set scaleBar.isShown=true;
		set scaleBar.lineWidth=1.0;
		set scaleBar.scaleRange=0.0;
		set tipLabels.colorAttribute="User selection";
		set tipLabels.displayAttribute="Names";
		set tipLabels.fontName="sansserif";
		set tipLabels.fontSize=8;
		set tipLabels.fontStyle=0;
		set tipLabels.isShown=true;
		set tipLabels.significantDigits=4;
		set trees.order=false;
		set trees.orderType="increasing";
		set trees.rooting=false;
		set trees.rootingType="User Selection";
		set trees.transform=false;
		set trees.transformType="cladogram";
	end;
	"""
	treewriter = open("tree.scalings.nexus", "w")
	println(treewriter, nexus)
	close(treewriter)

	return outputlist[1]
end



function initialise_tree_aa(rng::AbstractRNG, modelparams::ModelParams, nodelist::Array{TreeNode,1}, numcols::Int) 
	root = nodelist[1]
	reversenodelist = reverse(nodelist)
	node_aa = zeros(Int,length(nodelist),numcols)	
	node_h = zeros(Int,length(nodelist),numcols)
	for node in reversenodelist
		if isleafnode(node)
			for col=1:numcols				
				node_aa[node.nodeindex,col] = rand(rng,1:modelparams.alphabet)
				if 1 <= col <= length(node.data.protein.sites) && node.data.protein.sites[col].aa > 0
					node_aa[node.nodeindex,col] = node.data.protein.sites[col].aa
				end
			end

			logliks = ones(Float64, numcols, modelparams.numhiddenstates)*-Inf
			for h=1:modelparams.numhiddenstates
				logliks[1,h] = log(modelparams.hiddennodes[h].aa_node.probs[node_aa[node.nodeindex,1]])
			end
			for col=2:numcols
				for prevh=1:modelparams.numhiddenstates
					for h=1:modelparams.numhiddenstates
						ll = log(modelparams.hiddennodes[h].aa_node.probs[node_aa[node.nodeindex,col]])
						logliks[col,h] = logsumexp(logliks[col,h], logliks[col-1,prevh] + log(modelparams.transitionprobs[prevh,h]) + ll)
					end
				end
			end

			col = numcols
			node_h[node.nodeindex,col] = CommonUtils.sample(rng, exp.(logliks[col,:].-maximum(logliks[col,:])))
			while col > 1
				col -= 1
				node_h[node.nodeindex,col] = CommonUtils.sample(rng, exp.(logliks[col,:].-maximum(logliks[col,:])).*modelparams.transitionprobs[:,node_h[node.nodeindex,col+1]])
			end

			#=
			v = zeros(Float64, modelparams.numhiddenstates)
			for h=1:modelparams.numhiddenstates
				v[h] = 
			end
			node_h[node.nodeindex,col] = CommonUtils.sample(rng, v)=#
		else				
			for col=1:numcols
				randchildindex = node.children[rand(1:length(node.children))].nodeindex
				node_aa[node.nodeindex,col] = node_aa[randchildindex, col]
				#randchildindex = node.children[rand(1:length(node.children))].nodeindex
				node_h[node.nodeindex,col] = node_h[randchildindex, col]
			end
		end
	end

	paths = Array{Int,1}[]
	times = Array{Float64,1}[]
	for col=1:numcols
		state = node_h[1,col]
		push!(paths,Int[state,state])
		push!(times,Float64[0.0, 1.0])
	end	
	aapaths = Array{Int,1}[]
	aatimes = Array{Float64,1}[]
	for col=1:numcols
		aastate = node_aa[1,col]
		push!(aapaths,Int[aastate,aastate])
		push!(aatimes,Float64[0.0, 1.0])
	end	
	root.data.branchpath = BranchPath(paths,times)
	root.data.aabranchpath = BranchPath(aapaths,aatimes)
	ratespaths = Array{Int,1}[]
	ratestimes = Array{Float64,1}[]
	for col=1:numcols
		push!(ratespaths, Int[div(modelparams.numrates,2)])
		push!(ratestimes, Float64[0.0])
	end
	root.data.ratesbranchpath = BranchPath(ratespaths, ratestimes)

	closest = 1
	for z=1:length(modelparams.rates)
		if abs(modelparams.rates[z]-1.0) < abs(modelparams.rates[closest]-1.0)
			closest = z
		end
	end

	for node in nodelist
		if !isroot(node)			
			parentnode = get(node.parent)

			paths = Array{Int,1}[]
			times = Array{Float64,1}[]
			for col=1:numcols
				if modelparams.numhiddenstates == 1					
					push!(paths,Int[1])
					push!(times,Float64[0.0])
				else
					parentstate = node_h[parentnode.nodeindex,col]
					nodestate = node_h[node.nodeindex,col]
					V = getQandPt(modelparams, 0, 0, node_aa[parentnode.nodeindex,col], closest, 1.0)[2]
					path,time, success = modifiedrejectionsampling(rng, V*max(0.01,node.branchlength), parentstate, nodestate, modelparams)
					push!(paths,path)
					push!(times,time)
				end
			end

			node.data.branchpath = BranchPath(paths,times)

			aapaths = Array{Int,1}[]
			aatimes = Array{Float64,1}[]
			for col=1:numcols
				parentaastate = node_aa[parentnode.nodeindex,col]
				nodeaastate = node_aa[node.nodeindex,col]
				Q = getAAandPt(modelparams, 0, 0, node_h[parentnode.nodeindex,col], closest, 1.0)[2]
				aapath,aatime,success = modifiedrejectionsampling(rng, Q*max(0.01,node.branchlength), parentaastate, nodeaastate, modelparams)
				push!(aapaths,aapath)
				push!(aatimes,aatime)
			end
			node.data.aabranchpath = BranchPath(aapaths,aatimes)

			ratespaths = Array{Int,1}[]
			ratestimes = Array{Float64,1}[]

			for col=1:numcols
				push!(ratespaths, Int[closest])
				push!(ratestimes, Float64[0.0])
			end
			node.data.ratesbranchpath = BranchPath(ratespaths, ratestimes)
		end
	end
	for node in nodelist
		for col=1:numcols
			aa_h_to_joint!(modelparams, node, col)
		end
	end
end

function initialise_tree(rng::AbstractRNG, modelparams::ModelParams, inputroot::TreeNode, name_protein_dict::Dict{String,Tuple{Int64,Protein}}, numcols::Int; midpointroot::Bool=true)	
	binarize!(inputroot)
	nodelist = getnodelist(inputroot)
	for (index,node) in enumerate(nodelist)
		node.nodeindex = index
	end
	root = inputroot
	if midpointroot
		root = gettreefromnewick(getnewick(midpoint_root(nodelist)))
	end
	nodelist = getnodelist(root)
	for (index,node) in enumerate(nodelist)
		node.nodeindex = index
	end	

	for node in nodelist
		node.data = AugmentedNodeData(numcols)
		if haskey(name_protein_dict, node.name)
			proteinindex, protein = name_protein_dict[node.name]
			node.seqindex = proteinindex
			node.data.protein = protein			
		end
		node.data.inputbranchlength = node.branchlength
	end

	initialise_tree_aa(rng, modelparams, nodelist, numcols)

	return root,nodelist
end

function training_example_from_sequence_alignment(rng::AbstractRNG, modelparams::ModelParams, fastafile::String; newickfile=nothing, blindproteins::Array{String,1}=String[])
	sequences = AbstractString[]
    names = AbstractString[]

	FastaIO.FastaReader(fastafile) do fr
        for (desc, seq) in fr
            push!(names,desc)
            push!(sequences, seq)
        end
    end
    
    proteins = Protein[]
	name_protein_dict = Dict{String,Tuple{Int,Protein}}()
	for (p,(sequence,name)) in enumerate(zip(sequences,names))
	    protein = Protein(name)
	    for aa in sequence
		    site = SiteObservation()
		    if name in blindproteins
		    	site.aa = 0
		    else
		    	site.aa = CommonUtils.indexof(string(aa),aminoacids)
		    end
		    push!(protein.sites,site)
	    end
	    name_protein_dict[name] = (p, protein)
	    push!(proteins, protein)
	end

    if newickfile == nothing
	    newickstring, cachefile = Binaries.fasttreeaa(fastafile)
	    root = gettreefromnewick(newickstring)
	    root,nodelist = initialise_tree(rng, modelparams, root, name_protein_dict, length(sequences[1]))
	else
		root = gettreefromnewick(readlines(open(newickfile,"r"))[1])
	    root,nodelist = initialise_tree(rng, modelparams, root, name_protein_dict, length(sequences[1]),midpointroot=false)
	end

	return ( proteins, nodelist, sequences)
end

function training_example_from_json_family(rng::AbstractRNG, modelparams::ModelParams, json_family; blindproteins::Array{String,1}=String[], blindstructures::Array{String,1}=String[])
	root = gettreefromnewick(json_family["newick_tree"])
	binarize!(root)
	nodelist = getnodelist(root)
	for (index,node) in enumerate(nodelist)
		node.nodeindex = index
	end

	numcols = length(json_family["proteins"][1]["aligned_sequence"])

	proteins = Protein[]
	name_protein_dict = Dict{String,Tuple{Int,Protein}}()
	sequences = String[]
	for p=1:length(json_family["proteins"])
	    protein = Protein()
	    name = json_family["proteins"][p]["name"]
	    json_protein = json_family["proteins"][p]
	    push!(sequences, json_protein["aligned_sequence"])
	    for (aa,phi_psi,omega,bond_lengths,bond_angles) in zip(json_protein["aligned_sequence"], json_protein["aligned_phi_psi"], json_protein["aligned_omega"], json_protein["aligned_bond_lengths"], json_protein["aligned_bond_angles"])
	    	if name in blindproteins
	    		push!(protein.sites, SiteObservation())
		    else		    	
		        if name in blindstructures
		        	push!(protein.sites, SiteObservation(0,indexof(string(aa), aminoacids),-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0,-1000.0))
		        else
		        	phi = phi_psi[1]
		    		psi = phi_psi[2]
		        	push!(protein.sites, SiteObservation(0,indexof(string(aa), aminoacids),phi,omega,psi,bond_lengths[1],bond_lengths[2],bond_lengths[3],bond_angles[1],bond_angles[2],bond_angles[3]))
		        end
		    end
	    end
	    name_protein_dict[json_family["proteins"][p]["name"]] = (p, protein)
	    push!(proteins, protein)
	end

	root,nodelist = initialise_tree(rng, modelparams, root, name_protein_dict, numcols, midpointroot=false)
	if length(nodelist[1].children) == 1
		root = nodelist[1].children[1]
		root.parent = Nullable{TreeNode}()
		nodelist = TreeNode[root]
	end

	return (proteins, nodelist,json_family, sequences)
end

function getexitrate(node::TreeNode, inputcols::Array{Int,1}, modelparams::ModelParams)
	numcols = length(inputcols)
	exitrate = 0.0
	N = 0.0
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

		multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col])])
		hiddeniter = multi_iter.branchpathiterators[1]
		aaiter = multi_iter.branchpathiterators[2]	
		for it in multi_iter
			dt = (multi_iter.currtime-multi_iter.prevtime)
			
			changecol = selcol
			Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[1], aaiter.prevstates[1], node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
			exitrate +=  Qii*dt

			changecol = 0
			if multi_iter.branchpathindex == 1
				changecol = hiddeniter.mincol
				if changecol == selcol
					N += 1.0
				end
			else
				changecol = aaiter.mincol
				N += 1.0
			end
		end
	end

	return exitrate, N
end

function proposescalingfactor(rng::AbstractRNG, nodelist::Array{TreeNode,1}, inputcols::Array{Int,1}, modelparams::ModelParams)
	exitrate = 0.0
	N = 0.0
	for node in nodelist 
		if !isroot(node)
			numcols = length(inputcols)
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

				multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col])])
				hiddeniter = multi_iter.branchpathiterators[1]
				aaiter = multi_iter.branchpathiterators[2]
				#ratesiter = multi_iter.branchpathiterators[3]
				for it in multi_iter
					dt = (multi_iter.currtime-multi_iter.prevtime)
					
					changecol = selcol
					Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[1], aaiter.prevstates[1], node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
					exitrate +=  Qii*dt*node.branchlength

					changecol = 0
					if multi_iter.branchpathindex == 1
						changecol = hiddeniter.mincol
						if changecol == selcol
							N += 1.0
						end
					elseif multi_iter.branchpathindex == 2
						changecol = aaiter.mincol
						N += 1.0
					end

					
				end
			end
		end
	end
	alpha = N+1.0
	beta = -exitrate/modelparams.scalingfactor
	#=
	dist = Gamma(alpha, 1.0/beta)	
	scalingfactor = rand(dist)
	propratio = logpdf(dist, modelparams.scalingfactor)-logpdf(dist,scalingfactor)
	println(alpha,"\t",beta)
	println("SCALING ",modelparams.scalingfactor,"\t",scalingfactor)
	modelparams.scalingfactor = scalingfactor
	=#
	propratio = 0.0
	priorscale = 0.20
	priorshape = 1.0/priorscale
	priordist = Gamma(priorshape, priorscale)
	sampledist = Gamma(alpha, 1.0/beta)	
	t =	modelparams.scalingfactor
	if t > 0.0
		for z=1:1000
			samplet = rand(sampledist)
			deltall = logpdf(priordist,samplet) - logpdf(priordist,t)
			if exp(deltall) > rand(rng)
				t = samplet
			end
		end
	end
	println("SCALING ",modelparams.scalingfactor,"\t",t)
	modelparams.scalingfactor = t

	return propratio
end

function maxscalingfactor(rng::AbstractRNG, nodelist::Array{TreeNode,1}, inputcols::Array{Int,1}, modelparams::ModelParams)
	exitrate = 0.0
	N = 0.0
	for node in nodelist 
		if !isroot(node)
			numcols = length(inputcols)
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

				multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col])])
				hiddeniter = multi_iter.branchpathiterators[1]
				aaiter = multi_iter.branchpathiterators[2]
				#ratesiter = multi_iter.branchpathiterators[3]
				for it in multi_iter
					dt = (multi_iter.currtime-multi_iter.prevtime)
					
					changecol = selcol
					Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[1], aaiter.prevstates[1], node.data.ratesbranchpath.paths[col][end], node.data.ratesbranchpath.paths[col][end])
					exitrate +=  Qii*dt*node.branchlength

					changecol = 0
					if multi_iter.branchpathindex == 1
						changecol = hiddeniter.mincol
						if changecol == selcol
							N += 1.0
						end
					elseif multi_iter.branchpathindex == 2
						changecol = aaiter.mincol
						N += 1.0
					end

					
				end
			end
		end
	end
	alpha = N+1.0
	beta = -exitrate/modelparams.scalingfactor

	dist = Gamma(alpha, 1.0/beta)	
	scalingfactor = rand(dist)
	propratio = logpdf(dist, modelparams.scalingfactor)-logpdf(dist,scalingfactor)
	println(alpha,"\t",beta)
	println("SCALING ",modelparams.scalingfactor,"\t",scalingfactor)
	modelparams.scalingfactor = alpha/beta
	return propratio
end

function proposebranchlength(rng::AbstractRNG, node::TreeNode, cols::Array{Int,1}, modelparams::ModelParams)
	totalexitrate, N = getexitrate(node, cols, modelparams)
	
	alpha = N+1.0
	beta = -totalexitrate


	priorscale = 0.20
	#println("branch length: ", node.data.inputbranchlength)
	if node.data.inputbranchlength <= 0.0
		println(node.name,"\t",node.data.inputbranchlength)
	end
	priorshape = max(node.data.inputbranchlength, 1e-3)/priorscale
	priordist = Gamma(priorshape, priorscale)

	sampledist = Gamma(alpha, 1.0/beta)
	
	t =	node.branchlength
	if t > 0.0
		for z=1:500
			samplet = rand(sampledist)
			deltall = logpdf(priordist,samplet) - logpdf(priordist,t)
			if exp(deltall) > rand(rng)
				t = samplet
			end
		end
	end
	
	propratio = 0.0



	#propratio = logpdf(sampledist, node.branchlength)-logpdf(sampledist,t)



	#println(node.name,"\t",totalexitrate,"\t",N,"\t",node.data.inputbranchlength,"\t",mean(dist),"\t",t)
	#=
	t =  log(1.0 - rand(rng))/totalexitrate
	newll = log(-totalexitrate) + totalexitrate*t
	oldll = log(-totalexitrate) + totalexitrate*node.branchlength
	propratio = oldll-newll
	=#

	return t,propratio
end

function setrates(nodelist::Array{TreeNode,1}, col::Int, ratecat::Int)
	for node in nodelist
		node.data.ratesbranchpath.paths[col][end] = ratecat
	end
end


function rates_helper3(x::Array{Float64,1}, data::Array{Float64,1})
	ll = sum(logpdf(Gamma(x[1], 1.0/x[1]), data))
	println(x[1],"\t",ll)
	return ll
end

function optimizerates3(data::Array{Float64,1}, modelparams::ModelParams)
	opt = Opt(:LN_COBYLA,1)
    localObjectiveFunction = ((param, grad) -> rates_helper3(param, data))
    lower = ones(Float64, 1)*1e-2
    upper = ones(Float64, 1)*1e6
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, 100)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, Float64[modelparams.rate_alpha])
	modelparams.rate_alpha = minx[1]
end

function optimize_gamma(rng::AbstractRNG, trainingexamples, modelparams::ModelParams)	
	#=
	for (proteins,nodelist) in trainingexamples
		numcols = length(proteins[1])
		for col=1:numcols
			samplesiterates(rng, col, nodelist, modelparams)
		end
	end=#

	exitrate = 0.0
	N = 0.0
	data = Float64[]
	for (proteins,nodelist) in trainingexamples
		for node in nodelist 
			if !isroot(node)
				numcols = length(proteins[1])
				for col=1:numcols
					exitrate, N = getexitrate(node, Int[col], modelparams)
					exitrate = exitrate/modelparams.rates[node.data.ratesbranchpath.paths[col][end]]
					dist = Gamma(N+1.0, -1.0/exitrate)
					push!(data, rand(dist))
				end
			end
		end
	end

	#=
	var = Statistics.var(data)
	theta = var/mean(data)
	k = var/(theta*theta)


	println("Gamma")
	println(length(data))
	println(mean(data))
	println(var)
	println(k,"\t",theta)
	modelparams.rates = discretizegamma(k, theta, modelparams.numrates)
	println(modelparams.rates)=#
	optimizerates3(data,modelparams)
	modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)
	println(modelparams.rates)
end

function rates_likelihood(trainingexamples, modelparams::ModelParams)	
	ll = 0.0
	logfreqs = log.(modelparams.rate_freqs)
	for (proteins,nodelist) in trainingexamples
		numcols = length(proteins[1])
		for col=1:numcols
			oldr = nodelist[1].data.ratesbranchpath.paths[col][end]
			currll = -Inf
			for r=1:modelparams.numrates
				setrates(nodelist, col, r)
				currll = logsumexp(currll, logfreqs[r] + augmentedloglikelihood_site(nodelist, col, modelparams))
			end
			setrates(nodelist, col, oldr)
			ll += currll
		end
	end
	return ll
end	

function samplesiterates(rng::AbstractRNG, cols::Array{Int,1}, col::Int, nodelist::Array{TreeNode,1}, modelparams::ModelParams)
	likelihoods = zeros(Float64, modelparams.numrates)
	for r=1:modelparams.numrates
		setrates(nodelist, col, r)
		likelihoods[r] = augmentedloglikelihood(nodelist, cols, modelparams)
	end
	likelihoods = exp.(likelihoods .- maximum(likelihoods))
	likelihoods = likelihoods.*modelparams.rate_freqs
	likelihoods /= sum(likelihoods)
	setrates(nodelist, col, CommonUtils.sample(rng,likelihoods))
	return likelihoods
end	

function samplepaths_seperate_new(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false; samplehiddenstates::Bool=true, sampleaminoacids::Bool=true, dosamplesiterates::Bool=false, accept_everything::Bool=false)
	#=
	for node in nodelist
		joint_to_aa_h!(modelparams, node, col)
	end
	=#
	numcols = length(nodelist[1].data.branchpath.paths)
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
	hidden_accepted = false
	aa_accepted = false
	accepted_hidden = 0.0
	accepted_hidden_total = 0.0
	accepted_aa = 0.0
	accepted_aa_total = 0.0
	hiddentime = time()
	if samplehiddenstates
		if modelparams.numhiddenstates == 1
			accepted_hidden += 1.0
			hidden_accepted = true
		else
			#=
			temppath = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
			temptimes = Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]
			site1 = augmentedloglikelihood(nodelist, cols, modelparams)	
			rootliks = felsensteinresample_hidden(rng, proteins, nodelist, col, cols, col, modelparams)
			cont = backwardsampling(rng,nodelist, CommonUtils.sample(rng, rootliks), col, modelparams)
			prop1 = hidden_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
			#prop1 = proposallikelihood_stack(nodelist, col, modelparams, temppath, temptimes)	
			site2 = augmentedloglikelihood(nodelist, cols, modelparams)	
			prop2 = hidden_proposal_likelihood(nodelist, col, modelparams)
			#prop2 = proposallikelihood_stack(nodelist, col, modelparams)
			if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)
				accepted_hidden += 1.0
				hidden_accepted = true
			else
				if accept_everything && cont

				else
					for (index,node) in enumerate(nodelist)
						node.data.branchpath.paths[col] = copy(temppath[index])
						node.data.branchpath.times[col] = copy(temptimes[index])
					end
				end
			end=#
			
			temppath = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
			temptimes = Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]			
			site1 = augmentedloglikelihood(nodelist, cols, modelparams)
			rootliks = felsensteinresample_hidden(rng, proteins, nodelist, col, cols, col, modelparams)				
			prop1 = hidden_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
			for z=1:5
				cont = backwardsampling(rng,nodelist, CommonUtils.sample(rng, rootliks), col, modelparams)
				#cont = backwardsampling(rng,nodelist[1], CommonUtils.sample(rng, rootliks), col, modelparams)				
				site2 = augmentedloglikelihood(nodelist, cols, modelparams)	
				prop2 = hidden_proposal_likelihood(nodelist, col, modelparams)
				if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)
					hidden_accepted = true
					site1 = site2
					prop1 = prop2
					temppath = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
					temptimes = Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]
				else
					if accept_everything && cont
						site1 = site2
						prop1 = prop2
						temppath = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
						temptimes = Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]
					else
						for (index,node) in enumerate(nodelist)
							node.data.branchpath.paths[col] = copy(temppath[index])
							node.data.branchpath.times[col] = copy(temptimes[index])
						end
					end
				end
			end
			if hidden_accepted
				accepted_hidden += 1.0
			end
		end	
		accepted_hidden_total += 1.0
	end
	hiddentime = time() - hiddentime

	aatime = time()
	if sampleaminoacids
		cols = Int[col]
		temppath = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
		temptimes = Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]	
		site1 = augmentedloglikelihood(nodelist, cols, modelparams)
		rootliks = felsensteinresample_aa(rng, proteins, nodelist, cols, col, modelparams)		
		prop1 = aa_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
		for z=1:3
			cont = backwardsampling_aa(rng,nodelist, CommonUtils.sample(rng, rootliks), col, modelparams)
			site2 = augmentedloglikelihood(nodelist, cols, modelparams)
			prop2 = aa_proposal_likelihood(nodelist, col, modelparams)
			if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)			
				aa_accepted = true
				site1 = site2
				prop1 = prop2
				temppath = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
				temptimes = Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]			
			else
				if accept_everything && cont
					site1 = site2
					prop1 = prop2
					temppath = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
					temptimes = Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]			
				else
					for (index,node) in enumerate(nodelist)
						node.data.aabranchpath.paths[col] = copy(temppath[index])
						node.data.aabranchpath.times[col] = copy(temptimes[index])
					end
				end
			end
		end
		if aa_accepted
			accepted_aa += 1.0
		end
		accepted_aa_total += 1.0
	end
	aatime = time() - aatime

	if dosamplesiterates
		samplesiterates(rng, cols, col, nodelist, modelparams)
	end

	#=
	for node in nodelist
		aa_h_to_joint!(modelparams, node, col)
	end=#

	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total, hidden_accepted, aa_accepted, hiddentime, aatime
end

function cosine_similarity(v1::Array{Float64,1}, v2::Array{Float64,1})
	return dot(v1,v2)/(norm(v1)*norm(v2))
end

function calculateloglikelihood(modelparams::ModelParams, trainingexamples::Array{Tuple,1})
	observationll = 0.0
	augmentedll = 0.0
	for (proteins,nodelist) in trainingexamples
		numcols = length(proteins[1])
		augmentedll +=  augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
		observationll += observationloglikelihood(proteins, nodelist, modelparams)
	end
	return augmentedll, observationll
end

function calculateloglikelihood_array(modelparams::ModelParams, trainingexamples::Array{Tuple,1})
	augmented_array = Float64[]
	observation_array = Float64[]
	observationll = 0.0
	augmentedll = 0.0
	for (proteins,nodelist) in trainingexamples
		numcols = length(proteins[1])
		augtemp = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
		augmentedll +=  augtemp
		push!(augmented_array, augtemp)

		obstemp = observationloglikelihood(proteins, nodelist, modelparams)		
		observationll += obstemp
		push!(observation_array, obstemp)
	end
	return augmentedll, observationll, augmented_array, observation_array
end

function aaintegrated_loglikelihood(rng::AbstractRNG, modelparams::ModelParams, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1})
	sequencell = 0.0
	numcols = length(proteins[1])
	for col=1:numcols
		sequencell += felsensteinresample_aa(rng, proteins, nodelist, Int[col], col, modelparams, false)
	end
	return sequencell
end

function aaintegrated_loglikelihood_array(rng::AbstractRNG, modelparams::ModelParams, trainingexamples::Array{Tuple,1})
	sequencell = 0.0
	sequence_array = Float64[]
	for (proteins,nodelist) in trainingexamples
		temp = aaintegrated_loglikelihood(rng, modelparams, proteins, nodelist)
		sequencell += temp
		push!(sequence_array, temp)
	end
	return sequencell, sequence_array
end

function aaintegrated_loglikelihood(rng::AbstractRNG, modelparams::ModelParams, trainingexamples::Array{Tuple,1})
	sequencell = 0.0
	for (proteins,nodelist) in trainingexamples
		sequencell += aaintegrated_loglikelihood(rng, proteins, nodelist)
	end
	return sequencell
end

function count_aminoacid_substitutions(rng::AbstractRNG, modelparams::ModelParams, node::TreeNode)
	aajumps = 0.0
	for aapath in node.data.aabranchpath.paths
		aajumps += length(aapath) - 1.0
	end
	return aajumps
end	

function gethiddenarray(modelparams::ModelParams)
	x = Float64[]
	for h=1:modelparams.numhiddenstates
		for h2=1:modelparams.numhiddenstates
			if h != h2
				push!(x, modelparams.transitionrates[h,h2])
			end
		end
	end
	return x
end

function sethiddenarray(modelparams::ModelParams, x::Array{Float64,1})
	index = 1
	for h=1:modelparams.numhiddenstates
		for h2=1:modelparams.numhiddenstates
			if h != h2
				modelparams.transitionrates[h,h2] = x[index]
				index += 1
			end
		end
	end
	return modelparams
end

function hiddenratehelper(x::Array{Float64,1}, trainingexamples::Array{Tuple,1}, modelparams::ModelParams,extra)
	sethiddenarray(modelparams,x)
	augmentedll = 0.0
	for (proteins,nodelist) in trainingexamples		
		numcols = length(proteins[1])
		augmentedll += augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
	end
	return augmentedll
end

function optimizehiddenrates(trainingexamples::Array{Tuple,1}, modelparams::ModelParams, extra)
	numparams = modelparams.numhiddenstates*modelparams.numhiddenstates-modelparams.numhiddenstates
	opt = Opt(:LN_COBYLA, numparams)
    localObjectiveFunction = ((param, grad) -> hiddenratehelper(param, trainingexamples, modelparams,extra))
    lower = ones(Float64, numparams)*1e-2
    upper = ones(Float64, numparams)*1e4
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, 75)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, gethiddenarray(modelparams))
	return minx
end

function hiddenratescaling_helper(x::Array{Float64,1}, trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	modelparams.mu = x[1]
	augmentedll = 0.0
	for (proteins,nodelist) in trainingexamples	
		numcols = length(proteins[1])
		augmentedll += augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
	end
	ll = augmentedll + -1.0*x[1]
	return ll
end

function optimizehiddenscaling(trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	opt = Opt(:LN_COBYLA, 1)
    localObjectiveFunction = ((param, grad) -> hiddenratescaling_helper(param, trainingexamples, modelparams))
    lower = ones(Float64, 1)*1e-2
    upper = ones(Float64, 1)*1e6
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, 15)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, Float64[modelparams.mu])
	return Float64[modelparams.mu]
end

function rates_helper(x::Array{Float64,1}, trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	modelparams.mu = x[1]
	modelparams.hiddenmu = x[2]
	augmentedll = 0.0
	for (proteins,nodelist) in trainingexamples	
		numcols = length(proteins[1])
		augmentedll += augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
	end
	println(augmentedll,"\t",x)
	return augmentedll
end

function optimizerates(trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	opt = Opt(:LN_COBYLA,2)
    localObjectiveFunction = ((param, grad) -> rates_helper(param, trainingexamples, modelparams))
    lower = ones(Float64, 2)*1e-2
    upper = ones(Float64, 2)*1e6
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, 25)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, Float64[modelparams.mu,modelparams.hiddenmu])
	modelparams.mu = minx[1]
	modelparams.hiddenmu = minx[2]
end



function sampletraininginstances(iter::Int, rng::AbstractRNG, trainingexamples::Array{Tuple,1}, modelparams::ModelParams; maxsamplesperiter::Int=500, sitethreshold::Int=2, dosamplesiterates::Bool=true, samplebranchlengths::Bool=true, family_names::Array{String,1}, accept_everything::Bool=false)
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
			accepted = zeros(Int, numcols)			
			for i=1:maxsamplesperiter					
				randcols = shuffle(rng, Int[i for i=1:numcols])
				for col in randcols
					if accepted[col] < sitethreshold || i % 20 == 0 || (col > 1 && accepted[col-1] < sitethreshold) || (col < numcols && accepted[col+1] < sitethreshold) 
						a1,a2,a3,a4, hidden_accepted, aa_accepted, hiddentime, aatime = samplepaths_seperate_new(rng,col,proteins,nodelist, modelparams, dosamplesiterates=dosamplesiterates, accept_everything=accept_everything)
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
			elapsedtime = time()-starttime
			min_accepted = minimum(accepted)
			if trainingindex % 10 == 0
				println("min_accepted ", min_accepted, " mean is ", mean(accepted))
				println(trainingindex, "\t",totalhiddentime,"\t",totalaatime,"\t", elapsedtime)
			end
		end	
	end
	return trainingexamples,totalbranchlength_output,accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total,totalhiddentime,totalaatime
end