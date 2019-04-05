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

mutable struct Sample 
	name::String
	aasamples::Array{Array{Int,1},1}
	hiddensamples::Array{Array{Int,1},1}
	json_family::Dict{String,Any}
	modelparams::ModelParams

	function Sample(name::String, modelparams::ModelParams)
		new(name, Int[], Int[], Dict{String,Any}(), modelparams)
	end
end



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
            c1.parent = Nullable{TreeNode}(n.children[end])
            c2.parent = Nullable{TreeNode}(n.children[end])
            n.children[end].parent = Nullable{TreeNode}(n)
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
function reset_matrix_cache(modelparams::ModelParams)
	global countcachemisses
	global countcachehits
	println("RESET: CACHE HITS $(countcachehits) / $(countcachehits+countcachemisses) ($(countcachehits / (countcachehits+countcachemisses)))")
	countcachemisses = 0
	countcachehits = 0
	modelparams.matrixcache = Dict{Tuple{Int,Int,Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}()
end



function getAAandPt(modelparams::ModelParams, prev_hmm::Int, next_hmm::Int, h::Int, rate::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = (-1,-1,h,0)
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
		return Q, Q
	else
		return Q, absmat(real(V*Diagonal(exp.(D*t*sitescale))*Vi))
	end
end


function getQandPt(modelparams::ModelParams, prevh::Int, nexth::Int, aa::Int, rate::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = (-1,prevh, nexth, aa)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1

		Q = constructHiddenMatrix(modelparams, prevh, nexth, aa)
		
		
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
		return Q, Q
	else
		return Q, absmat(real(V*Diagonal(exp.(D*t*sitescale))*Vi))
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
		newpath, newtime, success = modifiedrejectionsampling(rng, R, state, b,(modelparams))
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
			samplepath, sampletimes, success = modifiedrejectionsampling(rng, R, path[z], path[z+1], (modelparams))
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

function felsensteinhelper(node::TreeNode, selcolin::Int, incols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
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
    for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		R,Pi = getQandPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol, 0), aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], dt)
		#println("R ",R)
		#println("Pi ",Pi)
		#println(get(hiddeniter.prevstates, selcol-1, 0),"\t", get(hiddeniter.prevstates, selcol, 0),"\t", aaiter.prevstates[1], "\t", node.data.ratesbranchpath.paths[aacol][end], "\t", dt)
		Pret *= Pi
		push!(node.data.branchpath.Rmatrices, R*dt)
		push!(node.data.branchpath.RmatricesX, R)
    	push!(node.data.branchpath.Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)
    	if index == 1
    		node.data.branchpath.R, node.data.branchpath.P2 =  getQandPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol, 0), aaiter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], node.branchlength)
    	end
    	index += 1
	end
	push!(dummytime,1.0)

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




function felsensteinresample(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, selcolin::Int, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams, sample::Bool=true, useoldsampling::Bool=false)
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
				lefttransprobs,Pmatrices_left,vs_left = felsensteinhelper(nodelist[leftchildindex], selcol, cols, aacol, likelihoods[leftchildindex,:], modelparams)
				righttransprobs,Pmatrices_right,vs_right = felsensteinhelper(nodelist[rightchildindex], selcol, cols, aacol, likelihoods[rightchildindex,:], modelparams)

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
		loglikelihoods = log.(likelihoods)
		for i=1:size(likelihoods,1)
			loglikelihoods[i,:] = loglikelihoods[i,:] .+ logm[i]
		end
		=#
		rootstate = CommonUtils.sample(rng,rootliks)
		rootnode.data.branchpath.paths[selcol] = Int[rootstate]
		rootnode.data.branchpath.times[selcol] = Float64[0.0]
		if useoldsampling
			backwardsampling_old(rng,nodelist[1], rootstate, selcol, modelparams)
		else
			return backwardsampling(rng,nodelist[1], rootstate, selcol, modelparams)
		end
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

function backwardsampling(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		P = child.data.branchpath.P
		R = child.data.branchpath.R
		liks = P[state,:].*child.data.branchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime, success = modifiedrejectionsampling(rng, R*child.branchlength, state, b,(modelparams))
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
end

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



function backwardsampling_old(rng::AbstractRNG,node::TreeNode, state::Int, selcol::Int,modelparams::ModelParams)
	for child in node
		path = Int[state]
		for (Pi,v) in zip(child.data.branchpath.Pmatrices, child.data.branchpath.vs)
			liks = Pi[path[end],:].*v			
			samplestate = CommonUtils.sample(rng,liks)
			push!(path,samplestate)			
		end

		newpath = Int[]
		newtime = Float64[]
		for z=1:length(path)-1
			dt = child.data.branchpath.time[z+1]-child.data.branchpath.time[z]
			samplepath, sampletimes, success = modifiedrejectionsampling(rng, child.data.branchpath.Rmatrices[z], path[z], path[z+1],(modelparams))
			if !success
				return false
			end
			append!(newpath,samplepath)
			append!(newtime,(sampletimes*dt) .+ child.data.branchpath.time[z])
		end

		newpath, newtime = removevirtualjumps(newpath, newtime)
		child.data.branchpath.paths[selcol] = newpath
		child.data.branchpath.times[selcol] = newtime
		success = backwardsampling_old(rng,child, path[end],selcol,modelparams)
		if !success
			return false
		end
	end
	return true
end

function hidden_proposal_likelihood_old(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams, paths::Array{Array{Int,1},1}=Array{Int,1}[], times::Array{Array{Float64,1},1}=Array{Float64,1}[])
	newpath = Array{Int,1}[]
	newtimes = Array{Float64,1}[]
	if length(paths) > 0
		for (nodeindex,node) in enumerate(nodelist)			
			node.data.branchpath.paths[aacol] = copy(paths[nodeindex])
			node.data.branchpath.times[aacol] = copy(times[nodeindex])
		end
		newpath  = Array{Int,1}[copy(node.data.branchpath.paths[aacol]) for node in nodelist]
		newtimes =  Array{Float64,1}[copy(node.data.branchpath.times[aacol]) for node in nodelist]
	end

	selcol = findfirst(x -> x == aacol, cols)
	#cols = copy(incols)
	#deleteat!(cols,selcol)
	h = nodelist[1].data.branchpath.paths[aacol][end]
	loglikelihood = log(modelparams.hiddennodes[h].aa_node.probs[nodelist[1].data.aabranchpath.paths[aacol][end]])
	for (nodeindex, node) in enumerate(nodelist)
		if !isroot(node)
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol])])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]
			Rindex = 1
    		for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
				loglikelihood += node.data.branchpath.RmatricesX[Rindex][hiddeniter.prevstates[selcol],hiddeniter.prevstates[selcol]]*dt
				if multi_iter.branchpathindex == 1 && hiddeniter.mincol == selcol
					loglikelihood += log(node.data.branchpath.RmatricesX[Rindex][hiddeniter.prevstates[selcol],hiddeniter.currstates[selcol]]*node.branchlength)
				else
					Rindex += 1
				end				
			end
		end
	end
	if length(paths) > 0
		for (nodeindex,node) in enumerate(nodelist)			
			node.data.branchpath.paths[aacol] = newpath[nodeindex]
			node.data.branchpath.times[aacol] = newtimes[nodeindex]
		end
	end
	return loglikelihood
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
    for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		R,Pi = getAAandPt(modelparams, 0, 0, hiddeniter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], dt)				
    	Pret *= Pi
		push!(node.data.aabranchpath.Rmatrices, R*dt)
		push!(node.data.aabranchpath.RmatricesX, R)
    	push!(node.data.aabranchpath.Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)
    	if index == 1
    		node.data.aabranchpath.R, node.data.aabranchpath.P2 = getAAandPt(modelparams, 0, 0, hiddeniter.prevstates[1], node.data.ratesbranchpath.paths[aacol][end], node.branchlength)
    	end
    	index += 1
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
        			# TODO add something here!!!!
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
		rootstate = CommonUtils.sample(rng,rootliks)
		rootnode.data.aabranchpath.paths[aacol] = Int[rootstate]
		rootnode.data.aabranchpath.times[aacol] = Float64[0.0]
		if useoldsampling
			return backwardsampling_aa_old(rng,nodelist[1], rootstate, aacol,modelparams)
		else
			return backwardsampling_aa(rng,nodelist[1], rootstate, aacol,modelparams)
		end
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

function backwardsampling_aa_old(rng::AbstractRNG,node::TreeNode, state::Int, selcol::Int, modelparams::ModelParams)
	for child in node
		path = Int[state]
		for (Pi,v) in zip(child.data.aabranchpath.Pmatrices, child.data.aabranchpath.vs)
			liks = Pi[path[end],:].*v
			samplestate = CommonUtils.sample(rng,liks)
			push!(path,samplestate)
		end

		newpath = Int[]
		newtime = Float64[]
		for z=1:length(path)-1
			dt = child.data.aabranchpath.time[z+1]-child.data.aabranchpath.time[z]
			samplepath, sampletimes, success = modifiedrejectionsampling(rng, child.data.aabranchpath.Rmatrices[z], path[z], path[z+1],(modelparams))
			if !success
				return false
			end
			append!(newpath,samplepath)
			append!(newtime,(sampletimes*dt) .+ child.data.aabranchpath.time[z])
		end

		newpath, newtime = removevirtualjumps(newpath, newtime)

		child.data.aabranchpath.paths[selcol] = newpath
		child.data.aabranchpath.times[selcol] = newtime
		success = backwardsampling_aa_old(rng,child, path[end],selcol,modelparams)
		if !success
			return false
		end
	end

	return true
end

function aa_proposal_likelihood_old(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams, paths::Array{Array{Int,1},1}=Array{Int,1}[], times::Array{Array{Float64,1},1}=Array{Float64,1}[])
	newpath = Array{Int,1}[]
	newtimes = Array{Float64,1}[]
	if length(paths) > 0
		for (nodeindex,node) in enumerate(nodelist)			
			node.data.aabranchpath.paths[aacol] = copy(paths[nodeindex])
			node.data.aabranchpath.times[aacol] = copy(times[nodeindex])
		end
		newpath  = Array{Int,1}[copy(node.data.aabranchpath.paths[aacol]) for node in nodelist]
		newtimes =  Array{Float64,1}[copy(node.data.aabranchpath.times[aacol]) for node in nodelist]
	end


	h = nodelist[1].data.branchpath.paths[aacol][end]
	loglikelihood = log(modelparams.hiddennodes[h].aa_node.probs[nodelist[1].data.aabranchpath.paths[aacol][end]])
	cols = Int[aacol]
	for (nodeindex, node) in enumerate(nodelist)
		if !isroot(node)
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]
			Rindex = 1
    		for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
				loglikelihood += node.data.aabranchpath.RmatricesX[Rindex][aaiter.prevstates[1],aaiter.prevstates[1]]*dt
				if multi_iter.branchpathindex == 2 && aaiter.mincol === 1	
					loglikelihood += log(node.data.aabranchpath.RmatricesX[Rindex][aaiter.prevstates[1],aaiter.currstates[1]]*node.branchlength)
				else
					Rindex += 1
				end				
			end
		end
	end

	if length(paths) > 0
		for (nodeindex,node) in enumerate(nodelist)			
			node.data.aabranchpath.paths[aacol] = newpath[nodeindex]
			node.data.aabranchpath.times[aacol] = newtimes[nodeindex]
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

function backwardsampling_aa(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		P = child.data.aabranchpath.P
		R = child.data.aabranchpath.R
		liks = P[state,:].*child.data.aabranchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime, success = modifiedrejectionsampling(rng, R*child.branchlength, state, b,(modelparams))
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
end

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

function samplesingle(rng::AbstractRNG, node::TreeNode, modelparams::ModelParams)
	""" use backwardsampling to efficiently sample a tree consisting of a single node """
	numcols = length(node.data.protein.sites)
	logliks = ones(Float64, numcols, modelparams.numhiddenstates)*-Inf	
	for h=1:modelparams.numhiddenstates
		logliks[1,h] = siteloglikelihood(node.data.protein.sites[1], h, node.data.protein.sites[1].aa, modelparams)
	end
	for col=2:numcols
		for prevh=1:modelparams.numhiddenstates
			for h=1:modelparams.numhiddenstates
				ll = siteloglikelihood(node.data.protein.sites[col], h, node.data.protein.sites[col].aa, modelparams)
				logliks[col,h] = logsumexp(logliks[col,h], logliks[col-1,prevh] + log(modelparams.transitionprobs[prevh,h]) + ll)
			end
		end
	end

	col = numcols
	sampledstates = zeros(Int, numcols)
	sampledstates[col] = CommonUtils.sample(rng, exp.(logliks[col,:].-maximum(logliks[col,:])))
	while col > 1
		col -= 1
		sampledstates[col] = CommonUtils.sample(rng, exp.(logliks[col,:].-maximum(logliks[col,:])).*modelparams.transitionprobs[:,sampledstates[col+1]])
	end

	for col=1:numcols
		node.data.branchpath.paths[col][end] = sampledstates[col]
	end	
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
	    Z = -Inf
		for aa=1:modelparams.alphabet
			v[aa] = siteloglikelihood(protein.sites[col], h, aa, modelparams)
			Z = CommonUtils.logsumexp(Z, v[aa])
		end
		return exp.(v .- Z)
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
		if aa != '-'
			site = SiteObservation()
			#site.aa = indexof(string(aa), aminoacids)
			site.h =  node.data.branchpath.paths[col][end]
			h = site.h
			hiddennode = modelparams.hiddennodes[h]
			#site.aa = CommonUtils.sample(rng, hiddennode.aa_node.probs)
			site.aa = node.data.aabranchpath.paths[col][end]
			if modelparams.hidden_conditional_on_aa
				#=
				site.phi = pimod(vonmisesrand(rng, hiddennode.phi_nodes[site.aa].dist))
				site.omega = pimod(vonmisesrand(rng, hiddennode.omega_nodes[site.aa].dist))
				site.psi = pimod(vonmisesrand(rng, hiddennode.psi_nodes[site.aa].dist))
				=#
			else
				#=
				site.phi = pimod(vonmisesrand(rng, hiddennode.phi_node.dist))
				site.omega = pimod(vonmisesrand(rng, hiddennode.omega_node.dist))
				site.psi = pimod(vonmisesrand(rng, hiddennode.psi_node.dist))
				=#
			end
			bond_lengths = rand(rng, hiddennode.bond_lengths_node.mvn)
			site.bond_length1 = bond_lengths[1]
			site.bond_length2 = bond_lengths[2]
			site.bond_length3 = bond_lengths[3]
			#=
			site.bond_angle1 = pimod(vonmisesrand(rng, hiddennode.bond_angle1_node.dist))
			site.bond_angle2 = pimod(vonmisesrand(rng, hiddennode.bond_angle2_node.dist))
			site.bond_angle3 = pimod(vonmisesrand(rng, hiddennode.bond_angle3_node.dist))
			println("J", col)=#
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
					V = getQandPt(modelparams, 0, 0, node_aa[parentnode.nodeindex,col], closest, 1.0)[1]
					path,time, success = modifiedrejectionsampling(rng, V*max(0.01,node.branchlength), parentstate, nodestate, nothing)
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
				Q = getAAandPt(modelparams, 0, 0, node_h[parentnode.nodeindex,col], closest, 1.0)[1]
				aapath,aatime,success = modifiedrejectionsampling(rng, Q*max(0.01,node.branchlength), parentaastate, nodeaastate, nothing)
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

function training_example_from_sequence_alignment(rng::AbstractRNG, modelparams::ModelParams, fastafile::String; newickfile=nothing, blindnodenames::Array{String,1}=String[])
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
		    if name in blindnodenames
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

function training_example_from_json_family(rng::AbstractRNG, modelparams::ModelParams, json_family; blindnodenames::Array{String,1}=String[], blindstructurenames::Array{String,1}=String[])
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
	    	if name in blindnodenames
	    		push!(protein.sites, SiteObservation())
		    else		    	
		        if name in blindstructurenames
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
	return (proteins, nodelist,json_family, sequences)
end

function getexitrate_slow(node::TreeNode, cols::Array{Int,1}, modelparams::ModelParams)
	numcols = length(cols)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]	
	exitrate = 0.0
	N = 0.0
	for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)
		
		for col=1:numcols
			changecol = col
			Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol])
			exitrate +=  Qii*dt
		end

		changecol = 0
		if multi_iter.branchpathindex == 1
			changecol = hiddeniter.mincol
		else
			changecol = aaiter.mincol
		end

		if changecol > 0
			N += 1.0
		end
	end

	
	d1, d2 = getexitrate_fast(node, cols, modelparams)
	println("exit ", exitrate,"\t", d1)
	println("N ", N,"\t", d2)


	return exitrate, N
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

	dist = Gamma(alpha, 1.0/beta)	
	scalingfactor = rand(dist)
	propratio = logpdf(dist, modelparams.scalingfactor)-logpdf(dist,scalingfactor)
	println(alpha,"\t",beta)
	println("SCALING ",modelparams.scalingfactor,"\t",scalingfactor)
	modelparams.scalingfactor = scalingfactor
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

	dist = Gamma(alpha, 1.0/beta)	
	t = rand(dist)
	propratio = logpdf(dist, node.branchlength)-logpdf(dist,t)
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

function samplepaths2(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams)
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
	accepted_hidden = 0.0
	accepted_hidden_total = 0.0
	accepted_aa = 0.0
	accepted_aa_total = 0.0
	
	useoldsampling = false

	temppath = Array{Int,1}[node.data.branchpath.paths[col] for node in nodelist]
	temptimes = Array{Float64,1}[node.data.branchpath.times[col] for node in nodelist]
	site1 = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, false, false)
	cont = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams)
	prop1 = hidden_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
	site2 = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, false, false)
	prop2 = hidden_proposal_likelihood(nodelist, col, modelparams)

	if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)
		accepted_hidden += 1.0
	else
		index = 1
		for node in nodelist
			node.data.branchpath.paths[col] = temppath[index]
			node.data.branchpath.times[col] = temptimes[index]
			index += 1
		end
	end
	accepted_hidden_total += 1.0

	temppath = Array{Int,1}[node.data.aabranchpath.paths[col] for node in nodelist]
	temptimes = Array{Float64,1}[node.data.aabranchpath.times[col] for node in nodelist]
	cont = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, true, true)

	if cont
		accepted_aa += 1.0
	else
		index = 1
		for node in nodelist
			node.data.aabranchpath.paths[col] = temppath[index]
			node.data.aabranchpath.times[col] = temptimes[index]
			if node.nodeindex != index
				println("ERROR :( ", node.nodeindex,"\t", index)
				exit()
			end
			index += 1
		end
	end
	accepted_aa_total += 1.0

	samplesiterates(rng, cols, col, nodelist, modelparams)

	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total
end	

function samplepaths4(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false)
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
	accepted_hidden = 0.0
	accepted_hidden_total = 0.0
	accepted_aa = 0.0
	accepted_aa_total = 0.0
	temppath = Array{Int,1}[node.data.branchpath.paths[col] for node in nodelist]
	temptimes = Array{Float64,1}[node.data.branchpath.times[col] for node in nodelist]
	tempaapath = Array{Int,1}[node.data.aabranchpath.paths[col] for node in nodelist]
	tempaatimes = Array{Float64,1}[node.data.aabranchpath.times[col] for node in nodelist]	
	site1 = augmentedloglikelihood(nodelist, cols, modelparams)
	cont = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams,true,true)	
	if cont
		cont = felsensteinresample_aa(rng, proteins, nodelist, cols,col, modelparams,true,true)
	end

	if cont
		newpath  = Array{Int,1}[node.data.branchpath.paths[col] for node in nodelist]
		newtimes =  Array{Float64,1}[node.data.branchpath.times[col] for node in nodelist]	
		newaapath  = Array{Int,1}[node.data.aabranchpath.paths[col] for node in nodelist]
		newaatimes =  Array{Float64,1}[node.data.aabranchpath.times[col] for node in nodelist]
		site2 = augmentedloglikelihood(nodelist, cols, modelparams)

		for (nodeindex, node) in enumerate(nodelist)
			node.data.branchpath.paths[col] = temppath[nodeindex]
			node.data.branchpath.times[col] = temptimes[nodeindex]
		end
		hiddenZ = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams,false)
		prop1 = hidden_proposal_likelihood_old(nodelist, cols, col, modelparams, temppath, temptimes)-hiddenZ
		hiddenAA = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, false)	
		prop1 += aa_proposal_likelihood_old(nodelist, col, modelparams, tempaapath, tempaatimes)-hiddenAA
		for (nodeindex, node) in enumerate(nodelist)
			node.data.branchpath.paths[col] = newpath[nodeindex]
			node.data.branchpath.times[col] = newtimes[nodeindex]
		end

		
		for (nodeindex, node) in enumerate(nodelist)
			node.data.aabranchpath.paths[col] = tempaapath[nodeindex]
			node.data.aabranchpath.times[col] = tempaatimes[nodeindex]
		end
		hiddenZ = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams,false)
		prop2 = hidden_proposal_likelihood_old(nodelist, cols, col, modelparams)-hiddenZ
		hiddenAA = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, false)
		prop2 += aa_proposal_likelihood_old(nodelist, col, modelparams)-hiddenAA
	end

	if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)
		accepted_hidden += 1.0
		for (nodeindex, node) in enumerate(nodelist)				
			node.data.branchpath.paths[col] = newpath[nodeindex]
			node.data.branchpath.times[col] = newtimes[nodeindex]
			node.data.aabranchpath.paths[col] = newaapath[nodeindex]
			node.data.aabranchpath.times[col] = newaatimes[nodeindex]	
		end
	else
		for (nodeindex, node) in enumerate(nodelist)
			node.data.branchpath.paths[col] = temppath[nodeindex]
			node.data.branchpath.times[col] = temptimes[nodeindex]
			node.data.aabranchpath.paths[col] = tempaapath[nodeindex]
			node.data.aabranchpath.times[col] = tempaatimes[nodeindex]
		end
	end
	accepted_hidden_total += 1.0

	samplesiterates(rng, cols, col, nodelist, modelparams)

	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total
end

# samplepaths_simultaneous
function samplepaths_simultaneous(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false; dosamplesiterates::Bool=false, accept_everything::Bool=false)
	for node in nodelist
		joint_to_aa_h!(modelparams, node, col)
	end
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
	accepted_hidden = 0.0
	accepted_hidden_total = 0.0
	accepted_aa = 0.0
	accepted_aa_total = 0.0
	hidden_accepted = false
	aa_accepted = false	
	temppath = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
	temptimes = Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]
	tempaapath = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
	tempaatimes = Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]	
	site1 = augmentedloglikelihood(nodelist, cols, modelparams)
	cont = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams)	
	if cont
		cont = felsensteinresample_aa(rng, proteins, nodelist, cols,col, modelparams)
	end

	site2 = 0.0
	if cont
		site2 = augmentedloglikelihood(nodelist, cols, modelparams)
		newpath  = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
		newtimes =  Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]	
		newaapath  = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
		newaatimes =  Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]

		for (nodeindex, node) in enumerate(nodelist)
			node.data.branchpath.paths[col] = copy(temppath[nodeindex])
			node.data.branchpath.times[col] = copy(temptimes[nodeindex])
		end
		
		hiddenZ = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams,false)
		#println("hidden_proposal_likelihood_old")
		prop1 = hidden_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)-hiddenZ		
		#prop1 = hidden_proposal_likelihood_old(nodelist, cols, col, modelparams, temppath, temptimes)-hiddenZ		
		aaZ = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, false)		
		prop1 += aa_proposal_likelihood(nodelist, col, modelparams, tempaapath, tempaatimes)-aaZ
		#prop1 += aa_proposal_likelihood_old(nodelist, col, modelparams, tempaapath, tempaatimes)-aaZ
		for (nodeindex, node) in enumerate(nodelist)
			node.data.branchpath.paths[col] = copy(newpath[nodeindex])
			node.data.branchpath.times[col] = copy(newtimes[nodeindex])
		end

		
		for (nodeindex, node) in enumerate(nodelist)
			node.data.aabranchpath.paths[col] = copy(tempaapath[nodeindex])
			node.data.aabranchpath.times[col] = copy(tempaatimes[nodeindex])
		end
		hiddenZ = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams,false)
		#println("hidden_proposal_likelihood_new")
		prop2 = hidden_proposal_likelihood(nodelist, col, modelparams)-hiddenZ
		#prop2 = hidden_proposal_likelihood_old(nodelist, cols, col, modelparams)-hiddenZ
		#println("done")
		aaZ = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, false)
		prop2 += aa_proposal_likelihood(nodelist, col, modelparams)-aaZ
		#prop2 += aa_proposal_likelihood_old(nodelist, col, modelparams)-aaZ
	end

	#println((site2-site1)+(prop1-prop2),"\t",prop1,"\t",prop2,"\t",site1,"\t",site2)
	if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)
		accepted_hidden += 1.0
		aa_accepted = true
		hidden_accepted = true
		for (nodeindex, node) in enumerate(nodelist)
			node.data.branchpath.paths[col] = copy(newpath[nodeindex])
			node.data.branchpath.times[col] = copy(newtimes[nodeindex])
			node.data.aabranchpath.paths[col] = copy(newaapath[nodeindex])
			node.data.aabranchpath.times[col] = copy(newaatimes[nodeindex])	
		end
	else
		if accept_everything && cont
			for (nodeindex, node) in enumerate(nodelist)
				node.data.branchpath.paths[col] = copy(newpath[nodeindex])
				node.data.branchpath.times[col] = copy(newtimes[nodeindex])
				node.data.aabranchpath.paths[col] = copy(newaapath[nodeindex])
				node.data.aabranchpath.times[col] = copy(newaatimes[nodeindex])	
			end
		else
			for (nodeindex, node) in enumerate(nodelist)
				node.data.branchpath.paths[col] = copy(temppath[nodeindex])
				node.data.branchpath.times[col] = copy(temptimes[nodeindex])
				node.data.aabranchpath.paths[col] = copy(tempaapath[nodeindex])
				node.data.aabranchpath.times[col] = copy(tempaatimes[nodeindex])
			end
		end
	end
	accepted_hidden_total += 1.0

	if dosamplesiterates
		samplesiterates(rng, cols, col, nodelist, modelparams)
	end

	for node in nodelist
		aa_h_to_joint!(modelparams, node, col)
	end

	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total,hidden_accepted,aa_accepted
end

#samplepaths_seperate
function samplepaths_seperate_old(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false; dosamplesiterates::Bool=false)
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
	#=
	if useoldsampling
		temppath = Array{Int,1}[node.data.branchpath.paths[col] for node in nodelist]
		temptimes = Array{Float64,1}[node.data.branchpath.times[col] for node in nodelist]
		site1 = augmentedloglikelihood_site(nodelist, col, modelparams)	
		cont = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams, useoldsampling)
		if cont
			accepted_hidden += 1.0
		else
			index = 1
			for node in nodelist
				node.data.branchpath.paths[col] = temppath[index]
				node.data.branchpath.times[col] = temptimes[index]
				index += 1
			end
		end
		accepted_hidden_total += 1.0
	else=#
		temppath = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
		temptimes = Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]
		site1 = augmentedloglikelihood(nodelist, cols, modelparams)	
		#println("backwardsampling")
		cont = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams)
		#println("hidden_proposal_likelihood_old")
		prop1 = hidden_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
		site2 = augmentedloglikelihood(nodelist, cols, modelparams)	
		#println("hidden_proposal_likelihood_new")
		prop2 = hidden_proposal_likelihood(nodelist, col, modelparams)
		#println((site2-site1)+(prop1-prop2),"\t",prop1,"\t",prop2,"\t",site1,"\t",site2)
		#println("done")
		if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)
			accepted_hidden += 1.0
			hidden_accepted = true
		else
			for (index,node) in enumerate(nodelist)
				node.data.branchpath.paths[col] = copy(temppath[index])
				node.data.branchpath.times[col] = copy(temptimes[index])
			end
		end

		#println("A",[node.data.branchpath.paths[col] for node in nodelist],"\t",[node.data.branchpath.times[col] for node in nodelist])
		#println("B",[copy(node.data.branchpath.paths[col]) for node in nodelist], [copy(node.data.branchpath.times[col]) for node in nodelist])
		#println("C ", site1, "\t", site2)
		#println("D ", prop1, "\t", prop2)
		accepted_hidden_total += 1.0
	#end

	
	if useoldsampling
		temppath = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
		temptimes = Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]
		cont = felsensteinresample_aa(rng, proteins, nodelist, Int[col],col, modelparams, true, useoldsampling)

		if cont

			aa_accepted = true
			accepted_aa += 1.0
		else
			for (index,node) in enumerate(nodelist)
				node.data.aabranchpath.paths[col] = copy(temppath[index])
				node.data.aabranchpath.times[col] = copy(temptimes[index])
			end
		end
		accepted_aa_total += 1.0
	else
		temppath = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
		temptimes = Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]
		
		site1 = augmentedloglikelihood(nodelist, cols, modelparams)		
		cont = felsensteinresample_aa(rng, proteins, nodelist, Int[col], col, modelparams)
		site2 = augmentedloglikelihood(nodelist, cols, modelparams)

		prop1 = aa_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)		
		prop2 = aa_proposal_likelihood(nodelist, col, modelparams)
		#println(prop1,"\t", prop2, "\t")
		#println(temppath, "\t", temptimes)
		#println(Array{Float64,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist],"\t",Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist])

		if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)			
			aa_accepted = true
			accepted_aa += 1.0
		else
			for (index,node) in enumerate(nodelist)
				node.data.aabranchpath.paths[col] = copy(temppath[index])
				node.data.aabranchpath.times[col] = copy(temptimes[index])
			end
		end
		accepted_aa_total += 1.0
	end

	#samplesiterates(rng, cols, col, nodelist, modelparams)

	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total, hidden_accepted, aa_accepted
end

#=
function sample_dataset(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false; dosamplesiterates::Bool=false, accept_everything::Bool=false)
	if lengths(proteins) == 1
		samplesingle(nodelist[1], modelparams)
	else
		samplepaths_seperate_new(rng, col, proteins, nodelist, modelparams, useoldsampling, dosamplesiterates=dosamplesiterates, accept_everything=accept_everything)
	end
end=#

function samplepaths_seperate_new(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false; dosamplesiterates::Bool=false, accept_everything::Bool=false)
	for node in nodelist
		joint_to_aa_h!(modelparams, node, col)
	end
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
	if modelparams.numhiddenstates == 1
		accepted_hidden += 1.0
		hidden_accepted = true
	else
		temppath = Array{Int,1}[copy(node.data.branchpath.paths[col]) for node in nodelist]
		temptimes = Array{Float64,1}[copy(node.data.branchpath.times[col]) for node in nodelist]
		site1 = augmentedloglikelihood(nodelist, cols, modelparams)	
		cont = felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams)
		prop1 = hidden_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
		site2 = augmentedloglikelihood(nodelist, cols, modelparams)	
		prop2 = hidden_proposal_likelihood(nodelist, col, modelparams)
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
		end
	end
	accepted_hidden_total += 1.0

	cols = Int[col]
	temppath = Array{Int,1}[copy(node.data.aabranchpath.paths[col]) for node in nodelist]
	temptimes = Array{Float64,1}[copy(node.data.aabranchpath.times[col]) for node in nodelist]	
	site1 = augmentedloglikelihood(nodelist, cols, modelparams)
	cont = felsensteinresample_aa(rng, proteins, nodelist, cols, col, modelparams)
	site2 = augmentedloglikelihood(nodelist, cols, modelparams)
	prop1 = aa_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
	prop2 = aa_proposal_likelihood(nodelist, col, modelparams)
	if cont && exp((site2-site1)+(prop1-prop2)) > rand(rng)			
		aa_accepted = true
		accepted_aa += 1.0
	else
		if accept_everything && cont

		else
			for (index,node) in enumerate(nodelist)
				node.data.aabranchpath.paths[col] = copy(temppath[index])
				node.data.aabranchpath.times[col] = copy(temptimes[index])
			end
		end
	end
	accepted_aa_total += 1.0

	if dosamplesiterates
		samplesiterates(rng, cols, col, nodelist, modelparams)
	end

	for node in nodelist
		aa_h_to_joint!(modelparams, node, col)
	end

	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total, hidden_accepted, aa_accepted
end

# samplepaths_simultaneous
function samplepaths(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false; dosamplesiterates::Bool=false, accept_everything::Bool=false)
	for node in nodelist
		aa_h_to_joint!(modelparams, node, col)
	end

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
	accepted_hidden = 0.0
	accepted_hidden_total = 0.0
	accepted_aa = 0.0
	accepted_aa_total = 0.0
	hidden_accepted = false
	aa_accepted = false	
	temppath = Array{Int,1}[copy(node.data.jointbranchpath.paths[col]) for node in nodelist]
	temptimes = Array{Float64,1}[copy(node.data.jointbranchpath.times[col]) for node in nodelist]
	tempproposalliks = Float64[node.data.proposallikelihood[col] for node in nodelist]
	site1 = augmentedloglikelihood(nodelist, cols, modelparams)
	cont = felsensteinresample_joint(rng, proteins, nodelist, col, cols,col, modelparams)	
	site2 = 0.0
	prop1 = 0.0
	prop2 = 0.0
	if cont
		prop2_n = joint_proposal_likelihood(nodelist, col, modelparams)

		site2 = augmentedloglikelihood(nodelist, cols, modelparams)
		newpath  = Array{Int,1}[copy(node.data.jointbranchpath.paths[col]) for node in nodelist]
		newtimes =  Array{Float64,1}[copy(node.data.jointbranchpath.times[col]) for node in nodelist]

		for (nodeindex, node) in enumerate(nodelist)
			node.data.jointbranchpath.paths[col] = copy(temppath[nodeindex])
			node.data.jointbranchpath.times[col] = copy(temptimes[nodeindex])
		end
		
		Z1 = felsensteinresample_joint(rng, proteins, nodelist, col, cols,col, modelparams,false)
		#prop1 = joint_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
		prop1_n = joint_proposal_likelihood(nodelist, col, modelparams)	
		for (nodeindex, node) in enumerate(nodelist)
			node.data.jointbranchpath.paths[col] = copy(newpath[nodeindex])
			node.data.jointbranchpath.times[col] = copy(newtimes[nodeindex])
		end
		#Z2 = felsensteinresample_joint(rng, proteins, nodelist, col, cols,col, modelparams,false)
		#prop2 = joint_proposal_likelihood(nodelist, col, modelparams)
		#=
		if rand(rng) < 0.01
			println(prop1,"\t",prop1_n,"\t",prop2,"\t",prop2_n)
		end
		=#
	end 

	if cont && exp((site2-site1)+(prop1_n-prop2_n)) > rand(rng)
		accepted_hidden += 1.0
		aa_accepted = true
		hidden_accepted = true
		for (nodeindex, node) in enumerate(nodelist)
			node.data.jointbranchpath.paths[col] = copy(newpath[nodeindex])
			node.data.jointbranchpath.times[col] = copy(newtimes[nodeindex])
		end
	else
		if cont && accept_everything

		else
			for (nodeindex, node) in enumerate(nodelist)
				node.data.jointbranchpath.paths[col] = copy(temppath[nodeindex])
				node.data.jointbranchpath.times[col] = copy(temptimes[nodeindex])
			end
		end
	end
	accepted_hidden_total += 1.0

	for node in nodelist
		joint_to_aa_h!(modelparams, node, col)
	end

	if dosamplesiterates
		samplesiterates(rng, cols, col, nodelist, modelparams)
	end
	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total,hidden_accepted,aa_accepted
end

function samplepaths_newnew(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams, useoldsampling::Bool=false; dosamplesiterates::Bool=false, accept_everything::Bool=false)
	for node in nodelist
		aa_h_to_joint!(modelparams, node, col)
	end

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
	accepted_hidden = 0.0
	accepted_hidden_total = 0.0
	accepted_aa = 0.0
	accepted_aa_total = 0.0
	hidden_accepted = false
	aa_accepted = false	
	temppath = Array{Int,1}[copy(node.data.jointbranchpath.paths[col]) for node in nodelist]
	temptimes = Array{Float64,1}[copy(node.data.jointbranchpath.times[col]) for node in nodelist]
	site1 = augmentedloglikelihood(nodelist, cols, modelparams)
	cont = felsensteinresample_joint(rng, proteins, nodelist, col, cols,col, modelparams)
	site2 = 0.0
	prop1 = 0.0
	prop2 = 0.0
	if cont
		site2 = augmentedloglikelihood(nodelist, cols, modelparams)
		newpath  = Array{Int,1}[copy(node.data.jointbranchpath.paths[col]) for node in nodelist]
		newtimes =  Array{Float64,1}[copy(node.data.jointbranchpath.times[col]) for node in nodelist]

		for (nodeindex, node) in enumerate(nodelist)
			node.data.jointbranchpath.paths[col] = copy(temppath[nodeindex])
			node.data.jointbranchpath.times[col] = copy(temptimes[nodeindex])
		end
		
		Z1 = felsensteinresample_joint(rng, proteins, nodelist, col, cols,col, modelparams,false)
		prop1 = sum(Float64[node.data.proposallikelihood[col] for node in nodelist])
		for (nodeindex, node) in enumerate(nodelist)
			node.data.branchpath.paths[col] = copy(newpath[nodeindex])
			node.data.branchpath.times[col] = copy(newtimes[nodeindex])
		end
		Z2 = felsensteinresample_joint(rng, proteins, nodelist, col, cols,col, modelparams,false)
		prop2 = sum(Float64[node.data.proposallikelihood[col] for node in nodelist])
	end 

	delta = (site2-site1)+(prop1-prop2)
	#println(cont,"\t",delta)
	if cont && exp(delta) > rand(rng)
		accepted_hidden += 1.0
		aa_accepted = true
		hidden_accepted = true
		for (nodeindex, node) in enumerate(nodelist)
			node.data.jointbranchpath.paths[col] = copy(newpath[nodeindex])
			node.data.jointbranchpath.times[col] = copy(newtimes[nodeindex])
		end
	else
		for (nodeindex, node) in enumerate(nodelist)
			node.data.jointbranchpath.paths[col] = copy(temppath[nodeindex])
			node.data.jointbranchpath.times[col] = copy(temptimes[nodeindex])
		end
	end
	accepted_hidden_total += 1.0

	for node in nodelist
		joint_to_aa_h!(modelparams, node, col)
	end

	if dosamplesiterates
		samplesiterates(rng, cols, col, nodelist, modelparams)
	end
	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total,hidden_accepted,aa_accepted
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