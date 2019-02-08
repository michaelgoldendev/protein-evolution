using Distributions
using LinearAlgebra
using Serialization
using JSON

push!(LOAD_PATH,string(@__DIR__,"/../../MolecularEvolution/src/"))
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
using Formatting
using NLopt
using SpecialFunctions

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

#=
function root(node::TreeNode)
	parents = []
	current = node
	while !isnull(current.parent)
		current = get(current.parent)
		push!(parents, current)
	end
	println([p.nodeindex for p in parents])

	bottom = parents[end]
	while length(parents) > 0
		pop!(parents)
		println("A",[b.nodeindex for b in bottom])
		println("B",bottom.nodeindex,"\t",parents[end].nodeindex)		
		exit()
		leftchild = bottom.children[1]
		push!(leftchild.children,bottom)		
		bottom.children = bottom.children[2:end]


		bottom.parent = leftchild
		newroot = TreeNode(0.1, "dummy$(length(parents))")
		push!(newroot.children,leftchild.children[1])
		push!(newroot.children,leftchild)
		leftchild.children[1].parent = Nullable{TreeNode}(newroot)
		leftchild.parent = Nullable{TreeNode}(newroot)
		leftchild.children = leftchild.children[2:end]
		bottom = newroot
	end
	return bottom
end=#

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
	#root(newroot)

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

mutable struct ModelParams
    alphabet::Int
    aminoacidQ::Array{Float64,2}
    aa_exchangeablities::Array{Float64,2}
    numrates::Int    
    rates::Array{Float64,1}
    rate_exchangeablities::Array{Float64,2}
    rate_freqs::Array{Float64,1}
    rate_mu::Float64
    rate_alpha::Float64
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
			push!(hiddennodes, HiddenNode())
		end
		matrixcache = Dict{Tuple{Int,Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}()
		numrates = 3
		rate_mu = 50.0
		rate_alpha = 5.0
		rates = discretizegamma(rate_alpha, 1.0/rate_alpha, numrates)
		rate_exchangeablities = ones(Float64, numrates, numrates)*50.0
		rate_freqs = ones(Float64, numrates)/numrates
        new(20,aminoacidQ,LGexchangeability,numrates, rates,rate_exchangeablities, rate_freqs, rate_mu, rate_alpha,numhiddenstates,initialprobs,transitioncounts,transitionprobs,transitionrates,hiddennodes,mu,hiddenmu,matrixcache,1.0,1.0)
    end
end

function estimate_hidden_transition_probs(modelparams::ModelParams)
	for h=1:modelparams.numhiddenstates
		modelparams.transitionprobs[h,:] = modelparams.transitioncounts[h,:] ./ sum(modelparams.transitioncounts[h,:])
	end
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
	modelparams.matrixcache = Dict{Tuple{Int,Int}, Tuple{Array{Float64,2},Array{Complex{Float64},2},Array{Complex{Float64},1},Array{Complex{Float64},2}}}()
end

function getRandPt(modelparams::ModelParams, prev_hmm::Int, next_hmm::Int, h::Int, aa::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = (-2,-2,h,aa)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1
		Q = constructRateMatrix(modelparams, prev_hmm, next_hmm,h, aa)
		
		
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

function getAAandPt(modelparams::ModelParams, prev_hmm::Int, next_hmm::Int, h::Int, ratecat::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = (-1,-1,h, ratecat)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1
		Q = constructAAMatrix(modelparams, prev_hmm, next_hmm, h, ratecat)
		
		
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


function getQandPt(modelparams::ModelParams, prevh::Int, nexth::Int, aa::Int, ratecat::Int, t::Float64)
	global countcachemisses
	global countcachehits

	key = (prevh, nexth, aa, ratecat)
	if !haskey(modelparams.matrixcache, key)		
		countcachemisses += 1

		Q = constructHiddenMatrix(modelparams, prevh, nexth, aa, ratecat)
		
		
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
	jointbranchpath::BranchPath
	ratesbranchpath::BranchPath
	dummy::Int
	protein::Protein
	inputbranchlength::Float64

	AugmentedNodeData(col::Int) = new(BranchPath(col),BranchPath(col),BranchPath(col),BranchPath(col),1, Protein(),1.0) 
	AugmentedNodeData(branchpath::BranchPath, aabranchpath::BranchPath, jointbranchpath::BranchPath, ratesbranchpath::BranchPath, dummy::Int) = new(branchpath, aabranchpath, jointbranchpath, ratesbranchpath, dummy, Protein(),1.0)
end

function felsensteinhelper(node::TreeNode, selcolin::Int, incols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, incols)
	cols = copy(incols)
	deleteat!(cols,selcol)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol]), BranchPathIterator(node.data.ratesbranchpath,Int[aacol])])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	ratesiter = multi_iter.branchpathiterators[3]
	P = Matrix{Float64}(I, modelparams.numhiddenstates, modelparams.numhiddenstates)
	Rmatrices = Array{Float64,2}[]
	Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	index = 1
    for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		R,Pi = getQandPt(modelparams, get(hiddeniter.prevstates, selcol-1, 0), get(hiddeniter.prevstates, selcol, 0), aaiter.prevstates[1], ratesiter.prevstates[1], dt)
		#R,Pi = getQandPt(modelparams, prevh, succh, aaiter.currstates[1], dt)
		P *= Pi
		push!(Rmatrices, R*dt)
    	push!(Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)
    	if index == 1
    		node.data.branchpath.R = R
    	end
    	index += 1
	end
	push!(dummytime,1.0)

    tempv = copy(v)
    pushfirst!(vs,v)
    for P in reverse(Pmatrices)
        tempv = P*tempv
    	pushfirst!(vs,copy(tempv))
    end
    popfirst!(vs)

    node.data.branchpath.Rmatrices = Rmatrices
    node.data.branchpath.Pmatrices = Pmatrices
    node.data.branchpath.vs = vs
    node.data.branchpath.time = dummytime
    node.data.branchpath.P = P
    return P,Pmatrices,vs
end




function felsensteinresample(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, selcolin::Int, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
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
	rootstate = CommonUtils.sample(rng,rootliks)
	rootnode.data.branchpath.paths[selcol] = Int[rootstate]
	rootnode.data.branchpath.times[selcol] = Float64[0.0]
	ll = backwardsampling(rng,nodelist[1], rootstate, selcol, modelparams)
	return log(freqs[rootstate])+ll
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
	h = rootnode.data.branchpath.paths[aacol][end]
	loglikelihood = log(freqs[h])
	nodeindex = 1
	for node in nodelist
		if !isroot(node)
			P = node.data.branchpath.P
			R = node.data.branchpath.R
			path = node.data.branchpath.paths[aacol]
			time = node.data.branchpath.times[aacol]
			if length(paths) > 0
				path = paths[nodeindex]
				time = times[nodeindex]
			end
			loglikelihood += log(P[path[1],path[end]])
			
			for i=1:length(path)-1
				dt = (time[i+1]-time[i])*node.branchlength
				loglikelihood += R[path[i],path[i]]*dt
				loglikelihood += R[path[i],path[i+1]]*node.branchlength
			end
			dt = (1.0 - time[end])*node.branchlength
			loglikelihood += R[path[end],path[end]]*dt
		end
		#=		h = node.data.branchpath.paths[col][end]
		if node.seqindex > 0 && aacol < length(proteins[node.seqindex].sites)
			loglikelihood += siteloglikelihood(proteins[node.seqindex].sites[aacol], h, modelparams)
		end=#
		nodeindex += 1
	end
	return loglikelihood
end

function backwardsampling(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		P = child.data.branchpath.P
		R = child.data.branchpath.R
		liks = P[state,:].*child.data.branchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime = modifiedrejectionsampling(rng, R*child.branchlength, state, b,(modelparams))
		child.data.branchpath.paths[aacol] = newpath
		child.data.branchpath.times[aacol] = newtime
		backwardsampling(rng,child, b,aacol, modelparams)
	end
	return 0.0
end

function felsensteinhelper_aa(node::TreeNode, cols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	cols = Int[aacol]
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.ratesbranchpath,cols)])
	hiddeniter = multi_iter.branchpathiterators[1]
	ratesiter = multi_iter.branchpathiterators[2]
	Rmatrices = Array{Float64,2}[]
	Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	P = Matrix{Float64}(I, modelparams.alphabet, modelparams.alphabet)
	index = 1
    for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		R,Pi = getAAandPt(modelparams, 0, 0, hiddeniter.prevstates[1], ratesiter.prevstates[1], dt)				
    	P *= Pi
		push!(Rmatrices, R*dt)
    	push!(Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)
    	if index == 1
    		node.data.aabranchpath.R = R
    	end
    	index += 1
	end
	push!(dummytime,1.0)

    tempv = copy(v)
    pushfirst!(vs,v)
    for P in reverse(Pmatrices)    	
        tempv = P*tempv
    	pushfirst!(vs,copy(tempv))
    end
    popfirst!(vs)

    node.data.aabranchpath.Rmatrices = Rmatrices
    node.data.aabranchpath.Pmatrices = Pmatrices
    node.data.aabranchpath.vs = vs
    node.data.aabranchpath.time = dummytime
    node.data.aabranchpath.P = P
    return P,Pmatrices,vs
end

function felsensteinresample_aa(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
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
				for a=1:modelparams.alphabet
					likelihoods[nodeindex,a] = 1.0
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
	rootstate = CommonUtils.sample(rng,rootliks)
	rootnode.data.aabranchpath.paths[aacol] = Int[rootstate]
	rootnode.data.aabranchpath.times[aacol] = Float64[0.0]
	#backwardsampling_aa(rng,nodelist[1], rootstate, selcol,likelihoods,print,modelparams)
	ll = backwardsampling_aa(rng,nodelist[1], rootstate, aacol,modelparams)
	#return log(modelparams.hiddennodes[h].aa_node.probs[rootstate]) + ll
	return log(modelparams.hiddennodes[h].aa_node.probs[rootstate]) + ll
end


function aa_proposal_likelihood(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams, paths::Array{Array{Int,1},1}=Array{Int,1}[], times::Array{Array{Float64,1},1}=Array{Float64,1}[])
	h = nodelist[1].data.branchpath.paths[aacol][end]
	loglikelihood = log(modelparams.hiddennodes[h].aa_node.probs[nodelist[1].data.aabranchpath.paths[aacol][end]])
	nodeindex = 1
	for node in nodelist
		if !isroot(node)
			P = node.data.aabranchpath.P
			R = node.data.aabranchpath.R
			path = node.data.aabranchpath.paths[aacol]
			time = node.data.aabranchpath.times[aacol]
			if length(paths) > 0
				path = paths[nodeindex]
				time = times[nodeindex]
			end
			loglikelihood += log(P[path[1],path[end]])
			
			for i=1:length(path)-1
				dt = (time[i+1]-time[i])*node.branchlength
				loglikelihood += R[path[i],path[i]]*dt
				loglikelihood += R[path[i],path[i+1]]*node.branchlength
			end
			dt = (1.0 - time[end])*node.branchlength
			loglikelihood += R[path[end],path[end]]*dt
		end
		#=		h = node.data.branchpath.paths[col][end]
		if node.seqindex > 0 && aacol < length(proteins[node.seqindex].sites)
			loglikelihood += siteloglikelihood(proteins[node.seqindex].sites[aacol], h, modelparams)
		end=#
		nodeindex += 1
	end
	return loglikelihood
end

function backwardsampling_aa(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		P = child.data.aabranchpath.P
		R = child.data.aabranchpath.R
		liks = P[state,:].*child.data.aabranchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime = modifiedrejectionsampling(rng, R*child.branchlength, state, b,(modelparams))
		child.data.aabranchpath.paths[aacol] = newpath
		child.data.aabranchpath.times[aacol] = newtime
		backwardsampling_aa(rng,child, b,aacol, modelparams)
	end
	return 0.0
end

function felsensteinhelper_rates(node::TreeNode, cols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	cols = Int[aacol]
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	#ratesiter = multi_iter.branchpathiterators[3]
	Rmatrices = Array{Float64,2}[]
	Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	P = Matrix{Float64}(I, modelparams.numrates, modelparams.numrates)
	index = 1
    for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength

		R,Pi =  getRandPt(modelparams, 0, 0, hiddeniter.prevstates[1], aaiter.prevstates[1], dt)		
    	P *= Pi
		push!(Rmatrices, R*dt)
    	push!(Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)
    	if index == 1
    		node.data.ratesbranchpath.R = R
    	end
    	index += 1
	end
	push!(dummytime,1.0)

    tempv = copy(v)
    pushfirst!(vs,v)
    for P in reverse(Pmatrices)    	
        tempv = P*tempv
    	pushfirst!(vs,copy(tempv))
    end
    popfirst!(vs)

    node.data.ratesbranchpath.Rmatrices = Rmatrices
    node.data.ratesbranchpath.Pmatrices = Pmatrices
    node.data.ratesbranchpath.vs = vs
    node.data.ratesbranchpath.time = dummytime
    node.data.ratesbranchpath.P = P
    return P,Pmatrices,vs
end

function felsensteinresample_rates(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
	likelihoods = ones(Float64, length(nodelist), modelparams.numrates)*-Inf
	logm = zeros(Float64,length(nodelist))

	stack = Int[1]
	while length(stack) > 0
		nodeindex = stack[end]
		node = nodelist[nodeindex]
		if isleafnode(node)
			v = rates_helper(node, aacol, modelparams)
			for a=1:modelparams.numrates
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
				lefttransprobs,Pmatrices_left,vs_left = felsensteinhelper_rates(nodelist[leftchildindex], cols, aacol, likelihoods[leftchildindex,:], modelparams)
				righttransprobs,Pmatrices_right,vs_right = felsensteinhelper_rates(nodelist[rightchildindex], cols, aacol, likelihoods[rightchildindex,:], modelparams)
        		likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])

        		m = maximum(likelihoods[nodeindex,:])
				likelihoods[nodeindex,:] = likelihoods[nodeindex,:] ./ m
				logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
        		pop!(stack)
        	end
        end
    end

	rootnode = nodelist[1]
	rootstate = CommonUtils.sample(rng,likelihoods[1,:].*modelparams.rate_freqs)
	rootnode.data.ratesbranchpath.paths[aacol] = Int[rootstate]
	rootnode.data.ratesbranchpath.times[aacol] = Float64[0.0]
	backwardsampling_rates(rng,nodelist[1], rootstate, aacol,modelparams)
end

function backwardsampling_rates(rng::AbstractRNG, node::TreeNode, state::Int, aacol::Int, modelparams::ModelParams)
	for child in node
		P = child.data.ratesbranchpath.P
		R = child.data.ratesbranchpath.R
		liks = P[state,:].*child.data.ratesbranchpath.vs[end]
		b = CommonUtils.sample(rng,liks)
		newpath, newtime = modifiedrejectionsampling(rng, R*child.branchlength, state, b,(modelparams))
		child.data.ratesbranchpath.paths[aacol] = newpath
		child.data.ratesbranchpath.times[aacol] = newtime
		backwardsampling_rates(rng,child, b,aacol, modelparams)
	end
	return 0.0
end

function rates_proposal_likelihood(nodelist::Array{TreeNode,1}, aacol::Int, modelparams::ModelParams, paths::Array{Array{Int,1},1}=Array{Int,1}[], times::Array{Array{Float64,1},1}=Array{Float64,1}[])
	loglikelihood = log(modelparams.rate_freqs[nodelist[1].data.ratesbranchpath.paths[aacol][end]])
	nodeindex = 1
	for node in nodelist
		if !isroot(node)
			P = node.data.ratesbranchpath.P
			R = node.data.ratesbranchpath.R
			path = node.data.ratesbranchpath.paths[aacol]
			time = node.data.ratesbranchpath.times[aacol]
			if length(paths) > 0
				path = paths[nodeindex]
				time = times[nodeindex]
			end
			loglikelihood += log(P[path[1],path[end]])
			
			for i=1:length(path)-1
				dt = (time[i+1]-time[i])*node.branchlength
				loglikelihood += R[path[i],path[i]]*dt
				loglikelihood += R[path[i],path[i+1]]*node.branchlength
			end
			dt = (1.0 - time[end])*node.branchlength
			loglikelihood += R[path[end],path[end]]*dt
			loglikelihood += log(rates_helper(node, aacol, modelparams)[path[end]])
		end
		#=		h = node.data.branchpath.paths[col][end]
		if node.seqindex > 0 && aacol < length(proteins[node.seqindex].sites)
			loglikelihood += siteloglikelihood(proteins[node.seqindex].sites[aacol], h, modelparams)
		end=#
		nodeindex += 1
	end
	return loglikelihood
end

function rates_helper(node::TreeNode, col::Int, modelparams::ModelParams)
	numcols = length(node.data.branchpath.paths)

	v = zeros(Float64, modelparams.numrates)
	for rateindex=1:modelparams.numrates
		loglikelihood = 0.0
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
		#ratesiter = multi_iter.branchpathiterators[3]

		for it in multi_iter
			dt = (multi_iter.currtime-multi_iter.prevtime)
			changecol = selcol
			Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol], rateindex, rateindex)
			loglikelihood += Qii*dt*node.branchlength

			changecol = 0
			if multi_iter.branchpathindex == 1
				changecol = hiddeniter.mincol
			elseif multi_iter.branchpathindex == 2
				changecol = aaiter.mincol
			end

			if changecol == selcol
				prevstatesh = hiddeniter.prevstates[changecol]
				currstatesh = prevstatesh
				prevstatesaa = aaiter.prevstates[changecol]
				currstatesaa = prevstatesaa
				if multi_iter.branchpathindex == 1
					currstatesh = hiddeniter.currstates[changecol]
				elseif multi_iter.branchpathindex == 2
					currstatesaa = aaiter.currstates[changecol]
				end

				Qhi = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), prevstatesh, currstatesh, prevstatesaa, currstatesaa, rateindex, rateindex)
				loglikelihood += log(Qhi*node.branchlength)
			end
		end
		v[rateindex] = loglikelihood
	end

	
	v = v .- maximum(v)
	v = exp.(v)
	v = v/sum(v)
	return v
end

function augmentedloglikelihood_site(nodelist::Array{TreeNode,1}, col::Int, modelparams::ModelParams)
	numcols = length(nodelist[1].data.branchpath.paths)
	loglikelihood =  log(modelparams.rate_freqs[nodelist[1].data.ratesbranchpath.paths[col][end]])
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

			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols), BranchPathIterator(node.data.ratesbranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]
			ratesiter = multi_iter.branchpathiterators[3]

			for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)
				changecol = selcol
				Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol], ratesiter.prevstates[changecol], ratesiter.prevstates[changecol])
				loglikelihood += Qii*dt*node.branchlength

				changecol = 0
				if multi_iter.branchpathindex == 1
					changecol = hiddeniter.mincol
				elseif multi_iter.branchpathindex == 2
					changecol = aaiter.mincol
				else
					changecol = ratesiter.mincol
				end

				if changecol == selcol
					prevstatesh = hiddeniter.prevstates[changecol]
					currstatesh = prevstatesh
					prevstatesaa = aaiter.prevstates[changecol]
					currstatesaa = prevstatesaa
					prevstatesrates = ratesiter.prevstates[changecol]
					currstatesrates = prevstatesrates
					if multi_iter.branchpathindex == 1
						currstatesh = hiddeniter.currstates[changecol]
					elseif multi_iter.branchpathindex == 2
						currstatesaa = aaiter.currstates[changecol]
					elseif multi_iter.branchpathindex == 3
						currstatesrates = ratesiter.currstates[changecol]
					end

					Qhi = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), prevstatesh, currstatesh, prevstatesaa, currstatesaa, prevstatesrates, currstatesrates)
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
		loglikelihood += log(modelparams.rate_freqs[nodelist[1].data.ratesbranchpath.paths[col][end]])
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

				multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols), BranchPathIterator(node.data.ratesbranchpath,cols)])
				hiddeniter = multi_iter.branchpathiterators[1]
				aaiter = multi_iter.branchpathiterators[2]
				ratesiter = multi_iter.branchpathiterators[3]

				for it in multi_iter
					dt = (multi_iter.currtime-multi_iter.prevtime)
					changecol = selcol
					Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol], ratesiter.prevstates[changecol], ratesiter.prevstates[changecol])
					loglikelihood += Qii*dt*node.branchlength

					changecol = 0
					if multi_iter.branchpathindex == 1
						changecol = hiddeniter.mincol
					elseif multi_iter.branchpathindex == 2
						changecol = aaiter.mincol
					else
						changecol = ratesiter.mincol
					end

					if changecol == selcol
						prevstatesh = hiddeniter.prevstates[changecol]
						currstatesh = prevstatesh
						prevstatesaa = aaiter.prevstates[changecol]
						currstatesaa = prevstatesaa
						prevstatesrates = ratesiter.prevstates[changecol]
						currstatesrates = prevstatesrates
						if multi_iter.branchpathindex == 1
							currstatesh = hiddeniter.currstates[changecol]
						elseif multi_iter.branchpathindex == 2
							currstatesaa = aaiter.currstates[changecol]
						elseif multi_iter.branchpathindex == 3
							currstatesrates = ratesiter.currstates[changecol]
						end

						Qhi = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), prevstatesh, currstatesh, prevstatesaa, currstatesaa, prevstatesrates, currstatesrates)
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

function getrateentry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, h::Int, aa::Int, prevrate::Int, currrate::Int)
	if prevrate != currrate
		return modelparams.scalingfactor*modelparams.rate_exchangeablities[prevrate,currrate]*modelparams.rate_freqs[currrate]*modelparams.hiddennodes[h].aa_node.probs[aa]
	elseif prevrate == currrate
		q = 0.0
		for r=1:modelparams.numrates
			if r != prevrate
				q -= getrateentry(modelparams,prevh_hmm,nexth_hmm, h, aa, prevrate, r)
			end
		end
		return q
	end
	#=
	if abs(prevrate - currrate) == 1
		#return modelparams.rate_mu*modelparams.hiddennodes[h].aa_node.probs[aa]
		modelparams.rate_exchangeablities[prevrate,currrate]*modelparams.hiddennodes[h].aa_node.probs[aa]
	elseif prevrate != currrate
		return 0.0
	elseif prevrate == currrate
		q = 0.0
		for r=1:modelparams.numrates
			if r != prevrate
				q -= getrateentry(modelparams,prevh_hmm,nexth_hmm, h, aa, prevrate, r)
			end
		end
		return q
	end=#
end

function getaaentry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, h::Int, prevaa::Int, curraa::Int, ratecat::Int)
	if prevaa != curraa
		#=
		prevprob = 1.0
		if prevh_hmm != 0
			prevprob = modelparams.transitionprobs[prevh_hmm,h]
		end
		nextprob = 1.0
		if nexth_hmm != 0
			nextprob = modelparams.transitionprobs[h,nexth_hmm]
		end
		return modelparams.mu*prevprob*nextprob*modelparams.aa_exchangeablities[prevaa,curraa]*modelparams.hiddennodes[h].aa_node.probs[curraa]
		=#
		return modelparams.scalingfactor*modelparams.rates[ratecat]*modelparams.mu*modelparams.aa_exchangeablities[prevaa,curraa]*modelparams.hiddennodes[h].aa_node.probs[curraa]
	else
		q = 0.0
		for aa=1:modelparams.alphabet
			if aa != prevaa
				q -= getaaentry(modelparams,prevh_hmm,nexth_hmm, h, prevaa, aa, ratecat)
			end
		end
		return q
	end
end

function gethiddenentry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, prevh::Int, currh::Int, aa::Int, ratecat::Int)
	if prevh != currh
		#=
		prevprob = 1.0
		if prevh_hmm != 0
			prevprob = modelparams.transitionprobs[prevh_hmm,currh]
		end
		nextprob = 1.0
		if nexth_hmm != 0
			nextprob = modelparams.transitionprobs[currh,nexth_hmm]
		end
		return modelparams.hiddenmu*prevprob*nextprob*modelparams.transitionrates[prevh,currh]*modelparams.hiddennodes[currh].aa_node.probs[aa]
		=#
		
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
		return modelparams.scalingfactor*modelparams.rates[ratecat]*sqrt((prevprob*nextprob)/(prevprob2*nextprob2))*modelparams.hiddenmu*modelparams.transitionrates[prevh,currh]*modelparams.hiddennodes[currh].aa_node.probs[aa]
		
	else
		q = 0.0
		for h=1:modelparams.numhiddenstates
			if h != prevh
				q -= gethiddenentry(modelparams, prevh_hmm, nexth_hmm, prevh, h, aa, ratecat)
			end
		end
		return q
	end
end

function entry(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, prevh::Int, currh::Int, prevaa::Int, curraa::Int, prevrate::Int, currrate::Int)
	if prevh == currh && prevaa == curraa && prevrate == currrate
		q = 0.0
		for h=1:modelparams.numhiddenstates
			for aa=1:modelparams.alphabet
				for r=1:modelparams.numrates
					if ((h != prevh) + (aa != prevaa) + (r != prevrate)) == 1
						q -= entry(modelparams, prevh_hmm, nexth_hmm, prevh, h, prevaa, aa, prevrate, r)
					end
				end
			end
		end
		return q
	elseif prevh == currh && prevrate == currrate
		return getaaentry(modelparams, prevh_hmm, nexth_hmm, prevh, prevaa, curraa, prevrate)
	elseif prevaa == curraa && prevrate == currrate
		return gethiddenentry(modelparams, prevh_hmm, nexth_hmm, prevh, currh, prevaa, prevrate)
	elseif prevaa == curraa && prevh == currh
		return getrateentry(modelparams, prevh_hmm, nexth_hmm, prevh, prevaa, prevrate, currrate)
	else
		return 0.0
	end
end

function constructRateMatrix(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, h::Int, aa::Int)
	Q = zeros(Float64, modelparams.numrates, modelparams.numrates)
	for r1=1:modelparams.numrates
        for r2=1:modelparams.numrates
			if r1 != r2				
				Q[r1,r2] = getrateentry(modelparams, prevh_hmm, nexth_hmm, h, aa, r1, r2)
				Q[r1,r1] -= Q[r1,r2]
			end
		end
	end

    return Q
end

function constructAAMatrix(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, h::Int, ratecat::Int)
	Q = zeros(Float64, modelparams.alphabet, modelparams.alphabet)
	for aa1=1:modelparams.alphabet
        for aa2=1:modelparams.alphabet
			if aa1 != aa2
				Q[aa1,aa2] = getaaentry(modelparams, prevh_hmm, nexth_hmm, h, aa1, aa2, ratecat)
				Q[aa1,aa1] -= Q[aa1,aa2]
			end
		end
	end

    return Q
end

function constructHiddenMatrix(modelparams::ModelParams, prevh_hmm::Int, nexth_hmm::Int, aa::Int, ratecat::Int)
	Q = zeros(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)
	for h1=1:modelparams.numhiddenstates
        for h2=1:modelparams.numhiddenstates
			if h1 != h2
				Q[h1,h2] = gethiddenentry(modelparams, prevh_hmm, nexth_hmm, h1, h2, aa, ratecat)
				Q[h1,h1] -= Q[h1,h2]
			end
		end
	end

    return Q
end

function siteloglikelihood(site::SiteObservation, h::Int, modelparams::ModelParams)
	ll = 0.0
	if site.aa > 0
		#ll += log(modelparams.hiddennodes[h].aa_node.probs[site.aa])
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
	    v = zeros(Float64, modelparams.numhiddenstates)
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
			#site.aa = indexof(string(aa), aminoacids)
			site.h =  node.data.branchpath.paths[col][end]
			h = site.h
			hiddennode = modelparams.hiddennodes[h]
			#site.aa = CommonUtils.sample(rng, hiddennode.aa_node.probs)
			site.aa = node.data.aabranchpath.paths[col][end]
			#=
			site.phi = pimod(vonmisesrand(rng, hiddennode.phi_node.dist))
			site.omega = pimod(vonmisesrand(rng, hiddennode.omega_node.dist))
			site.psi = pimod(vonmisesrand(rng, hiddennode.psi_node.dist))
			=#
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
	node_rates = zeros(Int,length(nodelist),numcols)
	#randomratecat1 = rand(1:modelparams.numrates)
	#randomratecat2 = max(1,min(modelparams.numrates, randomratecat1 + rand(1:3) - 2))
	randomratecat1  = div(modelparams.numrates,2)
	randomratecat2  = randomratecat1
	for node in reversenodelist
		for col=1:numcols
			if isleafnode(node)
				node_aa[node.nodeindex,col] = rand(rng,1:modelparams.alphabet)
				if 1 <= col <= length(node.data.protein.sites) && node.data.protein.sites[col].aa > 0
					node_aa[node.nodeindex,col] = node.data.protein.sites[col].aa
				end

				v = zeros(Float64, modelparams.numhiddenstates)
				for h=1:modelparams.numhiddenstates
					v[h] = modelparams.hiddennodes[h].aa_node.probs[node_aa[node.nodeindex,col]]
				end
				node_h[node.nodeindex,col] = CommonUtils.sample(rng, v)
				if rand(rng) < 0.5
					node_rates[node.nodeindex,col] = randomratecat1
				else
					node_rates[node.nodeindex,col] = randomratecat2
				end
			else
				randchildindex = node.children[rand(1:length(node.children))].nodeindex
				node_aa[node.nodeindex,col] = node_aa[randchildindex, col]
				#randchildindex = node.children[rand(1:length(node.children))].nodeindex
				node_h[node.nodeindex,col] = node_h[randchildindex, col]
				node_rates[node.nodeindex,col] = node_rates[randchildindex,col]
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
	ratepaths = Array{Int,1}[]
	ratetimes = Array{Float64,1}[]
	for col=1:numcols
		push!(ratepaths,Int[randomratecat1,randomratecat1])
		push!(ratetimes,Float64[0.0, 1.0])
	end	
	root.data.branchpath = BranchPath(paths,times)
	root.data.aabranchpath = BranchPath(aapaths,aatimes)
	root.data.ratesbranchpath = BranchPath(ratepaths,ratetimes)
	for node in nodelist
		if !isroot(node)			
			parentnode = get(node.parent)

			paths = Array{Int,1}[]
			times = Array{Float64,1}[]
			for col=1:numcols
				parentstate = node_h[parentnode.nodeindex,col]
				nodestate = node_h[node.nodeindex,col]
				V = getQandPt(modelparams, 0, 0, node_aa[parentnode.nodeindex,col], node_rates[parentnode.nodeindex,col], 1.0)[1]
				path,time = modifiedrejectionsampling(rng, V*max(0.01,node.branchlength), parentstate, nodestate, nothing)
				push!(paths,path)
				push!(times,time)
			end

			node.data.branchpath = BranchPath(paths,times)

			aapaths = Array{Int,1}[]
			aatimes = Array{Float64,1}[]
			for col=1:numcols
				parentaastate = node_aa[parentnode.nodeindex,col]
				nodeaastate = node_aa[node.nodeindex,col]
				Q = getAAandPt(modelparams, 0, 0, node_h[parentnode.nodeindex,col], node_rates[parentnode.nodeindex,col], 1.0)[1]
				aapath,aatime = modifiedrejectionsampling(rng, Q*max(0.01,node.branchlength), parentaastate, nodeaastate, nothing)
				push!(aapaths,aapath)
				push!(aatimes,aatime)
			end
			node.data.aabranchpath = BranchPath(aapaths,aatimes)

			ratepaths = Array{Int,1}[]
			ratestimes = Array{Float64,1}[]
			for col=1:numcols
				#push!(ratepaths, Int[div(modelparams.numrates,2)])
				#push!(ratestimes, Float64[0.0])
				parentratesstate = node_rates[parentnode.nodeindex,col]
				noderatesstate = node_rates[node.nodeindex,col]
				Q = getRandPt(modelparams, 0, 0, node_h[parentnode.nodeindex,col], node_aa[parentnode.nodeindex,col], 1.0)[1]
				#println(Q)
				ratepath,ratetime = modifiedrejectionsampling(rng, Q*max(0.01,node.branchlength), parentratesstate, noderatesstate, nothing)
				push!(ratepaths,ratepath)
				push!(ratestimes,ratetime)
			end
			node.data.ratesbranchpath = BranchPath(ratepaths,ratestimes)

			#=
			jointpaths = Array{Int,1}[]
			jointtimes = Array{Float64,1}[]
			for col=1:numcols
				parentstate = (node_h[parentnode.nodeindex,col]-1)*modelparams.numhiddenstates + node_aa[parentnode.nodeindex,col]
				nodestate = (node_h[node.nodeindex,col]-1)*modelparams.numhiddenstates + node_aa[node.nodeindex,col]
				Q = getJointPt(modelparams, 0, 0, 1.0)[1]
				jointpath,jointtime = modifiedrejectionsampling(rng, Q*max(0.001,node.branchlength), parentstate, nodestate, nothing)
				push!(jointpaths,jointpath)
				push!(jointtimes,jointtime)
				node.data.aabranchpath.paths[col] = Int[(a-1)%modelparams.alphabet+1 for a in jointpath]
				node.data.aabranchpath.times[col] = copy(jointtime)
				node.data.branchpath.paths[col] = Int[div(a-1,modelparams.alphabet)+1 for a in jointpath]
				node.data.branchpath.times[col] = copy(jointtime)
				node.data.aabranchpath.paths[col], node.data.aabranchpath.times[col] = removevirtualjumps(node.data.aabranchpath.paths[col], node.data.aabranchpath.times[col])
				node.data.branchpath.paths[col], node.data.branchpath.times[col] = removevirtualjumps(node.data.branchpath.paths[col], node.data.branchpath.times[col])
			end
			node.data.jointbranchpath = BranchPath(jointpaths,jointtimes)
			=#
			
			#=

			col = rand(1:numcols)
			cols = [col]
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols), BranchPathIterator(node.data.ratesbranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]	
			ratesiter = multi_iter.branchpathiterators[3]
			println("START")
			println(col,"\t",node.data.branchpath.paths[col],"\t", node.data.branchpath.times[col])
			println(col,"\t",node.data.aabranchpath.paths[col],"\t", node.data.aabranchpath.times[col])
			println(col,"\t",node.data.ratesbranchpath.paths[col],"\t", node.data.ratesbranchpath.times[col])
			for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)
				curriter = multi_iter.branchpathiterators[multi_iter.branchpathindex]
				println(multi_iter.prevtime,"\t",multi_iter.currtime,"\t",curriter.prevstates,"\t",curriter.prevtime,"\t",curriter.currstates,"\t",curriter.currtime,"\t",curriter.mincol)
			end
			println("DONE")
			println("")=#
		end
	end
end

function initialise_tree(rng::AbstractRNG, modelparams::ModelParams, inputroot::TreeNode, name_protein_dict::Dict{String,Tuple{Int64,Protein}}, numcols::Int; midpointroot::Bool=true, scalebranchlengths::Bool=true)	
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

	if scalebranchlengths
		println("Branch length scaling: ", modelparams.branchscalingfactor)
		#=
		for node in nodelist
			node.branchlength *= modelparams.branchscalingfactor
		end=#
		modelparams.scalingfactor = modelparams.branchscalingfactor
	end

	initialise_tree_aa(rng, modelparams, nodelist, numcols) 
	#=
	V = getQandPt(modelparams, 0, 0, 1,2, -1.0)[1]
	paths = Array{Int,1}[]
	times = Array{Float64,1}[]
	for col=1:numcols
		state = rand(rng,1:modelparams.numhiddenstates)
		push!(paths,Int[state,state])
		push!(times,Float64[0.0, 1.0])
	end	
	aapaths = Array{Int,1}[]
	aatimes = Array{Float64,1}[]
	for col=1:numcols
		aastate = rand(rng,1:modelparams.alphabet)
		push!(aapaths,Int[aastate,aastate])
		push!(aatimes,Float64[0.0, 1.0])
	end	
	root.data.branchpath = BranchPath(paths,times)
	root.data.aabranchpath = BranchPath(aapaths,aatimes)
	for node in nodelist
		if !isroot(node)			
			parentnode = get(node.parent)

			paths = Array{Int,1}[]
			times = Array{Float64,1}[]
			for col=1:numcols
				parentstate = parentnode.data.branchpath.paths[col][end]
				nodestate = rand(rng,1:modelparams.numhiddenstates)
				path,time = modifiedrejectionsampling(rng, V*0.5, parentstate, nodestate, nothing)
				push!(paths,path)
				push!(times,time)
			end

			aapaths = Array{Int,1}[]
			aatimes = Array{Float64,1}[]
			for col=1:numcols
				parentaastate = parentnode.data.aabranchpath.paths[col][end]
				nodeaastate = rand(rng,1:modelparams.alphabet)
				if 1 <= col <= length(node.data.protein.sites) && node.data.protein.sites[col].aa > 0
					nodeaastate = node.data.protein.sites[col].aa
				end
				aapath,aatime = modifiedrejectionsampling(rng, modelparams.aminoacidQ, parentaastate, nodeaastate, nothing)
				push!(aapaths,aapath)
				push!(aatimes,aatime)
			end

			node.data.branchpath = BranchPath(paths,times)
			node.data.aabranchpath = BranchPath(aapaths,aatimes)
		end
	end=#

	return root,nodelist
end

function training_example_from_sequence_alignment(rng::AbstractRNG, modelparams::ModelParams, fastafile::String; newickfile=nothing, blindnodenames::Array{String,1}=String[], scalebranchlengths::Bool=true)
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
	    root,nodelist = initialise_tree(rng, modelparams, root, name_protein_dict, length(sequences[1]), scalebranchlengths=scalebranchlengths)
	else
		root = gettreefromnewick(readlines(open(newickfile,"r"))[1])
	    root,nodelist = initialise_tree(rng, modelparams, root, name_protein_dict, length(sequences[1]),midpointroot=false, scalebranchlengths=scalebranchlengths)
	end

	return ( proteins, nodelist, sequences)
end

function training_example_from_json_family(rng::AbstractRNG, modelparams::ModelParams, json_family; blindnodenames::Array{String,1}=String[], scalebranchlengths::Bool=true)
	root = gettreefromnewick(json_family["newick_tree"])
	binarize!(root)
	nodelist = getnodelist(root)
	for (index,node) in enumerate(nodelist)
		node.nodeindex = index
	end

	numcols = length(json_family["proteins"][1]["aligned_sequence"])
	
	#=
	V = constructHiddenMatrix(modelparams, ones(Float64,modelparams.numhiddenstates), ones(Float64,modelparams.numhiddenstates))
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
	end=#

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
		    	phi = phi_psi[1]
		    	psi = phi_psi[2]
		        push!(protein.sites, SiteObservation(0,indexof(string(aa), aminoacids),phi,omega,psi,bond_lengths[1],bond_lengths[2],bond_lengths[3],bond_angles[1],bond_angles[2],bond_angles[3]))
		    end
	    end
	    name_protein_dict[json_family["proteins"][p]["name"]] = (p, protein)
	    push!(proteins, protein)
	end

	root,nodelist = initialise_tree(rng, modelparams, root, name_protein_dict, numcols, midpointroot=false, scalebranchlengths=scalebranchlengths)
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

		multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col]), BranchPathIterator(node.data.ratesbranchpath,Int[col])])
		hiddeniter = multi_iter.branchpathiterators[1]
		aaiter = multi_iter.branchpathiterators[2]
		ratesiter = multi_iter.branchpathiterators[3]
		for it in multi_iter
			dt = (multi_iter.currtime-multi_iter.prevtime)
			
			changecol = selcol
			Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[1], aaiter.prevstates[1], ratesiter.prevstates[1], ratesiter.prevstates[1])
			exitrate +=  Qii*dt

			changecol = 0
			if multi_iter.branchpathindex == 1
				changecol = hiddeniter.mincol
				if changecol == selcol
					N += 1.0
				end
			elseif multi_iter.branchpathindex == 2
				changecol = aaiter.mincol
				N += 1.0
			else
				changecol = ratesiter.mincol
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

				multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col]), BranchPathIterator(node.data.ratesbranchpath,Int[col])])
				hiddeniter = multi_iter.branchpathiterators[1]
				aaiter = multi_iter.branchpathiterators[2]
				ratesiter = multi_iter.branchpathiterators[3]
				for it in multi_iter
					dt = (multi_iter.currtime-multi_iter.prevtime)
					
					changecol = selcol
					Qii = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[1], aaiter.prevstates[1], ratesiter.prevstates[1], ratesiter.prevstates[1])
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
					elseif multi_iter.branchpathindex == 3
						changecol = ratesiter.mincol
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

function samplepaths(rng::AbstractRNG, col::Int,proteins,nodelist::Array{TreeNode,1}, modelparams::ModelParams)
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
	accepted_rates = 0.0
	accepted_rates_total = 0.0
	temppath = Array{Int,1}[node.data.branchpath.paths[col] for node in nodelist]
	temptimes = Array{Float64,1}[node.data.branchpath.times[col] for node in nodelist]
	site1 = augmentedloglikelihood_site(nodelist, col, modelparams)	
	felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams)
	prop1 = hidden_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
	site2 = augmentedloglikelihood_site(nodelist, col, modelparams)
	prop2 = hidden_proposal_likelihood(nodelist, col, modelparams)

	if exp((site2-site1)+(prop1-prop2)) > rand(rng)
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
	site1 = augmentedloglikelihood_site(nodelist, col, modelparams)	
	felsensteinresample_aa(rng, proteins, nodelist, cols,col, modelparams)
	prop1 = aa_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
	site2 = augmentedloglikelihood_site(nodelist, col, modelparams)
	prop2 = aa_proposal_likelihood(nodelist, col, modelparams)
	#println(site1,"\t",site2,"\t",prop1,"\t",prop2,"\t",site2-site1,"\t",prop1-prop2,"\t",site2-site1+prop1-prop2)

	if exp((site2-site1)+(prop1-prop2)) > rand(rng)
		accepted_aa += 1.0
	else
		index = 1
		for node in nodelist
			#println("a", length(node.data.aabranchpath.paths))
			#println("b", length(temppath))
			#println(length(nodelist))
			node.data.aabranchpath.paths[col] = temppath[index]
			node.data.aabranchpath.times[col] = temptimes[index]
			index += 1
		end
	end
	accepted_aa_total += 1.0

	#=
	for node in nodelist
		if !isroot(node)
			println(modelparams.rates)
			println(rates_helper(node, col, modelparams))
		end
	end=#
	temppath = Array{Int,1}[node.data.ratesbranchpath.paths[col] for node in nodelist]
	temptimes = Array{Float64,1}[node.data.ratesbranchpath.times[col] for node in nodelist]
	site1 = augmentedloglikelihood_site(nodelist, col, modelparams)	
	felsensteinresample_rates(rng, proteins, nodelist, cols,col, modelparams)
	prop1 = rates_proposal_likelihood(nodelist, col, modelparams, temppath, temptimes)
	site2 = augmentedloglikelihood_site(nodelist, col, modelparams)
	prop2 = rates_proposal_likelihood(nodelist, col, modelparams)
	#println(site1,"\t",site2,"\t",prop1,"\t",prop2,"\t",site2-site1,"\t",prop1-prop2,"\t",site2-site1+prop1-prop2)

	if exp((site2-site1)+(prop1-prop2)) > rand(rng)
		accepted_rates += 1.0
	else
		index = 1
		for node in nodelist
			node.data.ratesbranchpath.paths[col] = temppath[index]
			node.data.ratesbranchpath.times[col] = temptimes[index]
			index += 1
		end
	end
	accepted_rates_total += 1.0
	

	return accepted_hidden,accepted_hidden_total,accepted_aa,accepted_aa_total,accepted_rates,accepted_rates_total
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

function train(numhiddenstates::Int=5)
	#rng = MersenneTwister(10498012421321)
	#Random.seed!(1234)
	rng = MersenneTwister(15149874541631)
	Random.seed!(rand(rng,UInt))


	learnrates = true
	samplebranchlengths = true
	family_dir = "../data/families/"
	family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))

	modelparams = ModelParams(LGmatrix,numhiddenstates,1.0)
	maxloglikelihood = -Inf

	outputmodelname = string("_h.",modelparams.numhiddenstates,"_learnrates.",learnrates)
	modelfile = string("model", outputmodelname, ".model")

	loadfromcache = false
	if loadfromcache
		modelparams = Serialization.deserialize(open(modelfile, "r"))
	end

	println("Data files: ", length(family_files))
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
				training_example = (training_example[1], TreeNode[root], training_example[3], training_example[4])
			end
			push!(trainingexamples, training_example)
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
				cols = Int[]
				if col > 1
					push!(cols,col-1)
					selcol = 2
				end
				push!(cols,col)
				if col < numcols
					push!(cols,col+1)
				end
				felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams)
				felsensteinresample_aa(rng, proteins, nodelist, cols,col, modelparams)
			end
		end
	end	

	logwriter = open(string("trace", outputmodelname,".log"), "w")
	println(logwriter, "iter\ttotalll\tpathll\tobservationll")
	for iter=1:1000
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

		ratetransitionrate_events = ones(Float64, modelparams.numrates, modelparams.numrates)*0.01
		ratetransitionrate_times = ones(Float64, modelparams.numrates, modelparams.numrates)*0.001

		aatransitionrate_counts = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		aatransitionrate_totals = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		for aa=1:modelparams.alphabet
			aatransitionrate_counts[aa,aa] = 0.0
		end		

		speed_transition_counts =  ones(Float64, modelparams.numrates, modelparams.numrates)*0.01
		rate_counts = ones(Float64,modelparams.numrates)*0.1

		totalbranchlength_output = 0.0
		accepted_hidden = 0.0
		accepted_hidden_total = 0.0
		accepted_aa = 0.0
		accepted_aa_total = 0.0
		accepted_rates = 0.0
		accepted_rates_total = 0.0
		for (proteins,nodelist,json_family,sequences) in trainingexamples
			numcols = length(proteins[1])
			
			for i=1:2
				randcols = shuffle(rng, Int[i for i=1:numcols])
				for col in randcols
					a1,a2,a3,a4,a5,a6 = samplepaths(rng,col,proteins,nodelist, modelparams)
					accepted_hidden += a1
					accepted_hidden_total += a2
					accepted_aa += a3
					accepted_aa_total += a4
					accepted_rates += a5
					accepted_rates_total += a6
					if rand(rng) < 0.0001
						println("START PATHS")
						for node in nodelist
							println(node.data.ratesbranchpath.paths[col],"\t",node.data.ratesbranchpath.times[col])
						end
					end
				end



				if samplebranchlengths
					for node in nodelist
						if !isroot(node)
							t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
							node.branchlength = t

							#=
							println(iter," branches ", node.branchlength,"\t",t)				
							oldll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
							node.branchlength = t
							newll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
							println(iter," likeliho ", newll,"\t",oldll,"\t",propratio,"\t",newll-oldll+propratio)
							=#
							events = mean([length(p)-1.0 for p in node.data.branchpath.paths])
							aa_events = mean([length(p)-1.0 for p in node.data.aabranchpath.paths])
							#println(node.nodeindex,"\t",node.data.inputbranchlength,"\t",node.branchlength,"\t",events,"\t",aa_events)
						end
					end
				end
			end

			for col=1:numcols-1
				for node in nodelist
					branchiterator = BranchPathIterator(node.data.branchpath, Int[col,col+1])
					for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
						modelparams.transitioncounts[prevstates[1],prevstates[2]] += 1.0
					end
				end
			end

			for col=1:numcols-1
				for node in nodelist
					branchiterator = BranchPathIterator(node.data.ratesbranchpath, Int[col,col+1])
					for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
						speed_transition_counts[prevstates[1],prevstates[2]] += 1.0
					end
				end
			end

			for col=1:numcols
				for node in nodelist
					branchiterator = BranchPathIterator(node.data.ratesbranchpath, Int[col])
					for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
						if changecol > 0							
							rate_counts[currstates[changecol]] += 1.0
						end
					end
				end
			end
			modelparams.rate_freqs = rate_counts./sum(rate_counts)

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
		println("Acceptance:\t",accepted_hidden/accepted_hidden_total,"\t",accepted_aa/accepted_aa_total,"\t",accepted_rates/accepted_rates_total)

		for training_example in trainingexamples
			for node in training_example[2]
				if !isroot(node)
					totalbranchlength_output += node.branchlength
				end
			end
		end
		modelparams.branchscalingfactor = totalbranchlength_output/totalbranchlength_input

		for r=1:modelparams.numrates
			speed_transition_counts[r,:] = speed_transition_counts[r,:] ./ sum(speed_transition_counts[r,:])	
		end

		for h=1:modelparams.numhiddenstates	
			estimate_hidden_transition_probs(modelparams)
			estimate_categorical(modelparams.hiddennodes[h].aa_node, 1.0)
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
		println("SPEED ", speed_transition_counts)		
		println(modelparams.rate_freqs)

		if learnrates
			rate_events = 0.0
			rate_times = 0.0
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
							
							multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col]), BranchPathIterator(node.data.ratesbranchpath,Int[col])])
							hiddeniter = multi_iter.branchpathiterators[1]
							aaiter = multi_iter.branchpathiterators[2]
							ratesiter = multi_iter.branchpathiterators[3]
							for it in multi_iter
								dt = (multi_iter.currtime-multi_iter.prevtime)
								
								changecol = selcol

								for r=1:modelparams.numrates
									thisr = ratesiter.prevstates[1]
									if thisr != r
										rateentry = getrateentry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], aaiter.prevstates[1], thisr, r)
										
										#rate_times += rateentry*dt*node.branchlength/modelparams.rate_mu
										ratetransitionrate_times[thisr,r] += rateentry*dt*node.branchlength/modelparams.rate_exchangeablities[thisr, r]
									end
								end

								rate_entry = entry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[1], aaiter.prevstates[1], ratesiter.prevstates[1], ratesiter.prevstates[1])
								rate_cat_times[ratesiter.prevstates[1]] += -rate_entry*dt*node.branchlength/modelparams.rates[ratesiter.prevstates[1]]

								if multi_iter.branchpathindex == 1 && hiddeniter.mincol == selcol
									changecol = hiddeniter.mincol
									rate_cat_events[ratesiter.prevstates[1]] += 1.0									
								elseif multi_iter.branchpathindex == 2
									rate_cat_events[ratesiter.prevstates[1]] += 1.0
								elseif multi_iter.branchpathindex == 3
									ratetransitionrate_events[ratesiter.prevstates[1], ratesiter.currstates[1]] += 1.0
									rate_cat_events[ratesiter.prevstates[1]] += 1.0
								end
							end
						end
					end
				end
			end

			ratetransitionrates = zeros(Float64, modelparams.numrates, modelparams.numrates)
			for r1=1:modelparams.numrates
				ratetransitionrates[r1,r1] = 0.0
				for r2=1:modelparams.numrates
					if r1 != r2
						ratetransitionrates[r1,r2] = (ratetransitionrate_events[r1,r2]+ratetransitionrate_events[r2,r1])/(ratetransitionrate_times[r1,r2]+ratetransitionrate_times[r2,r1])
						ratetransitionrates[r1,r1] -= ratetransitionrates[r1,r2]
					end
				end
			end	

			modelparams.rates = rate_cat_events./rate_cat_times		
			modelparams.rate_exchangeablities = ratetransitionrates
			println("Rates\t", ratetransitionrates)
			println(modelparams.rates)
		end

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
							
							multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[col]), BranchPathIterator(node.data.ratesbranchpath,Int[col])])
							hiddeniter = multi_iter.branchpathiterators[1]
							aaiter = multi_iter.branchpathiterators[2]
							ratesiter = multi_iter.branchpathiterators[3]
							for it in multi_iter
								dt = (multi_iter.currtime-multi_iter.prevtime)
								
								changecol = selcol
								for h=1:modelparams.numhiddenstates
									thish = hiddeniter.prevstates[changecol]
									if thish != h
										hiddenrateentry = gethiddenentry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), thish, h, aaiter.prevstates[1], ratesiter.prevstates[1])
										#transitionrate_times[thish] += hiddenrateentry*dt*node.branchlength/modelparams.transitionrates[thish,h]
										transitionrate_totals[thish, h] += hiddenrateentry*dt*node.branchlength/modelparams.transitionrates[thish,h]
									end
								end

								for aa=1:modelparams.alphabet
									thisaa = aaiter.prevstates[1]
									if thisaa != aa
										aaentry = getaaentry(modelparams, get(hiddeniter.prevstates, changecol-1, 0), get(hiddeniter.prevstates, changecol+1, 0), hiddeniter.prevstates[changecol], thisaa, aa, ratesiter.prevstates[1])
										aatransitionrate_times[thisaa, aa] += aaentry*dt*node.branchlength/modelparams.aa_exchangeablities[thisaa, aa]
									end
								end

								if multi_iter.branchpathindex == 1 && hiddeniter.mincol == selcol
									changecol = hiddeniter.mincol
									transitionrate_counts[hiddeniter.prevstates[changecol], hiddeniter.currstates[changecol]] += 1.0
									transitionrate_events[hiddeniter.prevstates[changecol]] += 1.0
								elseif multi_iter.branchpathindex == 2
									aatransitionrate_events[aaiter.prevstates[1], aaiter.currstates[1]] += 1.0
								elseif multi_iter.branchpathindex == 3
									#ratetransitionrate_events[ratesiter.prevstates[1], ratesiter.currstates[1]] += 1.0
									#rate_cat_events[ratesiter.prevstates[1]] += 1.0
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
						#modelparams.transitionrates[h,h2] = transitionrate_counts[h,h2]/transitionrate_times[h]
						modelparams.transitionrates[h,h2] = (transitionrate_counts[h,h2]+transitionrate_counts[h2,h])/(transitionrate_totals[h,h2]+transitionrate_totals[h2,h])
						modelparams.transitionrates[h,h] -= modelparams.transitionrates[h,h2]
					end
				end
			end

			#println(transitionrate_counts)
			#println(transitionrate_totals)
			#exit()
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
		end

		#modelparams.rate_freqs = rate_counts./sum(rate_counts)

		#println(rate_cat_events./rate_cat_times)
		#optimizerates2(trainingexamples, modelparams)
		#optimizerates3(trainingexamples, modelparams)

		
		#
		#=
		if (iter) % 5 == 0
			optimizerates(trainingexamples, modelparams)
		end
		=#

		augmentedll_end,observationll_end = calculateloglikelihood(modelparams, trainingexamples)
		totalll_end = augmentedll_end+observationll_end
		

		println(logwriter, iter-1,"\t",totalll_end,"\t",augmentedll_end,"\t",observationll_end)
		flush(logwriter)
		reset_matrix_cache(modelparams)
		if 1 == 1 || totalll_end > maxloglikelihood
			maxloglikelihood = totalll_end			
			fout = open(modelfile, "w")
			Serialization.serialize(fout, modelparams)
			close(fout)
		else 
			fin = open(modelfile, "r")	
			modelparams = Serialization.deserialize(fin)
			close(fin)
		end
	end
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

function rates_helper2(x::Array{Float64,1}, trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	modelparams.rate_mu = x[1]
	modelparams.rate_alpha = x[2]
	modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)
	augmentedll = 0.0
	for (proteins,nodelist) in trainingexamples	
		numcols = length(proteins[1])
		augmentedll += augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
	end
	println(augmentedll,"\t",x,"\t",modelparams.rates)
	return augmentedll
end

function optimizerates2(trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	opt = Opt(:LN_COBYLA,2)
    localObjectiveFunction = ((param, grad) -> rates_helper2(param, trainingexamples, modelparams))
    lower = ones(Float64, 2)*1e-2
    upper = ones(Float64, 2)*1e6
    lower[2] = 1.0
    upper[2] = 1.0
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, 25)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, Float64[modelparams.rate_mu, modelparams.rate_alpha])
	modelparams.rate_mu = minx[1]
	modelparams.rate_alpha = minx[2]
end

function rates_helper3(x::Array{Float64,1}, trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	modelparams.rate_alpha = x[1]
	modelparams.rates = discretizegamma(modelparams.rate_alpha, 1.0/modelparams.rate_alpha, modelparams.numrates)
	augmentedll = 0.0
	for (proteins,nodelist) in trainingexamples	
		numcols = length(proteins[1])
		augmentedll += augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
	end
	augmentedll = augmentedll - 2.0*x[1]
	println(augmentedll,"\t",x,"\t",modelparams.rates)
	return augmentedll
end

function optimizerates3(trainingexamples::Array{Tuple,1}, modelparams::ModelParams)
	opt = Opt(:LN_COBYLA,1)
    localObjectiveFunction = ((param, grad) -> rates_helper3(param, trainingexamples, modelparams))
    lower = ones(Float64, 1)*1e-2
    upper = ones(Float64, 1)*1e6
    lower_bounds!(opt, lower)
    upper_bounds!(opt, upper)
    xtol_rel!(opt,1e-5)
    maxeval!(opt, 15)
    max_objective!(opt, localObjectiveFunction)
    (minf,minx,ret) = optimize(opt, Float64[modelparams.rate_alpha])
	modelparams.rate_alpha = minx[1]
end

using TraitAssociation

function infer(;modelfile="model_h.10_learnrates.true.model", samplebranchlengths::Bool=false)
	rng = MersenneTwister(10498012421321)
	Random.seed!(1234)

	#fin = open("model_h_25.model", "r")
	blindscores = Dict{String,Array{Float64,1}}()
	blindnodenames = String[]

	fin = open(modelfile, "r")	
	modelparams = Serialization.deserialize(fin)
	close(fin)

	println("Initialising...")
	
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
	newickfile= abspath("tree.mean.branch.consensus.nwk")		
	blindnodenames = String["6n41.pdb_1951"]
	#blindnodenames = String[]

	
	#fastafile = abspath("../data/hiv/curated6.fasta")
	#newickfile=abspath("../data/hiv/curated6.nwk")	
	#newickfile=abspath("tree.mean.branch.consensus.nwk")	
	#blindnodenames = String["B.US.1978.SF4.KJ704795"]
	#blindnodenames = String[]

	scalebranchlengths=true
	proteins,nodelist,sequences = training_example_from_sequence_alignment(rng, modelparams, fastafile, newickfile=newickfile, blindnodenames=blindnodenames, scalebranchlengths=scalebranchlengths)
	
	if length(blindnodenames) > 0
		LGreconstruction_score = TraitAssociation.pathreconstruction(fastafile,newickfile,blindnodenames,1)
		println("LGreconstruction_score ",LGreconstruction_score)
	end


	outputprefix = "tree"

	println("Initialisation finished.")

	inputnodelist = deepcopy(nodelist)
	inputtreewriter = open("$(outputprefix).input.nwk", "w")
	println(inputtreewriter, getnewick(nodelist[1]))
	close(inputtreewriter)

	

	numcols = length(proteins[1])
	mcmcwriter = open("$(outputprefix).log", "w")
	treewriter = open("$(outputprefix).mcmc.log", "w")
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
	subspersite_cache = Dict{Int,Array{Float64,1}}()
	branchlength_cache = Dict{Int,Array{Float64,1}}()
	for iter=1:10000

		if samplebranchlengths
			println("$(iter).1 Sampling branch lengths START")
			randindex = rand(2:length(nodelist))
			for node in nodelist
				if !isroot(node)	
					t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
					#println(iter," branches ", node.branchlength,"\t",t)				
					#oldll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
					node.branchlength = t
					#newll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
					#println(iter," likeliho ", newll,"\t",oldll,"\t",propratio,"\t",newll-oldll+propratio)
				end
			end
			println("$(iter).1 Sampling branch lengths DONE")
		elseif iter > 30
			#pre = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
			#propratio = proposescalingfactor(rng, nodelist, Int[col for col=1:numcols], modelparams)
			##post = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
			#println("LIKS ",(post-pre),"\t",propratio)
		end
		modelparams.scalingfactor = 1.0
		#println(totalexitrate,"\t",totalN)
		
		println("$(iter).2 Sampling sites START")
		randcols = shuffle(rng, Int[i for i=1:numcols])
		for col in randcols
			

			samplepaths(rng, col, proteins,nodelist, modelparams)
			#felsensteinresample(rng, proteins, nodelist, col, cols, col, modelparams)
			#felsensteinresample_aa(rng, proteins, nodelist, col, cols,col, modelparams)
		end
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
		observationll = observationloglikelihood(proteins, nodelist, modelparams)
		print(mcmcwriter, iter-1,"\t",augmentedll+observationll,"\t",augmentedll,"\t",observationll)
		for node in nodelist
			if !isroot(node)				
				cached_branchlengths = get(branchlength_cache, node.nodeindex, Float64[])
				push!(cached_branchlengths, node.branchlength)
				branchlength_cache[node.nodeindex] = cached_branchlengths

				print(mcmcwriter,"\t$(node.branchlength)")
			end
		end
		println("$(iter).4 Calculating likelihoods DONE")	

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


using PyPlot
function plot_nodes(modelfile)
	fin = open(modelfile, "r")
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
#plot_nodes("model_h.15_learnrates.true.model")
#infer()
