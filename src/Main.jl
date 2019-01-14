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
using Formatting

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



mutable struct ModelParams
    alphabet::Int
    aminoacidQ::Array{Float64,2}
    aa_exchangeablities::Array{Float64,2}
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
        new(20,aminoacidQ,LGexchangeability,numhiddenstates,initialprobs,transitionprobs,transitionrates,hiddennodes,mu,matrixcache)
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


function getAAandPt(modelparams::ModelParams, prevh::Int, nexth::Int, prevstateh::Int, currstateh::Int, t::Float64)
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
		Q = constructAAMatrix(modelparams, prevprobs, nextprobs, prevstateh, currstateh)
		
		#decomposition = eigen(Q)
		#D, V = decomposition.values, decomposition.vectors
		#Vi = inv(V)
		#modelparams.matrixcache[key] = (Q,V,D,Vi)
		return Q,exp(Q*t)

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

function getQandPt(modelparams::ModelParams, prevh::Int, nexth::Int, prevaa::Int, curraa::Int, t::Float64)
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
		Q = constructHiddenMatrix(modelparams, prevprobs, nextprobs, prevaa, curraa)
		
		#decomposition = eigen(Q)
		#D, V = decomposition.values, decomposition.vectors
		#Vi = inv(V)
		#modelparams.matrixcache[key] = (Q,V,D,Vi)
		return Q,exp(Q*t)

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
	dummy::Int
	protein::Protein

	AugmentedNodeData(col::Int) = new(BranchPath(col),BranchPath(col),1, Protein()) 
	AugmentedNodeData(branchpath::BranchPath, aabranchpath::BranchPath, dummy::Int) = new(branchpath, aabranchpath, dummy, Protein())
end

function gettransprobs_aa(node::TreeNode, selcolin::Int, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol])])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	P = Matrix{Float64}(I, modelparams.alphabet, modelparams.alphabet)
	Pmatrices = Array{Float64,2}[]
	dummypath = Int[]
	dummytime = Float64[]

	for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		prevh = 0
		succh = 0
		changecol = hiddeniter.mincol
		if changecol == 2
			prevh = hiddeniter.currstates[1]
		end
		if length(hiddeniter.currstates) == 3
			succh = hiddeniter.currstates[3]
		end	
		R,Pi = getAAandPt(modelparams, prevh, succh, hiddeniter.prevstates[selcol], hiddeniter.currstates[selcol], dt)

		P *= Pi
		push!(dummypath,0)
		push!(dummytime,multi_iter.prevtime)
	end
	return P
end

function gettransprobs(node::TreeNode, selcolin::Int, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol])])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	P = Matrix{Float64}(I, modelparams.numhiddenstates, modelparams.numhiddenstates)
	Pmatrices = Array{Float64,2}[]
	dummypath = Int[]
	dummytime = Float64[]

	for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		prevh = 0
		succh = 0
		changecol = hiddeniter.mincol
		if changecol == 2
			prevh = hiddeniter.currstates[1]
		end
		if length(hiddeniter.currstates) == 3
			succh = hiddeniter.currstates[3]
		end		
		R,Pi = getQandPt(modelparams, prevh, succh, aaiter.prevstates[1], aaiter.currstates[1], dt)

		P *= Pi
		push!(dummypath,0)
		push!(dummytime,multi_iter.prevtime)
	end
	return P
end

function felsensteinhelper_aa(node::TreeNode, selcolin::Int, cols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol])])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	P = Matrix{Float64}(I, modelparams.numhiddenstates, modelparams.numhiddenstates)
	Rmatrices = Array{Float64,2}[]
	Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]
	
    for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		prevh = 0
		succh = 0
		changecol = hiddeniter.mincol
		if changecol == 2
			prevh = hiddeniter.currstates[1]
		end
		if length(hiddeniter.currstates) == 3
			succh = hiddeniter.currstates[3]
		end		
		R,Pi = getAAandPt(modelparams, prevh, succh, hiddeniter.prevstates[selcol], hiddeniter.currstates[selcol], dt)

		push!(Rmatrices, R*dt)
    	push!(Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)
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
    return Pmatrices,vs
end

function felsensteinhelper(node::TreeNode, selcolin::Int, cols::Array{Int,1}, aacol::Int, v::Array{Float64,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,Int[aacol])])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]
	P = Matrix{Float64}(I, modelparams.numhiddenstates, modelparams.numhiddenstates)
	Rmatrices = Array{Float64,2}[]
	Pmatrices = Array{Float64,2}[]
	vs = Array{Float64,1}[]
	dummytime = Float64[]

    for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
		prevh = 0
		succh = 0
		changecol = hiddeniter.mincol
		if changecol == 2
			prevh = hiddeniter.currstates[1]
		end
		if length(hiddeniter.currstates) == 3
			succh = hiddeniter.currstates[3]
		end		
		R,Pi = getQandPt(modelparams, prevh, succh, aaiter.prevstates[1], aaiter.currstates[1], dt)

		push!(Rmatrices, R*dt)
    	push!(Pmatrices,Pi)
    	push!(dummytime,multi_iter.prevtime)
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
    return Pmatrices,vs
end

function felsensteinresample_aa(rng::AbstractRNG, proteins::Array{Protein,1}, nodelist::Array{TreeNode,1}, selcolin::Int, cols::Array{Int,1}, aacol::Int, modelparams::ModelParams)
	#selcol = findfirst(x -> x == selcolin, cols)
	selcol = selcolin
	likelihoods = ones(Float64, length(nodelist), modelparams.alphabet)*-Inf
	logm = zeros(Float64,length(nodelist))

	stack = Int[1]
	while length(stack) > 0
		nodeindex = stack[end]
		node = nodelist[nodeindex]
		if isleafnode(node)
			if 1 <= selcolin <= length(node.data.protein.sites) && node.data.protein.sites[selcolin].aa > 0
				for a=1:modelparams.alphabet
					likelihoods[nodeindex,a] = 0.0
				end
				likelihoods[nodeindex,node.data.protein.sites[selcolin].aa] = 1.0
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
                lefttransprobs = gettransprobs_aa(nodelist[leftchildindex], selcol, cols, aacol, modelparams)
				righttransprobs = gettransprobs_aa(nodelist[rightchildindex], selcol, cols, aacol, modelparams)

        		#likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])
        		likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])

				Pmatrices_left,vs_left = felsensteinhelper_aa(nodelist[leftchildindex], selcol, cols, aacol, likelihoods[leftchildindex,:], modelparams)
				Pmatrices_right,vs_right = felsensteinhelper_aa(nodelist[rightchildindex], selcol, cols, aacol, likelihoods[rightchildindex,:], modelparams)

				m = maximum(likelihoods[nodeindex,:])
				likelihoods[nodeindex,:] = likelihoods[nodeindex,:] ./ m
				logm[nodeindex] = log(m) + logm[leftchildindex] + logm[rightchildindex]
        		pop!(stack)
        	end
        end
    end

	rootnode = nodelist[1]
	len = length(rootnode.data.aabranchpath.paths)
	h = rootnode.data.branchpath.paths[selcol][1]
	rootliks = modelparams.hiddennodes[h].aa_node.probs.*likelihoods[1,:]
	rootstate = CommonUtils.sample(rng,rootliks)
	rootnode.data.aabranchpath.paths[selcol] = Int[rootstate]
	rootnode.data.aabranchpath.times[selcol] = Float64[0.0]
	print = false
	backwardsampling_aa(rng,nodelist[1], rootstate, selcol,likelihoods,print,modelparams)
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
                lefttransprobs = gettransprobs(nodelist[leftchildindex], selcol, cols, aacol, modelparams)
				righttransprobs = gettransprobs(nodelist[rightchildindex], selcol, cols, aacol, modelparams)

        		#likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])
        		likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:]).*observationlikelihood(node.data.protein, selcolin, modelparams)

				Pmatrices_left,vs_left = felsensteinhelper(nodelist[leftchildindex], selcol, cols, aacol, likelihoods[leftchildindex,:], modelparams)
				Pmatrices_right,vs_right = felsensteinhelper(nodelist[rightchildindex], selcol, cols, aacol, likelihoods[rightchildindex,:], modelparams)

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
		prevh = rootnode.data.branchpath.paths[selcol-1][1]
	end
	nexth = 0
	if selcol < len
		nexth = rootnode.data.branchpath.paths[selcol+1][1]
	end
	freqs = getinitialsiteprobs(modelparams, prevh, nexth)
	rootliks = freqs.*likelihoods[1,:]
	rootstate = CommonUtils.sample(rng,rootliks)
	rootnode.data.branchpath.paths[selcol] = Int[rootstate]
	rootnode.data.branchpath.times[selcol] = Float64[0.0]
	print = false
	backwardsampling(rng,nodelist[1], rootstate, selcol,likelihoods,print,modelparams)
end

function backwardsampling_aa(rng::AbstractRNG,node::TreeNode, state::Int, selcol::Int,likelihoods,print::Bool,modelparams::ModelParams)
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
			samplepath, sampletimes = modifiedrejectionsampling(rng, child.data.aabranchpath.Rmatrices[z], path[z], path[z+1],(modelparams))
			append!(newpath,samplepath)
			append!(newtime,(sampletimes*dt) .+ child.data.aabranchpath.time[z])
		end

		newpath, newtime = removevirtualjumps(newpath, newtime)

		child.data.aabranchpath.paths[selcol] = newpath
		child.data.aabranchpath.times[selcol] = newtime
		backwardsampling_aa(rng,child, path[end],selcol, likelihoods,print,modelparams)
	end
end

function backwardsampling(rng::AbstractRNG,node::TreeNode, state::Int, selcol::Int,likelihoods,print::Bool,modelparams::ModelParams)
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
			samplepath, sampletimes = modifiedrejectionsampling(rng, child.data.branchpath.Rmatrices[z], path[z], path[z+1],(modelparams))
			append!(newpath,samplepath)
			append!(newtime,(sampletimes*dt) .+ child.data.branchpath.time[z])
		end

		newpath, newtime = removevirtualjumps(newpath, newtime)

		if print && isleafnode(child) && child.data.branchpath.paths[selcol][end] != newpath[end]
			println(isleafnode(child),"\t",selcol)
			println(child.data.branchpath.Pmatrices)
			println(child.data.branchpath.vs)
			println("J", newpath,newtime)
			println("K", child.data.branchpath.paths[selcol], child.data.branchpath.times[selcol])
		end
		child.data.branchpath.paths[selcol] = newpath
		child.data.branchpath.times[selcol] = newtime
		backwardsampling(rng,child, path[end],selcol, likelihoods,print,modelparams)
	end
end


function augmentedloglikelihood(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, modelparams::ModelParams)
	numcols = length(cols)
	loglikelihood = 0.0

	for node in nodelist
		if !isroot(node)
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]	

			for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
				
				for col=1:numcols
					prevprobs = ones(Float64, modelparams.numhiddenstates)
					succprobs = ones(Float64, modelparams.numhiddenstates)
					changecol = col
					if changecol > 1
						prevh = hiddeniter.currstates[changecol-1]
						prevprobs = modelparams.transitionprobs[prevh,:]
					end
					if changecol < numcols
						succh = hiddeniter.currstates[changecol+1]
						succprobs = modelparams.transitionprobs[:,succh]
					end
					Qii = entry(modelparams, prevprobs, succprobs, hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol])
					loglikelihood += Qii*dt
				end

				prevprobs = ones(Float64, modelparams.numhiddenstates)
				succprobs = ones(Float64, modelparams.numhiddenstates)
				changecol = 0
				if multi_iter.branchpathindex == 1
					changecol = hiddeniter.mincol
				else
					changecol = aaiter.mincol
				end

				if changecol > 0
					if changecol > 1
						prevh = hiddeniter.currstates[changecol-1]
						prevprobs = modelparams.transitionprobs[prevh,:]
					end
					if changecol < numcols
						succh = hiddeniter.currstates[changecol+1]
						succprobs = modelparams.transitionprobs[:,succh]
					end
					Qhi = entry(modelparams, prevprobs, succprobs, hiddeniter.prevstates[changecol], hiddeniter.currstates[changecol], aaiter.prevstates[changecol], aaiter.currstates[changecol])
					loglikelihood += log(Qhi) 
				end
			end
		end
	end

	#=
	for node in nodelist
		if !isroot(node)
			branchiterator = BranchPathIterator(node.data.branchpath,cols)

			for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
				Qii = 0.0
				len = length(cols)
				for selcol=1:len
					prevstate = prevstates[selcol]
					prevh = 0
					if selcol > 1
						prevh = prevstates[selcol-1]
					end
					nexth = 0
					if selcol < len
						nexth = prevstates[selcol+1]
					end
					#Qii += getratematrixrow(prevstates, selcol, params, modelspecification,componentindex)[prevstate]					
					Qii += getQandPt(modelparams, prevh, nexth,0,0, -1.0)[1][prevstate,prevstate]
					#entry(modelparams, modelparams.transitionprobs[prevh,:], modelparams.transitionprobs[:,succh], prevh::Int, currh::Int, prevaa::Int, curraa::Int)
				end
				dt = (currtime-prevtime)*node.branchlength
				if changecol > 0
					prevstate = prevstates[changecol]
					currstate = currstates[changecol]

					prevh = 0
					if changecol > 1
						prevh = currstates[changecol-1]
					end
					nexth = 0
					if changecol < len
						nexth = currstates[changecol+1]
					end
					#Qhi = getratematrixrow(prevstates, changecol, params, modelspecification,componentindex)[currstate]
					Qhi = getQandPt(modelparams, prevh, nexth,0,0, -1.0)[1][prevstate,currstate]
					loglikelihood += log(Qhi)
				end
				loglikelihood +=  Qii*dt
			end
		end
	end=#
	return loglikelihood
end

function getinitialsiteprobs(modelparams::ModelParams, prevh::Int, nexth::Int)
	prevprobs = ones(Float64, modelparams.numhiddenstates)
	if prevh > 0
		prevprobs = modelparams.transitionprobs[prevh,:]
	end	
	nextprobs = ones(Float64, modelparams.numhiddenstates)
	if nexth > 0
		nextprobs = modelparams.transitionprobs[:,nexth]
	end
	logfreqs = zeros(Float64,modelparams.numhiddenstates)
	for a=1:modelparams.numhiddenstates
		logfreqs[a] = log(prevprobs[a]*nextprobs[a]) # TODO: include initial probs
	end
	freqs = exp.(logfreqs.-maximum(logfreqs))
	freqs /= sum(freqs)
end

function getaaentry(modelparams::ModelParams, h::Int, prevaa::Int, curraa::Int)
	if prevaa != curraa
		return modelparams.aa_exchangeablities[prevaa,curraa]*modelparams.hiddennodes[h].aa_node.probs[curraa]
	else
		q = 0.0
		for aa=1:modelparams.alphabet
			if aa != prevaa
				q -= getaaentry(modelparamss, h, prevaa, aa)
			end
		end
		return q
	end
end

function entry(modelparams::ModelParams, prevprobs::Array{Float64,1}, succprobs::Array{Float64,1}, prevh::Int, currh::Int, prevaa::Int, curraa::Int)
	if prevh == currh && prevaa == curraa
		q = 0.0
		for h=1:modelparams.numhiddenstates
			for aa=1:modelparams.alphabet
				if h != prevh || aa != prevaa
					q -= entry(modelparams, prevprobs, succprobs, prevh, h, prevaa, aa)
				end
			end
		end
		return q
	elseif prevh == currh
		return getaaentry(modelparams, prevh, prevaa, curraa)
	elseif prevaa == curraa
		return prevprobs[currh]*succprobs[currh]*modelparams.mu*modelparams.transitionrates[prevh,currh]
	else
		return prevprobs[currh]*succprobs[currh]*modelparams.mu*modelparams.transitionrates[prevh,currh]*getaaentry(modelparams, currh, prevaa,curraa)
	end
end

function constructAAMatrix(modelparams::ModelParams, prevprobs::Array{Float64,1}, succprobs::Array{Float64,1}, prevh::Int, currh::Int)
	Q = zeros(Float64, modelparams.alphabet, modelparams.alphabet)
	for aa1=1:modelparams.alphabet
        for aa2=1:modelparams.alphabet
			if aa1 != aa2
				Q[aa1,aa2] = entry(modelparams, prevprobs, succprobs, prevh, currh, aa1, aa2)
				Q[aa1,aa1] -= Q[aa1,aa2]
			end
		end
	end

    return Q
end

function constructHiddenMatrix(modelparams::ModelParams, prevprobs::Array{Float64,1}, succprobs::Array{Float64,1}, prevaa::Int, curraa::Int)
	#=
	for prevh=1:modelparams.numhiddenstates
		for currh=1:modelparams.numhiddenstates
			for prevaa=1:modelparams.alphabet
				for curraa=1:modelparams.alphabet
					println("$prevh $currh $prevaa $curraa ", entry(modelparams, prevprobs, succprobs, prevh, currh, prevaa, curraa))
				end
			end
		end
	end=#
	

	Q = zeros(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)
	for h1=1:modelparams.numhiddenstates
        for h2=1:modelparams.numhiddenstates
			if h1 != h2
				Q[h1,h2] = entry(modelparams, prevprobs, succprobs, h1, h2, prevaa, curraa)
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
	end

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

			#=
			node.data.branchpath.paths[1] = Int[1] 
			node.data.branchpath.paths[2] = Int[2]
			node.data.branchpath.paths[3] = Int[3]
			node.data.branchpath.times[1] = Float64[0.0] 
			node.data.branchpath.times[2] = Float64[0.0] 
			node.data.branchpath.times[3] = Float64[0.0] 
			node.data.aabranchpath.paths[2] = Int[20]
			node.data.aabranchpath.times[2] = Float64[0.0]

			println("M",length(paths))
			println("N",length(aapaths))

			it = BranchPathIterator(node.data.branchpath,Int[1,2,3])
			aait = BranchPathIterator(node.data.aabranchpath,Int[2])
			multit = MultiBranchPathIterator(BranchPathIterator[it,aait])

			for index in it.indices
				println(it.branch.paths[index],"\t", it.branch.times[index])
			end
			println(aait.branch.paths[2],"\t", aait.branch.times[2])
			for el in multit
				println("Z")
			end		
			println("HERE")=#
		end
	end

	return root,nodelist
end

function training_example_from_sequence_alignment(rng::AbstractRNG, modelparams::ModelParams, fastafile::String,newickfile::String=nothing)
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
		    if name == "6n41.pdb_1951"
		    	site.aa = 0
		    	println("here")
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

	return ( proteins, nodelist)
end

function training_example_from_json_family(rng::AbstractRNG, modelparams::ModelParams, json_family)
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
	for p=1:length(json_family["proteins"])
	    protein = Protein()
	    name = json_family["proteins"][p]["name"]
	    json_protein = json_family["proteins"][p]
	    for (aa,phi_psi,omega,bond_lengths,bond_angles) in zip(json_protein["aligned_sequence"], json_protein["aligned_phi_psi"], json_protein["aligned_omega"], json_protein["aligned_bond_lengths"], json_protein["aligned_bond_angles"])
	    	aaindex  = 0
	    	if name != "1yaaa"
	    		aaindex = indexof(string(aa), aminoacids)
		    end
	    	phi = phi_psi[1]
	    	psi = phi_psi[2]
	        push!(protein.sites, SiteObservation(0,aaindex,phi,omega,psi,bond_lengths[1],bond_lengths[2],bond_lengths[3],bond_angles[1],bond_angles[2],bond_angles[3]))
	    end
	    name_protein_dict[json_family["proteins"][p]["name"]] = (p, protein)
	    push!(proteins, protein)
	end

	root,nodelist = initialise_tree(rng, modelparams, root, name_protein_dict, numcols, midpointroot=false)
	return (proteins, nodelist,json_family)
end

#=
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
			Qii += constructHiddenMatrix(modelparams, prevprobs, nextprobs,1,2)[prevstate,prevstate]
		end
		dt = (currtime-prevtime)
		exitrate +=  Qii*dt
	end
	return exitrate
end
=#

function getexitrate(node::TreeNode, cols::Array{Int,1}, modelparams::ModelParams)
	numcols = length(cols)
	multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
	hiddeniter = multi_iter.branchpathiterators[1]
	aaiter = multi_iter.branchpathiterators[2]	
	exitrate = 0.0
	for it in multi_iter
		dt = (multi_iter.currtime-multi_iter.prevtime)
		
		for col=1:numcols
			prevprobs = ones(Float64, modelparams.numhiddenstates)
			succprobs = ones(Float64, modelparams.numhiddenstates)
			changecol = col
			if changecol > 1
				prevh = hiddeniter.currstates[changecol-1]
				prevprobs = modelparams.transitionprobs[prevh,:]
			end
			if changecol < numcols
				succh = hiddeniter.currstates[changecol+1]
				succprobs = modelparams.transitionprobs[:,succh]
			end
			Qii = entry(modelparams, prevprobs, succprobs, hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol])
			exitrate +=  Qii*dt
		end
	end
	#=
	println("EXIT", exitrate)
	if exitrate == 0.0		
		exit()
	end=#
	return exitrate
end

function getexitrates(node::TreeNode, cols::Array{Int,1}, modelparams::ModelParams)
	len = length(cols)
	exitrate = 0.0
	numcols = length(cols)
	loglikelihood = 0.0

	for node in nodelist
		if !isroot(node)
			multi_iter = MultiBranchPathIterator(BranchPathIterator[BranchPathIterator(node.data.branchpath,cols), BranchPathIterator(node.data.aabranchpath,cols)])
			hiddeniter = multi_iter.branchpathiterators[1]
			aaiter = multi_iter.branchpathiterators[2]	

			for it in multi_iter
				dt = (multi_iter.currtime-multi_iter.prevtime)*node.branchlength
				
				for col=1:numcols
					prevprobs = ones(Float64, modelparams.numhiddenstates)
					succprobs = ones(Float64, modelparams.numhiddenstates)
					changecol = col
					if changecol > 1
						prevh = hiddeniter.currstates[changecol-1]
						prevprobs = modelparams.transitionprobs[prevh,:]
					end
					if changecol < numcols
						succh = hiddeniter.currstates[changecol+1]
						succprobs = modelparams.transitionprobs[:,succh]
					end
					Qii = entry(modelparams, prevprobs, succprobs, hiddeniter.prevstates[changecol], hiddeniter.prevstates[changecol], aaiter.prevstates[changecol], aaiter.prevstates[changecol])
					loglikelihood += Qii*dt
				end

				prevprobs = ones(Float64, modelparams.numhiddenstates)
				succprobs = ones(Float64, modelparams.numhiddenstates)
				changecol = 0
				if multi_iter.branchpathindex == 1
					changecol = hiddeniter.mincol
				else
					changecol = aaiter.mincol
				end

				if changecol > 0
					if changecol > 1
						prevh = hiddeniter.currstates[changecol-1]
						prevprobs = modelparams.transitionprobs[prevh,:]
					end
					if changecol < numcols
						succh = hiddeniter.currstates[changecol+1]
						succprobs = modelparams.transitionprobs[:,succh]
					end
					Qhi = entry(modelparams, prevprobs, succprobs, hiddeniter.prevstates[changecol], hiddeniter.currstates[changecol], aaiter.prevstates[changecol], aaiter.currstates[changecol])
					loglikelihood += log(Qhi) 
				end
			end
		end
	end
end

function proposebranchlength(rng::AbstractRNG, node::TreeNode, cols::Array{Int,1}, modelparams::ModelParams)
	totalexitrate = getexitrate(node, cols, modelparams)
	t =  log(1.0 - rand(rng))/totalexitrate
	#println(totalexitrate,"\tT\t",t)
	newll = log(-totalexitrate) + totalexitrate*t
	oldll = log(-totalexitrate) + totalexitrate*node.branchlength
	propratio = oldll-newll

	#=
	priorll = -5.0*t + 5.0*node.branchlength
	if exp(priorll) > rand(rng)

	else 
		t = node.branchlength
	end=#

	return t,propratio
end

function cosine_similarity(v1::Array{Float64,1}, v2::Array{Float64,1})
	return dot(v1,v2)/(norm(v1)*norm(v2))
end

function train(numhiddenstates::Int=5)
	rng = MersenneTwister(10498012421321)
	Random.seed!(1234)

	learnrates = true
	samplebranchlengths = true

	family_dir = "../data/families/"
	family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))

	modelparams = ModelParams(LGmatrix,numhiddenstates,0.2)
	#modelparams = Serialization.deserialize(open("model_h.25_learnrates.true.model", "r"))

	outputmodelname = string("_h.",modelparams.numhiddenstates,"_learnrates.",learnrates)

	trainingexamples = Tuple[]
	for family_file in family_files[1:200]
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

	logwriter = open(string("trace", outputmodelname,".log"), "w")
	println(logwriter, "iter\ttotalll\tpathll\tobservationll")
	for iter=1:1000
		reset_matrix_cache(modelparams)
		aacounts = ones(Float64, modelparams.numhiddenstates, 20)*0.1
		transitioncounts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
		
		transitionrate_counts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
		transitionrate_totals = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.01
		for h=1:modelparams.numhiddenstates
			transitionrate_counts[h,h] = 0.0
		end

		aatransitionrate_counts = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.1
		aatransitionrate_totals = ones(Float64, modelparams.alphabet, modelparams.alphabet)*0.01
		for aa=1:modelparams.alphabet
			aatransitionrate_counts[aa,aa] = 0.0
		end

		for (proteins,nodelist,json_family) in trainingexamples
			numcols = length(proteins[1])

			if samplebranchlengths
				for node in nodelist
					if !isroot(node)
						t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
						node.branchlength = t
					end
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
				felsensteinresample(rng, proteins, nodelist, col, cols,col, modelparams)
				felsensteinresample_aa(rng, proteins, nodelist, col, cols,col, modelparams)
			end

			for col=1:numcols-1
				for node in nodelist
					branchiterator = BranchPathIterator(node.data.branchpath, Int[col,col+1])
					for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
						transitioncounts[prevstates[1],prevstates[2]] += 1.0
					end
				end
			end

			# count hidden transitions
			for col=1:numcols
				for node in nodelist
					if !isroot(node)
						branchiterator = BranchPathIterator(node.data.branchpath, Int[col])
						for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
							if changecol > 0
								transitionrate_counts[prevstates[changecol],currstates[changecol]] += 1.0/numcols
								dt = currtime-prevtime
								transitionrate_totals[prevstates[changecol],currstates[changecol]] += dt*node.branchlength/numcols
							end
						end
					end
				end
			end

			# count aa transitions
			for col=1:numcols
				for node in nodelist
					if !isroot(node)
						branchiterator = BranchPathIterator(node.data.aabranchpath, Int[col])
						for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
							if changecol > 0
								aatransitionrate_counts[prevstates[changecol],currstates[changecol]] += 1.0/numcols
								dt = currtime-prevtime
								aatransitionrate_totals[prevstates[changecol],currstates[changecol]] += dt*node.branchlength/numcols
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
			rowtotal = 0.0
			for h2=1:modelparams.numhiddenstates
				if h != h2
					rowtotal += transitionrate_counts[h,h2]
				else 
					transitionrate_counts[h,h2] = 0.0
				end
			end
			transitionrate_counts[h,:] /= rowtotal
		end

		#=
		transitionweights = ones(Float64,modelparams.numhiddenstates,modelparams.numhiddenstates)*-Inf		
		for h=1:modelparams.numhiddenstates
			maxweight = -Inf
			for h2=1:modelparams.numhiddenstates
				#beta = 10.0*modelparams.numhiddenstates
				beta = 5.0*iter
				aa_profile_similarity = cosine_similarity(modelparams.hiddennodes[h].aa_node.probs, modelparams.hiddennodes[h2].aa_node.probs)
				transitionweights[h,h2] = -beta*aa_profile_similarity
				maxweight = max(maxweight, transitionweights[h,h2])
			end
			transitionweights[h,:] = transitionweights[h,:] .- maxweight
		end
		transitionweights = exp.(transitionweights)=#

		#=
		for h=1:modelparams.numhiddenstates
			transitionrate_counts[h,:] /= sum(transitionrate_counts[h,:])
		end=#

		for h=1:modelparams.numhiddenstates
			for h2=1:modelparams.numhiddenstates
				modelparams.transitionprobs[h,h2] = transitioncounts[h,h2]
			end
			modelparams.transitionprobs[h,:] = modelparams.transitionprobs[h,:] ./ sum(transitioncounts[h,:])

			

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

		if learnrates				
			for h=1:modelparams.numhiddenstates
				modelparams.transitionrates[h,h] = 0.0
				for h2=1:modelparams.numhiddenstates
					if h != h2
						modelparams.transitionrates[h,h2] = (transitionrate_counts[h,h2]+transitionrate_counts[h2,h])*sum(transitionrate_totals[h,:]) / 2.0 / length(trainingexamples)
						modelparams.transitionrates[h2,h] = (transitionrate_counts[h,h2]+transitionrate_counts[h2,h])*sum(transitionrate_totals[h2,:]) / 2.0 / length(trainingexamples)
						modelparams.transitionrates[h,h] -= modelparams.transitionrates[h,h2]
					end
				end
			end
			modelparams.mu = 1.0
			println("HQ", modelparams.transitionrates)	
			aatransitionrates = zeros(Float64, modelparams.alphabet, modelparams.alphabet)
			for aa1=1:modelparams.alphabet
				aatransitionrates[aa1,aa1] = 0.0
				for aa2=1:modelparams.alphabet
					if aa1 != aa2
						aatransitionrates[aa1,aa2] = (aatransitionrate_totals[aa1,aa2]+aatransitionrate_totals[aa2,aa1])  / 2.0 / length(trainingexamples)
						#aatransitionrates[aa1,aa2] = (aatransitionrate_counts[aa1,aa2]+aatransitionrate_counts[aa2,aa1])*sum(aatransitionrate_totals[aa1,:]) / 2.0 / length(trainingexamples)
						#aatransitionrates[aa2,aa1] = (aatransitionrate_counts[aa1,aa2]+aatransitionrate_counts[aa2,aa1])*sum(aatransitionrate_totals[aa2,:]) / 2.0 / length(trainingexamples)
						aatransitionrates[aa1,aa1] -= aatransitionrates[aa1,aa2]
					end
				end
			end
			println("AQ",aatransitionrates)		
			println("LG",[sum(LGexchangeability[aa,:]) for aa=1:20])
			modelparams.aa_exchangeablities = aatransitionrates*20.0
		end
		#println("probs ", modelparams.transitionprobs)
		#println("rates ", transitionrate_counts)
		#println("rates ", transitionrate_totals)


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

		fout = open(string("model", outputmodelname, ".model"), "w")
		Serialization.serialize(fout, modelparams)
		close(fout)
	end
end

function count_aminoacid_substitutions(rng::AbstractRNG, modelparams::ModelParams, node::TreeNode)
	aajumps = 0.0
	for aapath in node.data.aabranchpath.paths
		aajumps += length(aapath) - 1.0
	end
	return aajumps
end	

function infer()
	rng = MersenneTwister(10498012421321)
	Random.seed!(1234)

	#fin = open("model_h_25.model", "r")
	fin = open("model_h.15_learnrates.true.model", "r")	
	modelparams = Serialization.deserialize(fin)
	close(fin)

	samplebranchlengths = true

	
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
	
	#proteins,nodelist = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/alignments/HCV_REF_2014_ns5b_PRO_curated.fasta"))
	#proteins,nodelist = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/alignments/hiv-curated-sel.fasta"))=#

	#proteins,nodelist = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/influenza_a/HA/H1N1_selection2.fas"), abspath("../data/influenza_a/HA/H1N1_selection2_rooted.fas.nwk"))
	#proteins,nodelist = training_example_from_sequence_alignment(rng, modelparams, abspath("../data/influenza_a/HA/selection4.fasta"), abspath("../data/influenza_a/HA/selection4.fasta.nwk"))
	
	inputnodelist = deepcopy(nodelist)
	inputtreewriter = open("tree.input.nwk", "w")
	println(inputtreewriter, getnewick(nodelist[1]))
	close(inputtreewriter)

	countmatches = 1e-10
	counttotal = 1e-10

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
	subspersite_cache = Dict{Int,Array{Float64,1}}()
	for iter=1:10000
		if samplebranchlengths
			randindex = rand(2:length(nodelist))
			for node in nodelist
				if !isroot(node)	
					t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
					node.branchlength = t
					#=				
					if node.nodeindex == randindex
						println(iter," branches ", node.branchlength,"\t",t)				
						oldll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
						
						newll = augmentedloglikelihood(nodelist, Int[col for col=1:numcols], modelparams)
						println(iter," likeliho ", newll,"\t",oldll,"\t",propratio)
					end=#
				end
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
			felsensteinresample(rng, proteins, nodelist, col, cols, col, modelparams)
			felsensteinresample_aa(rng, proteins, nodelist, col, cols,col, modelparams)
		end

		selnode = nodelist[2]
		println(selnode.name)		
		alignedsequence = json_family["proteins"][selnode.seqindex]["aligned_sequence"]
		sequence, phi_psi, omega, bond_angles, bond_lengths = protein_to_lists(sampletreenode(rng, selnode, modelparams, alignedsequence))
	    inputsequence = replace(alignedsequence, "-" => "")
	    for (aa1,aa2) in zip(sequence,inputsequence)
	    	if aa1 == aa2
	    		countmatches += 1.0
	    	end
	    	counttotal += 1.0
	    	println(iter,"\t",aa1,"\t",aa2,"\t",countmatches/counttotal)
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

		consensustreewriter = open("tree.mean.consensus.nwk", "w")
		for outputnode in outputnodelist
			if !isroot(outputnode)
				outputnode.branchlength = mean(subspersite_cache[outputnode.nodeindex][max(1,div(length(subspersite_cache[outputnode.nodeindex]),2)):end])
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
function plot_nodes()
	fin = open("model_h_25_backup.model", "r")
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
