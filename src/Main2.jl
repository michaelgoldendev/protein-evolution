using Distributions
using LinearAlgebra

push!(LOAD_PATH,string(@__DIR__,"/../MolecularEvolution/src/"))
using MolecularEvolution

push!(LOAD_PATH,@__DIR__)
using BranchPaths
using CommonUtils
using LG
using Random
using CTMCs

include("Data.jl")

secondarystructure = "HBEGITSC"
aminoacids = "ACDEFGHIKLMNPQRSTVWY"

mutable struct HiddenNode
    phi_mu::Float64
    phi_kappa::Float64
    psi_mu::Float64
    psi_kappa::Float64
    aadist::Array{Float64,1}

    function HiddenNode()
        new(0.0,1.0,0.0,1.0, ones(Float64,20)/20.0)
    end
end

mutable struct ModelParams
    alphabet::Int
    aminoacidQ::Array{Float64,2}
    numhiddenstates::Int
    initialprobs::Array{Float64,1}
    transitionprobs::Array{Float64,2}
    transitionrates::Array{Float64,2}

    function ModelParams(aminoacidQ::Array{Float64,2}, numhiddenstates::Int)
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
        new(20*numhiddenstates,aminoacidQ,numhiddenstates,initialprobs,transitionprobs,transitionrates)
    end
end


mutable struct AugmentedNodeData <: NodeData
	branchpath::BranchPath
	dummy::Int
	Rmatrices::Array{Array{Float64,2},1}
	Pmatrices::Array{Array{Float64,2},1}
	vs::Array{Array{Float64,1},1}
	time::Array{Float64,1}

	AugmentedNodeData(branchpath::BranchPath, dummy::Int) = new(branchpath, dummy, Array{Float64,2}[], Array{Float64,2}[], Array{Float64,1}[], Float64[])
end

function gettransprobs(node::TreeNode, selcolin::Int, cols::Array{Int,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	branchiterator = BranchPathIterator(node.data.branchpath,cols)
	P = Matrix{Float64}(I, modelparams.alphabet, modelparams.alphabet)
	Pmatrices = Array{Float64,2}[]
	dummypath = Int[]
	dummytime = Float64[]
	for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
		prevprobs = ones(Float64,modelparams.numhiddenstates)
		if changecol == 2
			prevh = div(currstates[1]-1, 20) + 1
			prevprobs = modelparams.transitionprobs[prevh,:]
		end
		succprobs = ones(Float64,modelparams.numhiddenstates)
		if length(currstates) == 3
			succh = div(currstates[1]-1, 20) + 1
			succprobs = modelparams.transitionprobs[:,succh]
		end
		R = constructJointMatrix(modelparams, prevprobs, succprobs)
		dt = (currtime-prevtime)*node.branchlength
		Pi = exp(R*dt)
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
		prevprobs = ones(Float64,modelparams.numhiddenstates)
		if changecol == 2
			prevh = div(currstates[1]-1, 20) + 1
			prevprobs = modelparams.transitionprobs[prevh,:]
		end
		succprobs = ones(Float64,modelparams.numhiddenstates)
		if length(currstates) == 3
			succh = div(currstates[1]-1, 20) + 1
			succprobs = modelparams.transitionprobs[:,succh]
		end
		R = constructJointMatrix(modelparams, prevprobs, succprobs)
        #R = constructJointMatrix(modelparams, ones(Float64,modelparams.numhiddenstates), ones(Float64,modelparams.numhiddenstates))
		dt = (currtime-prevtime)*node.branchlength
    	Pi = exp(R*dt)
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


function felsensteinresample(rng::AbstractRNG, proteins::Array{Array{BranchState,1},1}, nodelist::Array{TreeNode,1}, selcolin::Int, cols::Array{Int,1}, modelparams::ModelParams)
	selcol = findfirst(x -> x == selcolin, cols)
	likelihoods = ones(Float64, length(nodelist), modelparams.alphabet)*-Inf
	logm = zeros(Float64,length(nodelist))

	stack = Int[1]
	while length(stack) > 0
		nodeindex = stack[end]
		node = nodelist[nodeindex]
		if isleafnode(node)
            v = observationlikelihood(proteins[node.seqindex], cols[selcol], modelparams)
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

        		likelihoods[nodeindex, :] = (lefttransprobs*likelihoods[leftchildindex,:]).*(righttransprobs*likelihoods[rightchildindex,:])

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
	for a=1:modelparams.alphabet
		rootseq[selcol] = a
		#logfreqs[a] = getInitialLogLikelihood(params, modelspecification, rootseq)
		logfreqs[a] = 1.0/modelparams.alphabet
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

function augmentedlikelihood(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, params::ModelParams)
	loglikelihood = 0.0
	for node in nodelist
		if !isroot(node)
			branchiterator = BranchPathIterator(node.data.branchpath,cols)

			for (prevstates,prevtime,currstates,currtime,changecol) in branchiterator
				Qii = 0.0
				len = length(cols)
				for selcol=1:len
					prevstate = prevstates[selcol]
					Qii += getratematrixrow(prevstates, selcol, params, modelspecification,componentindex)[prevstate]
				end
				dt = (currtime-prevtime)*node.branchlength
				if changecol > 0
					prevstate = prevstates[changecol]
					currstate = currstates[changecol]
					Qhi = getratematrixrow(prevstates, changecol, params, modelspecification,componentindex)[currstate]
					if isnan(Qhi) || Qhi <= 0.0
						println(Qhi)
						println("NANANANAA")
						println(params.Q)
						println(getratematrixrow(prevstates, changecol, params, modelspecification,componentindex))
						println(changecol,"\t",prevstate,"\t",currstate)
						println(node.data.branchpath.paths[changecol])
						println(node.data.branchpath.times[changecol])
					end
					loglikelihood += log(Qhi)
				end
				loglikelihood +=  Qii*dt
			end
		end
	end
	return loglikelihood
end

function constructJointMatrix(modelparams::ModelParams, prevprobs::Array{Float64,1}, succprobs::Array{Float64,1})
    aminoacidQ = modelparams.aminoacidQ
    Q = zeros(Float64, 20*modelparams.numhiddenstates, 20*modelparams.numhiddenstates)
    for h1=1:modelparams.numhiddenstates
        for h2=1:modelparams.numhiddenstates
            for aa1=1:20
                for aa2=1:20
                    i = (h1-1)*20 + aa1
                    j = (h2-1)*20 + aa2
                    if i != j
                        if h1 != h2 && aa1 != aa2
                            Q[i,j] = prevprobs[h2]*succprobs[h2]*aminoacidQ[aa1,aa2]*modelparams.transitionrates[h1,h2]
                        elseif h1 == h2 && aa1 != aa2
                            Q[i,j] = prevprobs[h2]*succprobs[h2]*aminoacidQ[aa1,aa2]
                        end
                        Q[i,i] -= Q[i,j]
                    end
                end
            end
        end
    end
    return Q
end

pair = loadpairs("data/etdbn_training_data.txt")[1]

function observationlikelihood(protein::Array{BranchState,1}, col::Int, modelparams::ModelParams)
    v = zeros(Float64, modelparams.alphabet)
    if protein[col].aa > 0
        for h=1:modelparams.numhiddenstates
            v[(h-1)*20 + protein[col].aa] = 1.0
        end
    else
        for h=1:modelparams.numhiddenstates
            for aa=1:20
                v[(h-1)*20 + aa]
            end
        end
    end
    return v
end

modelparams = ModelParams(LGmatrix, 2)
V = constructJointMatrix(modelparams, ones(Float64,modelparams.numhiddenstates), ones(Float64,modelparams.numhiddenstates))


rng = MersenneTwister(10498012421321)
Random.seed!(1234)
rootstate = rand(rng, 1:modelparams.alphabet)



root = TreeNode(1.0, "root")
root.nodeindex = 1
paths = Array{Int,1}[]
times = Array{Float64,1}[]
for i=1:length(pair[1])
	state = (rand(rng,1:modelparams.numhiddenstates)-1)*20 + max(1, pair[1][i].aa)
	push!(paths,Int[state,state])
	push!(times,Float64[0.0, 1.0])
end
root.data = AugmentedNodeData(BranchPath(paths,times), 1)


left = TreeNode(1.0, "left")
left.nodeindex = 2
left.seqindex = 1
paths = Array{Int,1}[]
times = Array{Float64,1}[]
for i=1:length(pair[1])
	p1_state = (rand(rng,1:modelparams.numhiddenstates)-1)*20 + max(1, pair[1][i].aa)
	p2_state = (rand(rng,1:modelparams.numhiddenstates)-1)*20 + max(1, pair[2][i].aa)
	path,time = modifiedrejectionsampling(rng, V*2.0, p1_state, p1_state, nothing)
	push!(paths,path)
	push!(times,time)
end

left.data = AugmentedNodeData(BranchPath(paths,times), 1)

right = TreeNode(1.0, "right")
right.nodeindex = 3
right.seqindex = 2
paths = Array{Int,1}[]
times = Array{Float64,1}[]
for i=1:length(pair[1])
	p1_state = (rand(rng,1:modelparams.numhiddenstates)-1)*20 + max(1, pair[1][i].aa)
	p2_state = (rand(rng,1:modelparams.numhiddenstates)-1)*20 + max(1, pair[2][i].aa)
	path,time = modifiedrejectionsampling(rng, V*2.0, p1_state, p2_state, nothing)
	push!(paths,path)
	push!(times,time)
end
right.data = AugmentedNodeData(BranchPath(paths,times), 1)

addchild(root, left)
addchild(root, right)
nodelist = TreeNode[]
push!(nodelist,root)
push!(nodelist,left)
push!(nodelist,right)

initialprobs = rand(Dirichlet(modelparams.alphabet,1.0))

proteins = Array{BranchState,1}[]
push!(proteins,pair[1])
push!(proteins,pair[2])

numcols = length(proteins[1])
for col=1:numcols
	println(col)
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
	felsensteinresample(rng, proteins, nodelist, col, cols, modelparams)
end
