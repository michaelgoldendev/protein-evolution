using Distributions
using LinearAlgebra
using Serialization
push!(LOAD_PATH,string(@__DIR__,"/../../MolecularEvolution/src/"))
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

function pimod(angle::Float64)
  theta = mod2pi(angle)
  if theta > pi
    return theta -2.0*pi
  else
    return theta
  end
end

mutable struct VonMisesNode
	rx::Float64
	ry::Float64
	N::Float64
	mu::Float64
	kappa::Float64
	dist::VonMises
	data::Array{Float64,1}

	function VonMisesNode()
		mu = 0.0
		kappa = 1.0
		new(0.0, 0.0, 0.0, mu, kappa, VonMises(mu, kappa), Float64[])
	end
end

function estimatevonmises(vonmises_node::VonMisesNode)
	if length(vonmises_node.data) > 0
		vonmises_node.rx = 0.0
		vonmises_node.ry = 0.0
		for theta in vonmises_node.data
			vonmises_node.rx += cos(theta)
			vonmises_node.ry += sin(theta)
		end
		vonmises_node.N = length(vonmises_node.data)
		vonmises_node.data = Float64[]
	end

	rx = vonmises_node.rx
	ry = vonmises_node.ry
	n = vonmises_node.N

	c = rx / n
	s = ry / n
	rho = sqrt(c*c + s*s)

	if s > 0
		mu = acos(c / rho)
	else
		mu = 2.0*pi - acos(c / rho)
	end


	if rho < 2.0/3.0
		kappa = rho * ((2.0 - rho*rho) / (1.0 - rho*rho))
	else
		kappa = (rho + 1.0) / (4.0 * rho * (1 - rho))
	end

	vonmises_node.mu = mu
	vonmises_node.kappa = kappa
	vonmises_node.dist = VonMises(mu,kappa)
end

mutable struct HiddenNode
	phi_node::VonMisesNode
	psi_node::VonMisesNode
	omega_node::VonMisesNode
    aadist::Array{Float64,1}

    function HiddenNode()
        new(VonMisesNode(),VonMisesNode(),VonMisesNode(), ones(Float64,20)/20.0)
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
        new(numhiddenstates,aminoacidQ,numhiddenstates,initialprobs,transitionprobs,transitionrates,hiddennodes,mu)
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
			prevh = currstates[1]
			prevprobs = modelparams.transitionprobs[prevh,:]
		end
		succprobs = ones(Float64,modelparams.numhiddenstates)
		if length(currstates) == 3
			succh = currstates[3]
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
	#selcol = findfirst(x -> x == selcolin, cols)
	selcol = selcolin
	likelihoods = ones(Float64, length(nodelist), modelparams.alphabet)*-Inf
	logm = zeros(Float64,length(nodelist))

	stack = Int[1]
	while length(stack) > 0
		nodeindex = stack[end]
		node = nodelist[nodeindex]
		if isleafnode(node)
            v = observationlikelihood(proteins[node.seqindex], selcolin, modelparams)
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
		logfreqs[a] = log(prevprobs[a]*nextprobs[a])
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

function observationloglikelihood(proteins::Array{Array{BranchState,1},1}, nodelist::Array{TreeNode,1}, modelparams::ModelParams)
	ll = 0.0
	for node in nodelist
		if isleafnode(node)
			for col=1:length(node.data.branchpath.paths)
				h = node.data.branchpath.paths[col][end]
				site = proteins[node.seqindex][col]
				if site.aa > 0
					ll += log(modelparams.hiddennodes[h].aadist[site.aa])
				end
				if site.phi > -100.0
					ll += logpdf(modelparams.hiddennodes[h].phi_node.dist, site.phi)
				end
				if site.psi > -100.0
					ll += logpdf(modelparams.hiddennodes[h].psi_node.dist, site.psi)
				end
			end
		end
	end
	return ll
end

function augmentedlikelihood(nodelist::Array{TreeNode,1}, cols::Array{Int,1}, modelparams::ModelParams)
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

function observationlikelihood(protein::Array{BranchState,1}, col::Int, modelparams::ModelParams)
    v = zeros(Float64, modelparams.alphabet)
	for h=1:modelparams.numhiddenstates
		if protein[col].aa > 0
			v[h] = modelparams.hiddennodes[h].aadist[protein[col].aa]
		else
			for aa=1:20
				v[h] = modelparams.hiddennodes[h].aadist[aa]
			end
		end

		if protein[col].phi > -100.0
			v[h] *= pdf(modelparams.hiddennodes[h].phi_node.dist, protein[col].phi)
		end
		if protein[col].psi > -100.0
			v[h] *= pdf(modelparams.hiddennodes[h].psi_node.dist, protein[col].psi)
		end
	end
	return v
end


function createtrainingexample(pair)
	V = constructJointMatrix(modelparams, ones(Float64,modelparams.numhiddenstates), ones(Float64,modelparams.numhiddenstates))

	root = TreeNode(1.0, "root")
	root.nodeindex = 1
	paths = Array{Int,1}[]
	times = Array{Float64,1}[]
	for i=1:length(pair[1])
		state = rand(rng,1:modelparams.numhiddenstates)
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
		p1_state = rand(rng,1:modelparams.numhiddenstates)
		p2_state = rand(rng,1:modelparams.numhiddenstates)
		path,time = modifiedrejectionsampling(rng, V*0.5, p1_state, p1_state, nothing)
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
		p1_state = rand(rng,1:modelparams.numhiddenstates)
		p2_state = rand(rng,1:modelparams.numhiddenstates)
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

	proteins = Array{BranchState,1}[]
	push!(proteins,pair[1])
	push!(proteins,pair[2])

	return (proteins, nodelist)
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

function vonmisesrand(rng::AbstractRNG, vonmises::VonMises)
	dist = VonMises(0.0, vonmises.κ)
	u = rand(rng)
	v = -1.0
	lower = -pi
	upper = pi
	mid = (upper-lower)/2.0 + lower
	for i=1:10000
		mid = (upper-lower)/2.0 + lower
		v = cdf(dist, mid)
		if abs(v-u) < 1e-10
			break
		elseif v < u
			lower = mid
		else
			upper = mid
		end
	end
	if abs(v-u) > 0.001
		println("sample ",v,"\t",u,"\t",cdf(vonmises,mid),"\t",vonmises,"\t",lower,"\t",upper)
	end
	return mod2pi(mid + vonmises.μ)
end

function sampletip(rng::AbstractRNG, node::TreeNode, modelparams::ModelParams)
	numcols = length(node.data.branchpath.paths)
	dihedralangles = Tuple{Float64,Float64}[]
	for col=1:numcols
		phiz = vonmisesrand(rng, modelparams.hiddennodes[node.data.branchpath.paths[col][end]].phi_node.dist)
		psiz = vonmisesrand(rng, modelparams.hiddennodes[node.data.branchpath.paths[col][end]].psi_node.dist)
		push!(dihedralangles, (phiz,psiz))
	end
	return dihedralangles
end

pairs = loadpairs("../data/etdbn_training_data.txt")

numhiddenstates = 10
modelparams = ModelParams(LGmatrix,numhiddenstates,0.2)

#=
fin = open("model_h_15.model", "r")
modelparams = Serialization.deserialize(fin)
close(fin)=#

rng = MersenneTwister(10498012421321)
Random.seed!(1234)
rootstate = rand(rng, 1:modelparams.alphabet)

trainingexamples = Tuple[]
for p=1:600
	push!(trainingexamples, createtrainingexample(pairs[p]))
end

logwriter = open(string("trace", modelparams.numhiddenstates,".log"), "w")
println(logwriter, "iter","\t","loglikelihood")
for iter=1:1000
	aacounts = ones(Float64, modelparams.numhiddenstates, 20)*0.1
	transitioncounts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
	transitionratecounts = ones(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)*0.1
	transitionratetotals = zeros(Float64, modelparams.numhiddenstates, modelparams.numhiddenstates)
	for h=1:modelparams.numhiddenstates
		transitionratecounts[h,h] = 0.0
	end

	for (proteins,nodelist) in trainingexamples
		numcols = length(proteins[1])
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

					#=
					dt = (currtime-prevtime)*node.branchlength
					println(prevstates[1], "\t", prevstates[2], "\t", dt
					=#
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

		for col=1:numcols
			for node in nodelist
				for p=1:length(node.data.branchpath.paths[col])
					h = node.data.branchpath.paths[col][p]
					if p == length(node.data.branchpath.paths[col]) && isleafnode(node)
						if proteins[node.seqindex][col].aa > 0
							aacounts[h,proteins[node.seqindex][col].aa] += 1.0
							push!(modelparams.hiddennodes[h].phi_node.data, proteins[node.seqindex][col].phi)
							push!(modelparams.hiddennodes[h].psi_node.data, proteins[node.seqindex][col].psi)
						end
					else
						sampleaa = CommonUtils.sample(rng, modelparams.hiddennodes[h].aadist)
						aacounts[h,sampleaa] += 1.0
					end
				end
			end
		end


		for node in nodelist
			if !isroot(node)
				t,propratio = proposebranchlength(rng, node, Int[col for col=1:numcols], modelparams)
				node.branchlength = t
				#println(t,"\t",propratio)
			end
		end
	end

	totalll = 0.0
	for (proteins,nodelist) in trainingexamples
		numcols = length(proteins[1])
		#ll = augmentedlikelihood(nodelist, Int[col for col=1:numcols], modelparams)
		ll = observationloglikelihood(proteins, nodelist, modelparams)
		totalll += ll
	end
	println(logwriter, iter,"\t",totalll)
	flush(logwriter)

	for h=1:modelparams.numhiddenstates
		for aa=1:20
			modelparams.hiddennodes[h].aadist[aa] = aacounts[h,aa]
		end
		modelparams.hiddennodes[h].aadist /= sum(modelparams.hiddennodes[h].aadist)

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

		estimatevonmises(modelparams.hiddennodes[h].phi_node)
		estimatevonmises(modelparams.hiddennodes[h].psi_node)
		println(iter,"\t",h,"\t",modelparams.hiddennodes[h].phi_node.mu,"\t",modelparams.hiddennodes[h].phi_node.kappa)
		println(iter,"\t",h,"\t",modelparams.hiddennodes[h].psi_node.mu,"\t",modelparams.hiddennodes[h].psi_node.kappa)

		for aa=1:20
			println(iter,"\t",h,"\t",aminoacids[aa],"\t",modelparams.hiddennodes[h].aadist[aa])
		end
	end
	println("probs ", modelparams.transitionprobs)
	println("rates ", transitionratecounts)
	println("rates ", transitionratetotals)

	proteins,nodelist = trainingexamples[1]
	phipsisample = sampletip(rng, nodelist[2], modelparams)
	for col=1:length(phipsisample)
		node = nodelist[2]
		h = node.data.branchpath.paths[col][end]
		println(pimod(phipsisample[col][1]),"\t", proteins[1][col].phi,"\t",modelparams.hiddennodes[h].phi_node.dist,"\t",pimod(phipsisample[col][2]),"\t", proteins[1][col].psi,"\t",modelparams.hiddennodes[h].psi_node.dist)
	end
	#save(string("model_h_",modelparams.numhiddenstates,"_",iter,".model"), "modelparams", modelparams)
	fout = open(string("model_h_",modelparams.numhiddenstates, ".model"), "w")
	Serialization.serialize(fout, modelparams)
	close(fout)



	#println(aacounts./sum(aacounts,2))
end
