using Distributions
using LinearAlgebra
using Serialization
using JSON

push!(LOAD_PATH,string(@__DIR__,"/../../../dev/MolecularEvolution/src/"))
using MolecularEvolution

push!(LOAD_PATH,string(@__DIR__))
using EMNodes
using PDBBuilder
push!(LOAD_PATH,@__DIR__)
using BranchPaths
using CommonUtils
using LG
using Random
using CTMCs

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

function siteloglikelihood(site::BranchState, h::Int)
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
	if site.bondlength1 > -100.0 && site.bondlength2 > -100.0 && site.bondlength3 > -100.0
		ll += logpdf(modelparams.hiddennodes[h].bond_lengths_node.mvn, Float64[site.bondlength1, site.bondlength2, site.bondlength3])
	end=#
	return ll
end

function observationlikelihood(protein::Array{BranchState,1}, col::Int, modelparams::ModelParams)
    v = zeros(Float64, modelparams.alphabet)
	for h=1:modelparams.numhiddenstates
		v[h] = siteloglikelihood(protein[col], h)
	end
	return exp.(v .- maximum(v))
end


function observationloglikelihood(proteins::Array{Array{BranchState,1},1}, nodelist::Array{TreeNode,1}, modelparams::ModelParams)
	ll = 0.0
	for node in nodelist
		if isleafnode(node)
			for col=1:length(node.data.branchpath.paths)
				h = node.data.branchpath.paths[col][end]				 
				ll += siteloglikelihood(proteins[node.seqindex][col], h)
			end
		end
	end
	return ll
end

function sampletip(rng::AbstractRNG, node::TreeNode, modelparams::ModelParams)
	numcols = length(node.data.branchpath.paths)
	dihedralangles = Tuple{Float64,Float64,Float64}[]
	for col=1:numcols
		phiz = vonmisesrand(rng, modelparams.hiddennodes[node.data.branchpath.paths[col][end]].phi_node.dist)
		psiz = vonmisesrand(rng, modelparams.hiddennodes[node.data.branchpath.paths[col][end]].psi_node.dist)		
		omegaz = vonmisesrand(rng, modelparams.hiddennodes[node.data.branchpath.paths[col][end]].omega_node.dist)
		push!(dihedralangles, (phiz,psiz,omegaz))
	end
	return dihedralangles
end

function sampletip2(rng::AbstractRNG, node::TreeNode, modelparams::ModelParams, aligned_sequence::String)
	numcols = length(node.data.branchpath.paths)
	sequence = ""
    phi_psi = Tuple{Float64,Float64}[] 
    omega = Float64[]
    bond_angles = Tuple{Float64,Float64,Float64}[]
    bond_lengths = Array{Float64,1}[]
	for col=1:numcols
		aa = aligned_sequence[col]
		if aa != '-'
			sequence = string(sequence, aa)
			h = node.data.branchpath.paths[col][end]
			#println("A")
			phiz = pimod(vonmisesrand(rng, modelparams.hiddennodes[h].phi_node.dist))
			#println("B")
			psiz = pimod(vonmisesrand(rng, modelparams.hiddennodes[h].psi_node.dist))
			push!(phi_psi, (phiz,psiz))
			#println("C")
			omegaz = pimod(vonmisesrand(rng, modelparams.hiddennodes[h].omega_node.dist))
			push!(omega,omegaz)
			#println("D")
			bond_angle1z = pimod(vonmisesrand(rng, modelparams.hiddennodes[h].bond_angle1_node.dist))
			#println("E")
			bond_angle2z = pimod(vonmisesrand(rng, modelparams.hiddennodes[h].bond_angle2_node.dist))
			#println("F")
			bond_angle3z = pimod(vonmisesrand(rng, modelparams.hiddennodes[h].bond_angle3_node.dist))
			push!(bond_angles, (bond_angle1z, bond_angle2z, bond_angle3z))
			#println("G")
			push!(bond_lengths, rand(rng, modelparams.hiddennodes[h].bond_lengths_node.mvn))
			#println("H")
		end
	end
	return sequence, phi_psi, omega, bond_angles, bond_lengths
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

function training_example_from_json_family(rng::AbstractRNG, modelparams, json_family)
	root = gettreefromnewick(json_family["newick_tree"])
	binarize!(root)
	nodelist = getnodelist(root)
	seqindex = 1
	for (index,node) in enumerate(nodelist)
		if isleafnode(node)
			node.seqindex = seqindex
			seqindex += 1
		end
		node.nodeindex = index
	end

	V = constructJointMatrix(modelparams, ones(Float64,modelparams.numhiddenstates), ones(Float64,modelparams.numhiddenstates))

	numcols = length(json_family["proteins"][1]["aligned_sequence"])

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

	proteins = Array{BranchState,1}[]
	for p=1:length(json_family["proteins"])
	    protein = BranchState[]
	    json_protein = json_family["proteins"][p]
	    for (aa,phi_psi,omega,bond_lengths,bond_angles) in zip(json_protein["aligned_sequence"], json_protein["aligned_phi_psi"], json_protein["aligned_omega"], json_protein["aligned_bond_lengths"], json_protein["aligned_bond_angles"])
	    	phi = phi_psi[1]
	    	psi = phi_psi[2]
	        push!(protein, BranchState(0,indexof(string(aa), aminoacids),phi,omega,psi,bond_lengths[1],bond_lengths[2],bond_lengths[3],bond_angles[1],bond_angles[2],bond_angles[3]))
	    end
	    push!(proteins, protein)
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

rng = MersenneTwister(10498012421321)
Random.seed!(1234)

family_files = filter(f -> endswith(f,".fam"), readdir("../data/families/"))

numhiddenstates = 5
modelparams = ModelParams(LGmatrix,numhiddenstates,0.2)

trainingexamples = Tuple[]
for family_file in family_files[1:30]
	full_path = abspath(joinpath("../data/families/", family_file))
	json_family = JSON.parse(open(full_path, "r"))
	if 2 <= length(json_family["proteins"]) <= 1e10
		training_example = training_example_from_json_family(rng, modelparams, json_family)
		push!(trainingexamples, training_example)
		println(json_family["newick_tree"])
	end
end
#=
fin = open("model_h_15.model", "r")
modelparams = Serialization.deserialize(fin)
close(fin)=#




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
					site = proteins[node.seqindex][col]
					if site.aa > 0
						modelparams.hiddennodes[h].aa_node.counts[site.aa] += 1.0
						push!(modelparams.hiddennodes[h].phi_node.data, site.phi)
						push!(modelparams.hiddennodes[h].omega_node.data, site.omega)
						push!(modelparams.hiddennodes[h].psi_node.data, site.psi)
						push!(modelparams.hiddennodes[h].bond_angle1_node.data, site.bond_angle1)
						push!(modelparams.hiddennodes[h].bond_angle2_node.data, site.bond_angle2)
						push!(modelparams.hiddennodes[h].bond_angle3_node.data, site.bond_angle3)
						add_point(modelparams.hiddennodes[h].bond_lengths_node, Float64[site.bondlength1, site.bondlength2, site.bondlength3])
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

		#=
		for aa=1:20
			modelparams.hiddennodes[h].aadist[aa] = aacounts[h,aa]
		end
		modelparams.hiddennodes[h].aadist /= sum(modelparams.hiddennodes[h].aadist)
		=#
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
    sequence, phi_psi, omega, bond_angles, bond_lengths = sampletip2(rng, nodelist[2], modelparams, json_family["proteins"][1]["aligned_sequence"])

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


	#save(string("model_h_",modelparams.numhiddenstates,"_",iter,".model"), "modelparams", modelparams)
	fout = open(string("model_h_",modelparams.numhiddenstates, ".model"), "w")
	Serialization.serialize(fout, modelparams)
	close(fout)



	#println(aacounts./sum(aacounts,2))
end
