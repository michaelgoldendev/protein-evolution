module CTMCs
	using LinearAlgebra
	using Random

	push!(LOAD_PATH,@__DIR__)
	using CommonUtils

	push!(LOAD_PATH,string(@__DIR__,"/../../MolecularEvolution/src/"))
	using MolecularEvolution

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

	function recursivesampling(rng::AbstractRNG, Q::Array{Float64,2}, a::Int, b::Int, t::Float64=1.0)
		decomposition = eigen(Q)
	  	D, V = decomposition.values, decomposition.vectors
	  	Vi = inv(V)
	  	return recursivesampling(rng, Q, a, b, complex(V), complex(D), complex(Vi), t)
	end

	Pcache = Array{Float64,2}[zeros(Float64,1,1) for level=1:40]
	function recursivesampling(rng::AbstractRNG, Q::Array{Float64,2}, a::Int, b::Int, V::Array{Complex{Float64},2}, D::Array{Complex{Float64},1}, Vi::Array{Complex{Float64},2}, t::Float64=1.0)		
	  	for level=1:40
	  		Pcache[level][1,1] = 0.0
	  	end 
	  	paths,times,count = recursivesampling(Pcache, rng, Q, V, D, Vi, a, b, t)
	  	#=
	  	if count >= 400
	  		println("RECURSIONS: ", count)
	  		println("T",t)
	  		println(paths,"\t",times)
	  	end=#

	  	return paths,times
	end

	function recursivesampling(Pcache::Array{Array{Float64,2},1}, rng::AbstractRNG, Q::Array{Float64,2}, V::Array{Complex{Float64},2}, D::Array{Complex{Float64},1}, Vi::Array{Complex{Float64},2}, a::Int, b::Int, t::Float64, startt::Float64=0.0, endt::Float64=1.0, level::Int=1)
		levelstop = 40
		delta = (endt-startt)/2.0
		midt = startt + delta
		count = 1
		#=
		P = exp(Q*delta)
		v = P[a,:].*P[:,b]
		m = sample(rng, v)
		=#

		#v = V[a,:].*exp.(D*delta).*Vi[:,b]
		if Pcache[level][1,1] == 0.0
			Pcache[level] = absmat(real(V*Diagonal(exp.(D*t*delta))*Vi))
		end
		P = Pcache[level]
		v = P[a,:].*P[:,b]
		m = sample(rng, v)

		cont = true
		if a == m
			samplet =  log(1.0 - rand(rng))/(Q[m,m]*t*delta)
			if samplet > delta
				cont = false
			end
		end
		subpaths1 = Int[a]
		subtimes1 = Float64[startt]
		if cont && level < levelstop
			subpaths1, subtimes1, count1 = recursivesampling(Pcache,rng, Q, V, D, Vi, a, m, t, startt, midt, level+1)
			count += count1
		end

		cont = true
		if b == m
			samplet =  log(1.0 - rand(rng))/(Q[m,m]*t*delta)
			if samplet > delta
				cont = false
			end
		end
		count += 1
		subpaths2 = Int[m]
		subtimes2 = Float64[midt]
		if cont  && level < levelstop
			subpaths2, subtimes2,count2 = recursivesampling(Pcache,rng, Q, V, D, Vi, m, b, t, midt, endt, level+1)
			count += count2
		end

		paths = Int[subpaths1[1]]
		times = Float64[subtimes1[1]]
		for i=2:length(subpaths1)
			if subpaths1[i] != paths[end]
				push!(paths, subpaths1[i])
				push!(times, subtimes1[i])
			end
		end
		for i=1:length(subpaths2)
			if subpaths2[i] != paths[end]
				push!(paths, subpaths2[i])
				push!(times, subtimes2[i])
			end
		end

		#=
		if count >= 100000
			println("P",P)
			println("delta:",delta,"\t",t)
			println("V",v)
		end=#
		return paths,times,count
	end

	function approximatesampling(rng::AbstractRNG, Q::Array{Float64,2}, t::Float64, a::Int, b::Int, divisions=100000)
		alphabet = size(Q,1)
		P = exp(Q*(t/divisions))
		logP = log.(P)
		forward = ones(Float64,divisions,alphabet)*-Inf
		forward[1,a] = 0.0
		for i=2:divisions-1
			for prevc=1:alphabet
				for c=1:alphabet
					forward[i,c] = logsumexp(forward[i,c], forward[i-1,prevc]+logP[prevc,c])
				end
			end
		end
		for i=1:divisions
			forward[i,:] = exp.(forward[i,:].-maximum(forward[i,:]))
		end

		states = zeros(Int,divisions)
		states[end] = b
		lastc = b
		for i=divisions-1:-1:1
			probs = forward[i,:].*P[:,lastc]
			lastc = sample(rng, probs)
			states[i] = lastc
		end
		indices = Int[]
		lastc = 0
		for i=1:divisions
			if states[i] != lastc
				lastc = states[i]
				push!(indices,i)
			end
		end
		path = Int[states[1]]
		time = Float64[0.0]
		for i=2:length(indices)
			timestart = (indices[i]-1.0)/divisions
			timeend = indices[i]/divisions
			delta = timeend-timestart
			push!(path, states[indices[i]])
			push!(time, timestart+rand(rng)*delta)
		end
		return path,time
	end

	

	basepairings = Int[0,0,0,2,0,0,1,0,0,1,0,3,2,0,3,0]

	export MatrixExponential
	mutable struct MatrixExponential
		Q::Array{Float64,2}
		D::Array{Float64,1}
		V::Array{Float64,2}
		Vi::Array{Float64,2}
		Pts::Array{Array{Float64,2},1}

		function MatrixExponential(Q::Array{Float64,2}, nodelist::Array{TreeNode,1})
			decomposition = eigen(Q)
		  	D, V = decomposition.values, decomposition.vectors
		  	Vi = inv(V)
			Pts = Array{Float64,2}[]
			for node in nodelist
				if isroot(node)
					push!(Pts, eye(size(Q,1)))
				else
					push!(Pts, expm(node.branchlength*Q))
					#push!(Pts, V*Diagonal(exp(D*node.branchlength))*Vi)
				end
			end
			new(Q,D,V,Vi,Pts)
		end
	end

	export CTMC
	mutable struct CTMC
		eqfreqs::Array{Float64, 1}
		logeqfreqs::Array{Float64, 1}
		S::Array{Float64,2}
		Q::Array{Float64,2}
		D::Array{Float64,1}
		V::Array{Float64,2}
		Vi::Array{Float64,2}
		t::Float64
		Pt::Array{Float64,2}
		logPt::Array{Float64,2}
		enabled::Bool

		function CTMC(eqfreqs::Array{Float64,1}, S::Array{Float64,2}, t::Float64)
			n = length(eqfreqs)
			#Q = Diagonal(eqfreqs)*S

			Q = zeros(Float64, n, n)
			for i=1:n
				for j=1:n
					Q[i,j] = S[i,j]*eqfreqs[j]
				end
			end
			for i=1:n
				Q[i,i] = 0.0
				for j=1:n
					if i != j
						Q[i,i] -= Q[i,j]
					end
				end
			end
			decomposition = eigen(Q)
		  	D, V = decomposition.values, decomposition.vectors
		  	Vi = inv(V)
			Pt = V*Diagonal(exp(D*t))*Vi
			logPt = log(Pt)
			return new(eqfreqs, log(eqfreqs), S, Q, D,V, Vi, t, Pt, logPt, true)
		end
	end

	export gtr
	function gtr(q1, q2, q3, q4, q5, q6, piv::Array{Float64,1})
		S = zeros(Float64,4,4)
		S[1,2] = q1
		S[1,3] = q2
		S[1,4] = q3

		S[2,1] = q1
		S[2,3] = q4
		S[2,4] = q5

		S[3,1] = q2
		S[3,2] = q4
		S[3,4] = q6

		S[4,1] = q3
		S[4,2] = q5
		S[4,3] = q6

		for i=1:size(S,1)
			for j=1:size(S,2)
				S[i,j] *= piv[j]
			end
		end

		Q = copy(S)
		for i=1:size(Q,1)
			diagelem = 0.0
			for j=1:size(Q,2)
				if i != j
					diagelem -= S[i,j]
				end
			end
			Q[i,i] = diagelem
		end
		return Q
	end

	export calculatemusefreqs
	function calculatemusefreqs(obsFreqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64)
		dinucfreqs = zeros(Float64,16)

		piGC = obsFreqs[2]*obsFreqs[3]
		piAT = obsFreqs[1]*obsFreqs[4]
		piGT = obsFreqs[3]*obsFreqs[4]

		kappa = 1.0/(1.0 + 2.0*(piAT*(lambdaAT*lambdaAT-1.0)) + 2.0*(piGC*(lambdaGC*lambdaGC-1.0)) + 2.0*(piGT*(lambdaGT*lambdaGT-1.0)))
		for h=1:4
			for v=1:4
				idx = (h-1)*4+v
				if basepairings[idx] == 1
					dinucfreqs[idx] = kappa*lambdaGC*lambdaGC*obsFreqs[h]*obsFreqs[v]
				elseif basepairings[idx] == 2
					dinucfreqs[idx] = kappa*lambdaAT*lambdaAT*obsFreqs[h]*obsFreqs[v]
				elseif basepairings[idx] == 3
					dinucfreqs[idx] = kappa*lambdaGT*lambdaGT*obsFreqs[h]*obsFreqs[v]
				else
					dinucfreqs[idx] = kappa*obsFreqs[h]*obsFreqs[v]
				end
			end
		end
		return dinucfreqs
	end

	export logproductmusefreqs
	function logproductmusefreqs(obsFreqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64, dinuccounts::Array{Float64,1})
		piGC = obsFreqs[2]*obsFreqs[3]
		piAT = obsFreqs[1]*obsFreqs[4]
		piGT = obsFreqs[3]*obsFreqs[4]

		kappa = 1.0/(1.0 + 2.0*(piAT*(lambdaAT*lambdaAT-1.0)) + 2.0*(piGC*(lambdaGC*lambdaGC-1.0)) + 2.0*(piGT*(lambdaGT*lambdaGT-1.0)))
		ll = 0.0
		for h=1:4
			for v=1:4
				idx = (h-1)*4+v
				if basepairings[idx] == 1
					ll += dinuccounts[idx]*log(kappa*lambdaGC*lambdaGC*obsFreqs[h]*obsFreqs[v])
				elseif basepairings[idx] == 2
					ll += dinuccounts[idx]*log(kappa*lambdaAT*lambdaAT*obsFreqs[h]*obsFreqs[v])
				elseif basepairings[idx] == 3
					ll += dinuccounts[idx]*log(kappa*lambdaGT*lambdaGT*obsFreqs[h]*obsFreqs[v])
				else
					ll += dinuccounts[idx]*log(kappa*obsFreqs[h]*obsFreqs[v])
				end
			end
		end
		return ll
	end

	export muse
	"""
	Returns a \$16 \\times 16\$ matrix.
	"""
	function muse(lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64,q1::Float64, q2::Float64, q3::Float64, q4::Float64, q5::Float64, q6::Float64, obsfreqs::Array{Float64,1}=Float64[0.25,0.25,0.25,0.25], siteRate1::Float64=1.0, siteRate2::Float64=1.0)
		gtrQ = gtr(q1,q2,q3,q4,q5,q6,ones(Float64,4))
		dinucfreqs = calculatemusefreqs(obsfreqs,lambdaGC,lambdaAT,lambdaGT)
		Q = zeros(Float64, 16, 16)
		for h=1:16
			for v=h+1:16
				if v != h
					fromNuc = 0
					toNuc   = 0
					if div(h-1,4)+1 == div(v-1,4)+1
						toNuc   = ((v-1) % 4) + 1
						fromNuc = ((h-1) % 4) + 1
					elseif ((v-1)%4) + 1 == ((h-1) %4)+1
						toNuc   = div(v-1,4) + 1
						fromNuc = div(h-1,4) + 1
					end
					if fromNuc > 0
						rateMult  = 1.0
						rateMult2 = 1.0

						if basepairings[v] == 1
							rateMult  *= lambdaGC
							rateMult2 /= lambdaGC
						end
						if basepairings[h] == 1
							rateMult /= lambdaGC
							rateMult2 *= lambdaGC
						end


						if basepairings[v] == 2
							rateMult  *= lambdaAT
							rateMult2 /= lambdaAT
						end
						if basepairings[h] == 2
							rateMult /= lambdaAT
							rateMult2 *= lambdaAT
						end

						if basepairings[v] == 3
							rateMult  *= lambdaGT
							rateMult2 /= lambdaGT
						end
						if basepairings[h] == 3
							rateMult /= lambdaGT
							rateMult2 *= lambdaGT
						end

						if fromNuc != toNuc
							rateMult  *= gtrQ[fromNuc,toNuc]
							rateMult2 *= gtrQ[fromNuc,toNuc]
						end

						if div(h-1,4) == div(v-1,4)
							rateMult  *= siteRate2
							rateMult2 *= siteRate2
						elseif (v-1) %4  == (h-1) %4
							rateMult  *= siteRate1
							rateMult2 *= siteRate1
						end

						Q[h,v] = rateMult*obsfreqs[toNuc]
						Q[v,h] = rateMult2*obsfreqs[fromNuc]
					end
				end
			end
		end

		for i=1:16
			Q[i,i] = 0.0
			for j=1:16
				if i != j
					Q[i,i] -= Q[i,j]
				end
			end
		end

		return Q
	end

	export MuseModel
	mutable struct MuseModel
		alphabet::Int
		obsfreqs::Array{Float64,1}
		freqs::Array{Float64,1}
		logfreqs::Array{Float64,1}
		Q::Array{Float64,2}
		lambdaGC::Float64
		lambdaAT::Float64
		lambdaGT::Float64
		q1::Float64
		q2::Float64
		q3::Float64
		q4::Float64
		q5::Float64
		q6::Float64
		siteRate1::Float64
		siteRate2::Float64

		function MuseModel(obsfreqs::Array{Float64,1}, lambdaGC::Float64, lambdaAT::Float64, lambdaGT::Float64, q1::Float64, q2::Float64, q3::Float64, q4::Float64, q5::Float64, q6::Float64, siteRate1::Float64, siteRate2::Float64)
			freqs = calculatemusefreqs(obsfreqs,lambdaGC,lambdaAT,lambdaGT)
			logfreqs = log(freqs)
			Q =  muse(lambdaGC, lambdaAT, lambdaGT, q1, q2, q3, q4, q5, q6, obsfreqs, siteRate1, siteRate2)
			new(16, obsfreqs, freqs, logfreqs, Q, lambdaGC, lambdaAT, lambdaGT, q1, q2, q3, q4, q5, q6, siteRate1, siteRate2)
		end
	end

	export PathSampleGenerator
	mutable struct PathSampleGenerator
		rng::AbstractRNG
		alphabet::Int
		numnodes::Int
		paths::Array{Array{Array{Int,1},1},1}
		times::Array{Array{Array{Float64,1},1},1}
		Q::Array{Float64,2}
		maxcachesize::Int

		function PathSampleGenerator(rng::AbstractRNG, nodelist::Array{TreeNode,1}, Q::Array{Float64,2}, maxcachesize::Int=20)
			alphabet = size(Q,1)
			paths = Array{Array{Int,1},1}[]
			times = Array{Array{Float64,1},1}[]
			for a=1:alphabet
				for b=1:alphabet
					for node in nodelist
						push!(paths, Array{Int,1}[])
						push!(times, Array{Float64,1}[])
					end
				end
			end
			new(rng, alphabet, length(nodelist), paths, times, Q, maxcachesize)
		end
	end

	export modifiedrejectionsampling2
	function  modifiedrejectionsampling2(generator::PathSampleGenerator, node::TreeNode, x0::Int, xt::Int)
		Q = generator.Q*node.branchlength
		rng = generator.rng
		path = Int[]
		times = Float64[]
		R = copy(Q)
		for i=1:generator.alphabet
			R[i,i] = 0.0
		end
		S = copy(R)
		for i=1:size(Q,1)
			s = 0.0
			for j=1:size(Q,2)
				s += S[i,j]
			end
			for j=1:size(Q,2)
				S[i,j] /= s
			end
		end

		count = 0
		while true
			path = Int[x0]
			times = Float64[0.0]
			totalt = 0.0
			if x0 != xt
				r = rand(rng)
				totalt += log(1.0-r*(1.0-exp(Q[path[end],path[end]])))/Q[path[end],path[end]]
				push!(path, sample(rng, S[path[end],:]))
				push!(times, totalt)
			end

			while true
				r = rand(rng)
				samplet = log(1.0-r)/Q[path[end],path[end]]
				totalt += samplet
				if totalt > 1.0
					break
				end
				push!(path, sample(rng, S[path[end],:]))
				push!(times, totalt)
				count += 1
				if count > 1000000
					println("ERRROR! ", Q, "\t", x0, "\t",xt,"\t", node.branchlength,"\t",generator.Q)
					exit()
				end
			end

			if path[end] == xt
				break
			end

			index = (node.nodeindex-1)*generator.alphabet*generator.alphabet + (x0-1)*generator.alphabet + path[end]
			if length(generator.paths[index]) < generator.maxcachesize
				push!(generator.paths[index], path)
				push!(generator.times[index], times)
			end
		end

		return path, times
	end

	#=
	export modifiedrejectionsampling
	function  modifiedrejectionsampling(rng::AbstractRNG, Q::Array{Float64,2}, x0::Int, xt::Int, extradata)
		maxiters1 = 1000
		maxiters2 = 2000
		path = Int[]
		times = Float64[]
		S = copy(Q)
		for i=1:size(S,1)
			S[i,i] = 0.0
			S[i,:] ./= sum(S[i,:])
		end

		count = 0
		success = true
		while success
			path = Int[x0]
			times = Float64[0.0]
			totalt = 0.0
			if x0 != xt
				r = rand(rng)
				totalt += log(1.0-r*(1.0-exp(Q[path[end],path[end]])))/Q[path[end],path[end]]
				push!(path, sample(rng, S[path[end],:]))
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
				push!(path, sample(rng, S[path[end],:]))
				push!(times, totalt)
				count += 1
			end
			if path[end] == xt
				break
			end

		end

		#println("H")
		if !success
			path,times = recursivesampling(rng, Q, x0, xt)
			#path, times = approximatesampling(rng, Q, 1.0, x0, xt)
			return path,times,true
		else
			return path, times, success
		end
	end=#

	export removevirtualjumps
	function removevirtualjumps(path::Array{Int,1}, time::Array{Float64,1})
		newpath = Int[]
		newtime = Float64[]
		prevstate = 0
		for (newstate,t) in zip(path,time)
			if prevstate != newstate
				push!(newpath, newstate)
				push!(newtime,t)
				prevstate = newstate
			end
		end
		return newpath,newtime
	end


	ts = Float64[]
	function gettindex(t::Float64)
		retindex = 0
		rettime = 0.0
		for (ind,v) in enumerate(ts)
			if abs(v-t) < 0.2*t
				retindex = ind
				rettime = v
				break
			end
		end
		if retindex == 0
			push!(ts,t)
			retindex = length(ts)
			rettime = t 
		end
		return retindex,rettime
	end

	pathdict = Dict{Tuple{Int,Int,Int,Int,Int,Int},Array{Any,1}}()
	export samplepath
	function samplepath(rng::AbstractRNG, h::Int, tin::Float64, R::Array{Float64,2}, aa1::Int, aa2::Int)
		tindex, t = gettindex(tin)
		key = (0,0,h,tindex,aa1,aa2)
		if haskey(pathdict, key) && length(pathdict[key]) > 0
			println("CACHE HIT")
			t,path,time = pop!(pathdict[key])
			return t,path,time,true
		else
			cache = []
			path,time,success = modifiedrejectionsampling3(rng,pathdict,key, R,t, aa1, aa2)
			return t,path,time,success
		end	
	end

	export resetcache
	function resetcache()
		global pathdict
		pathdict = Dict{Tuple{Int,Int,Int,Int},Array{Any,1}}()
		global ts
		ts = Float64[]
	end
end