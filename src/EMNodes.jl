module EMNodes
	using Distributions
	using Random
	using LinearAlgebra
	using NLopt

	push!(LOAD_PATH,@__DIR__)
	using BivariateVonMises
	using WrappedUnivariateOU

	export SiteObservation
	mutable struct SiteObservation
		h::Int
		aa::Int
		phi::Float64		
		omega::Float64
		psi::Float64
		bond_length1::Float64
		bond_length2::Float64
		bond_length3::Float64
		bond_angle1::Float64
		bond_angle2::Float64
		bond_angle3::Float64
		phiobserved::Bool
		psiobserved::Bool

		function SiteObservation()
			new(0, 0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, false, false)
		end

		function SiteObservation(h::Int, aa::Int, phi::Float64, omega::Float64, psi::Float64)
			new(h,aa,phi,omega,psi, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, phi > -100.0, psi > -100.0)
		end

		function SiteObservation(h::Int, aa::Int, phi::Float64, omega::Float64, psi::Float64, bond_length1::Float64, bond_length2::Float64, bond_length3::Float64, bond_angle1::Float64, bond_angle2::Float64, bond_angle3::Float64)
			new(h,aa,phi,omega,psi, bond_length1, bond_length2, bond_length3, bond_angle1, bond_angle2, bond_angle3, phi > -100.0, psi > -100.0)
		end
	end

	export Protein
	mutable struct Protein
		name::String
		sites::Array{SiteObservation,1}

		function Protein()
			new("", SiteObservation[])
		end

		function Protein(numcols::Int)
			new("", SiteObservation[SiteObservation() for col=1:numcols])
		end

		function Protein(name::String)
			new(name, SiteObservation[])
		end
	end

	function Base.length(protein::Protein)
		return length(protein.sites)
	end

	export CategoricalNode
	mutable struct CategoricalNode 
		numcats::Int 
		probs::Array{Float64,1}
		counts::Array{Float64,1}
		alphas::Array{Float64,1}

		function CategoricalNode(freqs::Array{Float64,1})
			numcats = length(freqs)
			new(numcats, freqs, ones(Float64, numcats)*0.1, ones(Float64,numcats))
		end
	end

	export estimate_categorical
	function estimate_categorical(categorical_node::CategoricalNode, beta::Float64=1.0)
		categorical_node.probs = categorical_node.counts / sum(categorical_node.counts)
		categorical_node.probs = categorical_node.probs.^categorical_node.alphas.^beta
		categorical_node.probs /= sum(categorical_node.probs)
		#categorical_optimizeprior(categorical_node)
		categorical_node.counts =  ones(Float64, categorical_node.numcats)*0.1
	end

	function categorical_loglik_and_prior(categorical_node::CategoricalNode, params::Array{Float64,1})
		s = sum(params)
		if s > 1.0
			ll = -s*1e40
			##println("FF", ll)
			return -Inf
		elseif s < 1e-4
			ll = -(1.0-s)*1e40
			##println("EE", ll)
			return -Inf
		else
			alpha = 1e-10
			d = Dirichlet(alpha*ones(Float64,20))
			probs = Float64[]
			append!(probs, params)
			push!(probs, 1.0-s)
			logprobs = log.(probs)
			ll = sum(logprobs.*categorical_node.counts) + logpdf(d,probs)
			println("W\t",ll,"\t",probs)
			return ll
		end
	end

	function categorical_optimizeprior(categorical_node::CategoricalNode)
		opt = Opt(:LN_COBYLA,19)
		localObjectiveFunction = ((param, grad) -> categorical_loglik_and_prior(categorical_node, param))
		lower = ones(Float64, 19)*1e-4
		upper = ones(Float64, 19)
		upper_bounds!(opt, upper)
		xtol_rel!(opt,1e-5)
		maxeval!(opt, 2000)
		max_objective!(opt, localObjectiveFunction)
		initialparams = max.(1e-4, categorical_node.probs[1:19])
		initialprobs = copy(categorical_node.probs)
		initialloglik = categorical_loglik_and_prior(categorical_node, initialparams)
		(minf,minx,ret) = optimize(opt, initialparams)
		probs = Float64[]
		append!(probs,minx)
		push!(probs, 1.0 - sum(minx))
		categorical_node.probs = probs
		println(minf,"\t",initialloglik)
		println("INIT ", initialprobs)
		println("FINAL ", probs)
	end

	export MultivariateNode
	mutable struct MultivariateNode
		mvn::MvNormal
		data::Array{Array{Float64,1},1}

		function MultivariateNode()
			new(MvNormal([1.0,1.0,1.0]), Array{Float64,1}[])
		end
	end

	export add_point
	function add_point(mv_node::MultivariateNode, point::Array{Float64,1})
		cont = true
		for p in point
			if p <= -100.0
				cont = false
				break
			end
		end
		if cont
			push!(mv_node.data, point)
		end
	end

	export estimate_multivariate_node
	function estimate_multivariate_node(mv_node::MultivariateNode)	
		println("N=",length(mv_node.data))
		try
			matrix = zeros(Float64, length(mv_node.data[1]), length(mv_node.data))
			for col=1:size(matrix,2)
				for (row,el) in enumerate(mv_node.data[col])
					matrix[row,col] = el
				end
			end

			mvn_suffstats = suffstats(MvNormal, matrix)
			mv_node.mvn =  fit_mle(MvNormal, mvn_suffstats)
		catch Exception 
			mv_node.mvn = MvNormal(Float64[1.32, 1.47,1.53], Matrix{Float64}(I,3,3)*0.1)
		end
	end

	export VonMisesNode
	mutable struct VonMisesNode
		rx::Float64
		ry::Float64
		N::Float64
		mu::Float64
		kappa::Float64
		dist::VonMises
		data::Array{Float64,1}
		kappa_prior::ContinuousUnivariateDistribution
		maxkappa::Float64

		function VonMisesNode(kappa_prior_exp_rate::Float64=100.0)
			mu = 0.0
			kappa = 1e-5
			new(0.0, 0.0, 0.0, mu, kappa, VonMises(mu, kappa), Float64[], Exponential(kappa_prior_exp_rate), 5000.0)
		end
	end

	export estimatevonmises
	function estimatevonmises(vonmises_node::VonMisesNode)
		if length(vonmises_node.data) > 0
			vonmises_node.rx = 0.0
			vonmises_node.ry = 0.0
			vonmises_node.N = 0.0
			for theta in vonmises_node.data
				if theta > -100.0
					vonmises_node.rx += cos(theta)
					vonmises_node.ry += sin(theta)
					vonmises_node.N += 1.0
				end
			end		
		end

		if vonmises_node.N >= 2.0
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
			vonmises_node.dist = VonMises(mu, min(vonmises_node.maxkappa, kappa))
		end
		optimizeprior(vonmises_node)
		vonmises_node.data = Float64[]

	end

	function loglik_and_prior(vonmises_node::VonMisesNode, data::Array{Float64,1}, kappa::Float64)
		loglik = loglikelihood(VonMises(vonmises_node.mu, kappa), data) + logpdf(vonmises_node.kappa_prior, kappa)
		#println(loglik,"\t",kappa)
		return loglik
	end

	function optimizeprior(vonmises_node::VonMisesNode)
		#println("start value: ", vonmises_node.kappa, "\t", vonmises_node.N)
		if vonmises_node.N == 0
			vonmises_node.kappa = 1e-5
		else
			opt = Opt(:LN_COBYLA,1)
			data = filter(x -> x > -100.0, vonmises_node.data)
			localObjectiveFunction = ((param, grad) -> loglik_and_prior(vonmises_node, data, param[1]))
			lower = ones(Float64, 1)*0.1
			upper = ones(Float64, 1)*vonmises_node.maxkappa
			lower_bounds!(opt, lower)
			upper_bounds!(opt, upper)
			xtol_rel!(opt,1e-5)
			maxeval!(opt, 1000)
			max_objective!(opt, localObjectiveFunction)
			(minf,minx,ret) = optimize(opt, Float64[max(0.1, min(vonmises_node.kappa, vonmises_node.maxkappa))])
			vonmises_node.kappa = minx[1]
		end
		vonmises_node.dist = VonMises(vonmises_node.mu, vonmises_node.kappa)
		#println("end value: ", vonmises_node.kappa)
	end

	export HiddenNode
	mutable struct HiddenNode
		aa_node::CategoricalNode
		phi_node::VonMisesNode
		psi_node::VonMisesNode
		phipsi_node::BivariateVonMisesNode
		omega_node::VonMisesNode
		bond_angle1_node::VonMisesNode
		bond_angle2_node::VonMisesNode
		bond_angle3_node::VonMisesNode
		bond_lengths_node::MultivariateNode
		phi_nodes::Array{VonMisesNode,1}
		psi_nodes::Array{VonMisesNode,1}
		phipsi_nodes::Array{BivariateVonMisesNode,1}
		omega_nodes::Array{VonMisesNode,1}

		#=
		obsinvkappa::Float64
		diffusionrate::Float64
		diffusionnode::WrappedNormalOUNode
		diffusionnode_uncorrelated::WrappedNormalOUNode=#
		phi_diffusion_node::WrappedUnivariateOUNode
		psi_diffusion_node::WrappedUnivariateOUNode

	    #aadist::Array{Float64,1}

	    function HiddenNode(aafreqs::Array{Float64}=ones(Float64,20)*0.05)
	    	phi_nodes = VonMisesNode[VonMisesNode() for aa=1:20]
	    	psi_nodes = VonMisesNode[VonMisesNode() for aa=1:20]
	    	phipsi_nodes = BivariateVonMisesNode[BivariateVonMisesNode() for aa=1:20]
	    	omega_nodes = VonMisesNode[VonMisesNode() for aa=1:20]
	    	
	    	#=
	    	wn = WrappedNormalOUNode(1)
	    	WrappedNormalOU.set_parameters(wn, zeros(Float64,2), Float64[0.5,0.5,0.0], Float64[1.0,1.0])
	    	set_parameters(wn, 1.0)
	    	
	    	wnuncorrelated = WrappedNormalOUNode(1)
	    	WrappedNormalOU.set_parameters(wnuncorrelated, WrappedNormalOU.get_parameters(wn))
	    	set_parameters(wnuncorrelated, 1.0)
	    	
	        new(CategoricalNode(aafreqs),VonMisesNode(),VonMisesNode(),BivariateVonMisesNode(),VonMisesNode(),VonMisesNode(),VonMisesNode(),VonMisesNode(), MultivariateNode(), phi_nodes, psi_nodes, phipsi_nodes, omega_nodes,1e-4, 0.1, wn, wnuncorrelated)=#
	        new(CategoricalNode(aafreqs),VonMisesNode(),VonMisesNode(),BivariateVonMisesNode(),VonMisesNode(),VonMisesNode(),VonMisesNode(),VonMisesNode(), MultivariateNode(), phi_nodes, psi_nodes, phipsi_nodes, omega_nodes, WrappedUnivariateOUNode(), WrappedUnivariateOUNode())
	    end
	end
end