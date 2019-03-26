module EMNodes
	using Distributions
	using Random
	using LinearAlgebra
	using NLopt

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

		function SiteObservation()
			new(0, 0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0)
		end

		function SiteObservation(h::Int, aa::Int, phi::Float64, omega::Float64, psi::Float64)
			new(h,aa,phi,omega,psi, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0, -1000.0)
		end

		function SiteObservation(h::Int, aa::Int, phi::Float64, omega::Float64, psi::Float64, bond_length1::Float64, bond_length2::Float64, bond_length3::Float64, bond_angle1::Float64, bond_angle2::Float64, bond_angle3::Float64)
			new(h,aa,phi,omega,psi, bond_length1, bond_length2, bond_length3, bond_angle1, bond_angle2, bond_angle3)
		end
	end

	export Protein
	mutable struct Protein
		name::String
		sites::Array{SiteObservation,1}

		function Protein()
			new("", SiteObservation[])
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
		categorical_node.counts =  ones(Float64, categorical_node.numcats)*0.1
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

		function VonMisesNode()
			mu = 0.0
			kappa = 1e-5
			new(0.0, 0.0, 0.0, mu, kappa, VonMises(mu, kappa), Float64[], Exponential(2.0))
		end
	end

	export estimatevonmises
	function estimatevonmises(vonmises_node::VonMisesNode)
		if length(vonmises_node.data) > 0
			vonmises_node.rx = 0.0
			vonmises_node.ry = 0.0
			vonmises_node.N = 0
			for theta in vonmises_node.data
				if theta > -100.0
					vonmises_node.rx += cos(theta)
					vonmises_node.ry += sin(theta)
					vonmises_node.N += 1
				end
			end		
		end

		if vonmises_node.N >= 2
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
			vonmises_node.dist = VonMises(mu, min(700.0, kappa))
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
		println("start value: ", vonmises_node.kappa, "\t", vonmises_node.N)
		if vonmises_node.N == 0
			vonmises_node.kappa = 1e-5
		else
			opt = Opt(:LN_COBYLA,1)
			data = filter(x -> x > -100.0, vonmises_node.data)
			localObjectiveFunction = ((param, grad) -> loglik_and_prior(vonmises_node, data, param[1]))
			lower = ones(Float64, 1)*1e-5
			upper = ones(Float64, 1)*700.0
			lower_bounds!(opt, lower)
			upper_bounds!(opt, upper)
			xtol_rel!(opt,1e-5)
			maxeval!(opt, 1000)
			max_objective!(opt, localObjectiveFunction)
			(minf,minx,ret) = optimize(opt, Float64[min(vonmises_node.kappa, 700.0)])
			vonmises_node.kappa = minx[1]
		end
		vonmises_node.dist = VonMises(vonmises_node.mu, vonmises_node.kappa)
		println("end value: ", vonmises_node.kappa)
	end

	export vonmisesrand
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

	export HiddenNode
	mutable struct HiddenNode
		aa_node::CategoricalNode
		phi_node::VonMisesNode
		psi_node::VonMisesNode
		omega_node::VonMisesNode
		bond_angle1_node::VonMisesNode
		bond_angle2_node::VonMisesNode
		bond_angle3_node::VonMisesNode
		bond_lengths_node::MultivariateNode
		phi_nodes::Array{VonMisesNode,1}
		psi_nodes::Array{VonMisesNode,1}
		omega_nodes::Array{VonMisesNode,1}
	    #aadist::Array{Float64,1}

	    function HiddenNode(aafreqs::Array{Float64}=ones(Float64,20)*0.05)
	    	phi_nodes = VonMisesNode[VonMisesNode() for aa=1:20]
	    	psi_nodes = VonMisesNode[VonMisesNode() for aa=1:20]
	    	omega_nodes = VonMisesNode[VonMisesNode() for aa=1:20]
	        new(CategoricalNode(aafreqs),VonMisesNode(),VonMisesNode(),VonMisesNode(),VonMisesNode(),VonMisesNode(),VonMisesNode(), MultivariateNode(), phi_nodes, psi_nodes, omega_nodes)
	    end
	end
end