module WrappedUnivariateOU
	using Distributions
	push!(LOAD_PATH,@__DIR__)
	using CommonUtils
	using AngleUtils
	using LinearAlgebra

	export WrappedUnivariateOUNode
	mutable struct WrappedUnivariateOUNode
		mu::Float64
		sigma::Float64
		alpha::Float64
		t::Float64
		branchlength::Float64
		obserror::Float64
		maxk::Int
		statdist::Normal
		transdist::Normal

		statdistvm::VonMises
		transdistvm::VonMises

		function WrappedUnivariateOUNode(mu::Float64=0.0, sigma::Float64=0.4, alpha::Float64=0.1, t::Float64=1.0, branchlength::Float64=1.0, obserror::Float64=0.01, maxk::Int=1)
			statdist = Normal(0.0, sqrt(sigma*sigma/(2.0*alpha)))
			transdist = Normal(0.0, sqrt((sigma*sigma/(2.0*alpha))*(1.0 - exp(-2.0*alpha*t))))
			statdistvm = VonMises(mu,1.0/sigma)
			transdistvm = VonMises(mu,1.0/(t*alpha))
			new(mu, sigma, alpha, t, branchlength, obserror, maxk, statdist, transdist, statdistvm, transdistvm)
		end
	end

	export get_bounds
	function get_bounds(wn::WrappedUnivariateOUNode)
		lower = Float64[0.0, 0.0, 0.0, 0.0]
		upper = Float64[2.0*pi, 1.0, 50.0, 0.1]
		return lower, upper
	end
	
	export get_parameters
	function get_parameters(wn::WrappedUnivariateOUNode)
		return Float64[wn.mu, wn.sigma, wn.alpha, wn.obserror]
	end

	export get_transformed_parameters
	function get_transformed_parameters(wn::WrappedUnivariateOUNode)
		transsigma = wn.sigma*wn.sigma/(2.0*wn.alpha)
		return Float64[wn.mu, transsigma, wn.alpha, wn.obserror]
	end

	export set_transformed_parameters
	function set_transformed_parameters(wn::WrappedUnivariateOUNode, params::Array{Float64,1})
		transsigma = params[2]
		sigma = sqrt(transsigma*params[3]*2.0)
		set_parameters(wn, Float64[params[1], sigma, params[3], params[4]])
	end

	export set_parameters
	function set_parameters(wn::WrappedUnivariateOUNode, params::Array{Float64,1})
		wn.mu = mod2pi(params[1])
		wn.sigma = params[2]
		wn.alpha = params[3]
		wn.obserror = params[4]
		wn.statdist = Normal(0.0, sqrt(wn.sigma*wn.sigma/(2.0*wn.alpha)))
		wn.transdist = Normal(0.0, sqrt((wn.sigma*wn.sigma/(2.0*wn.alpha))*(1.0 - exp(-2.0*wn.alpha*wn.t))))
		#wn.statdistvm = VonMises(wn.mu, min(2000.0, wn.sigma))
		#wn.transdistvm = VonMises(0.0, min(2000.0, 1.0/(wn.t*wn.alpha)))
	end

	export set_time
	function set_time(wn::WrappedUnivariateOUNode, t::Float64)
		if wn.t != t
			wn.t = t
			wn.transdist = Normal(0.0, sqrt((wn.sigma*wn.sigma/(2.0*wn.alpha))*(1.0 - exp(-2.0*wn.alpha*wn.t))))
			#wn.transdistvm = VonMises(0.0, min(2000.0, 1.0/(wn.t*wn.alpha)))
		end
	end

	export set_branchlength
	function set_branchlength(wn::WrappedUnivariateOUNode, branchlength::Float64, includeobserror::Bool=false)
		wn.branchlength = branchlength
		if includeobserror
			set_time(wn, branchlength+wn.obserror)
		else
			set_time(wn, branchlength)
		end
	end

	export logstat
	function logstat(wn::WrappedUnivariateOUNode, x::Float64)		
		theta = mod2pi(x)
		
		logconst = -Inf
		for m=-wn.maxk:wn.maxk
			logconst = logsumexp(logconst, logpdf(wn.statdist, theta - wn.mu + m*2.0*pi))
		end
		return logconst
		
		#return logpdf(wn.statdistvm, theta)
	end

	export logtpd
	function logtpd(wn::WrappedUnivariateOUNode, x0::Float64, xt::Float64)
		theta0 = mod2pi(x0)
		thetat = mod2pi(xt)
		
		#wmarray, logconst = w(wn.mu, sqrt(wn.sigma*wn.sigma/(2.0*wn.alpha)), theta0, wn.maxk)
		logll = -Inf
		wmlogconstant = -Inf
		for m=-wn.maxk:wn.maxk
			mut = wn.mu + (theta0 - wn.mu + m*2.0*pi)*exp(-wn.alpha*wn.t)
			wm = logpdf(wn.statdist, theta0 - wn.mu + m*2.0*pi)
			wmlogconstant = logsumexp(wmlogconstant, wm)
			wrappedll = -Inf
			for k=-wn.maxk:wn.maxk
				wrappedll = logsumexp(wrappedll, logpdf(wn.transdist, thetat - mut + k*2.0*pi))
			end
			logll = logsumexp(logll, wrappedll+wm)
			#logll = logsumexp(logll, logpdf(wn.transdist, thetat-mut)+wm)
		end
		
		#=
		angledist = abs(x0-xt)
		if angledist >= 0.355421202
			logll -= 4.60517018598
		elseif angledist >= 0.177007686
			logll -= 2.30258509299
		end=#

		return logll - wmlogconstant		
		#return logpdf(wn.transdistvm, thetat-theta0)
	end

	#=
	export gettransitionmatrix
	function gettransitionmatrix(wn::WrappedUnivariateOUNode, numcats::Int, t::Float64=wn.t)
		set_time(wn,t)
		mat = zeros(Float64, numcats, numcats)
		logconsts = zeros(Float64, numcats)
		for x=1:numcats
			anglex = indextoangle(x, numcats)
			for y=1:numcats
				angley = indextoangle(y, numcats)
				mat[x,y] = logtpd(wn, anglex, angley)
			end
			logconsts[x] = maximum(mat[x,:])
			mat[x,:] = exp.(mat[x,:] .- logconsts[x])
		end
		return mat, logconsts
	end=#

	export getlogstatrow
	function getlogstatrow(wn::WrappedUnivariateOUNode, numcats::Int)
		v = zeros(Float64, numcats)
		for z=1:numcats
			v[z] = logstat(wn, indextoangle(z, numcats))
		end
		return v
	end

	export gettransitionmatrixrow
	function gettransitionmatrixrow(wn::WrappedUnivariateOUNode, numcats::Int, branchlength::Float64, includeobserror::Bool, row::Int)
		set_branchlength(wn, branchlength, includeobserror)
		if wn.t == 0.0
			return Matrix{Float64}(I,numcats,numcats)[row,:], 0.0
		else
			rowvec = zeros(Float64, numcats)
			anglex = indextoangle(row, numcats)
			for y=1:numcats
				angley = indextoangle(y, numcats)
				rowvec[y] = logtpd(wn, anglex, angley)
			end
			logconst = maximum(rowvec)
			rowvec = exp.(rowvec.-logconst)
			for v in rowvec
				if isnan(v)
					println("NAN")
					println(wn.t,"\t",branchlength,"\t",includeobserror,"\t",row)
					println(rowvec)
				end
			end
			return rowvec, logconst
		end
	end

	export gettransitionmatrix
	function gettransitionmatrix(wn::WrappedUnivariateOUNode, numcats::Int, branchlength::Float64, includeobserror::Bool)
		set_branchlength(wn, branchlength, includeobserror)
		if wn.t == 0.0
			return Matrix{Float64}(I,numcats,numcats), 0.0
		else
			mat = zeros(Float64, numcats, numcats)
			for x=1:numcats
				anglex = indextoangle(x, numcats)
				for y=1:numcats
					angley = indextoangle(y, numcats)
					mat[x,y] = logtpd(wn, anglex, angley)
				end
			end
			logconst = maximum(mat)
			return exp.(mat.-logconst), logconst
		end
	end

	export logjoint
	function logjoint(wn::WrappedUnivariateOUNode, x0::Float64, xt::Float64)
		return logstat(wn, x0) + logtpd(wn, x0, xt)
	end

	export estimateparameters
	function estimateparameters(data::Array{Float64,1})
		rx = 0.0
		ry = 0.0
		N = 0.0
		for theta in data
			rx += cos(mod2pi(theta))
			ry += sin(mod2pi(theta))
			N += 1.0
		end

		c = rx / N
		s = ry / N
		rho = sqrt(c*c + s*s)

		mu = 0.0
		kappa = 0.0
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
		return mu, kappa
	end

	#=
	function w(mu::Float64, sigma::Float64, x::Float64, maxk::Int, dist::Normal=Normal(0.0,sigma))		
		logliks = zeros(Float64,2*maxk+1)
		logconst = -Inf
		for m=-maxk:maxk
			index = m+maxk+1
			logliks[index] = logpdf(dist, x - mu + m*2.0*pi)
			logconst = logsumexp(logconst, logliks[index])
		end
		return logliks .- logconst, logconst
	end
	function logoustat(wn::WrappedUnivariateOUNode, x0::Float64)
		mean = wn.mu
		sigma = wn.sigma*wn.sigma/(2.0*wn.alpha)
	end

	function logoutpd(wn::WrappedUnivariateOUNode, x0::Float64, xt::Float64)
		mean = wn.mu + (x0 - wn.mu)*exp(-wn.alpha*wn.t)
		sigma = (wn.sigma*wn.sigma/(2.0*wn.alpha))*(1.0 - exp(-2.0*wn.alpha*wn.t))
	end=#
end
#println(WrappedUnivariateOU.estimateparameters(Float64[0.1,0.3,0.111]))
#=
wn = WrappedUnivariateOU.WrappedUnivariateOUNode()
WrappedUnivariateOU.set_parameters(wn, Float64[0.0, 1.1, 0.1])
println(WrappedUnivariateOU.get_parameters(wn))
transparams = WrappedUnivariateOU.get_transformed_parameters(wn)
println(transparams)
WrappedUnivariateOU.set_transformed_parameters(wn, transparams)
println(WrappedUnivariateOU.get_parameters(wn))=#
#=
wn = WrappedUnivariateOU.WrappedUnivariateOUNode()
WrappedUnivariateOU.set_transformed_parameters(wn,Float64[0.0, 1.0, 1.0])
WrappedUnivariateOU.set_time(wn, 1.0)
println(WrappedUnivariateOU.gettransitionmatrix(wn, 20))=#
#=
for x=-2*pi:0.1:2*pi
	#println(x, "\t", WrappedUnivariateOU.logstat(wn, x))
	println(x, "\t", WrappedUnivariateOU.logtpd(wn, 1.0*pi, x))
end=#
#=
println(wn)
println(WrappedUnivariateOU.logstat(wn, 0.1))
println(WrappedUnivariateOU.logtpd(wn, 0.1,0.2))=#