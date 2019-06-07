module WrappedUnivariateOU
	using Distributions
	push!(LOAD_PATH,@__DIR__)
	using CommonUtils
	using AngleUtils

	export WrappedUnivariateOUNode
	mutable struct WrappedUnivariateOUNode
		mu::Float64
		sigma::Float64
		alpha::Float64
		t::Float64
		maxk::Int
		statdist::Normal
		transdist::Normal

		function WrappedUnivariateOUNode(mu::Float64=0.0, sigma::Float64=2.0, alpha::Float64=1.0, t::Float64=1.0, maxk::Int=1)
			statdist = Normal(0.0, sqrt(sigma*sigma/(2.0*alpha)))
			transdist = Normal(0.0, sqrt((sigma*sigma/(2.0*alpha))*(1.0 - exp(-2.0*alpha*t))))
			new(mu, sigma, alpha, t, maxk, statdist, transdist)
		end
	end

	export get_bounds
	function get_bounds(wn::WrappedUnivariateOUNode)
		lower = Float64[-4.0*pi, 0.0, 0.0]
		upper = Float64[4.0*pi, 1000.0, 1000.0]
		return lower, upper
	end
	
	export get_parameters
	function get_parameters(wn::WrappedUnivariateOUNode)
		return Float64[wn.mu, wn.sigma, wn.alpha]
	end

	export get_transformed_parameters
	function get_transformed_parameters(wn::WrappedUnivariateOUNode)
		transsigma = wn.sigma*wn.sigma/(2.0*wn.alpha)
		return Float64[wn.mu, transsigma, wn.alpha]
	end

	export set_transformed_parameters
	function set_transformed_parameters(wn::WrappedUnivariateOUNode, params::Array{Float64,1})
		transsigma = params[2]
		transalpha = params[3]
		sigma = sqrt(transsigma*transalpha*2.0)
		set_parameters(wn, Float64[params[1], sigma, params[3]])
	end

	export set_parameters
	function set_parameters(wn::WrappedUnivariateOUNode, params::Array{Float64,1})
		mu = mod2pi(params[1])
		sigma = params[2]
		alpha = params[3]
		wn.mu = mu
		wn.sigma = sigma
		wn.alpha = alpha
		wn.statdist = Normal(0.0, sqrt(sigma*sigma/(2.0*alpha)))
		wn.transdist = Normal(0.0, sqrt((sigma*sigma/(2.0*alpha))*(1.0 - exp(-2.0*alpha*wn.t))))
	end

	export set_time
	function set_time(wn::WrappedUnivariateOUNode, t::Float64)
		if wn.t != t
			wn.t = t
			wn.transdist = Normal(0.0, sqrt((wn.sigma*wn.sigma/(2.0*wn.alpha))*(1.0 - exp(-2.0*wn.alpha*wn.t))))
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
	end

	export logtpd
	function logtpd(wn::WrappedUnivariateOUNode, x0::Float64, xt::Float64)
		theta0 = mod2pi(x0)
		thetat = mod2pi(xt)
		#wmarray, logconst = w(wn.mu, sqrt(wn.sigma*wn.sigma/(2.0*wn.alpha)), theta0, wn.maxk)
		logll = -Inf
		wmlogconstant = -Inf
		for m=-wn.maxk:wn.maxk			
			index = m+wn.maxk+1
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
		return logll - wmlogconstant
	end

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
	end

	#=
	export gettransitionmatrix
	function gettransitionmatrix(wn::WrappedUnivariateOUNode, numcats::Int, t::Float64=wn.t)
		set_time(wn,t)
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
	end=#

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

wn = WrappedUnivariateOU.WrappedUnivariateOUNode()
WrappedUnivariateOU.set_parameters(wn,Float64[0.0, 1.0, 1.0])
WrappedUnivariateOU.set_time(wn, 1.0)
println(WrappedUnivariateOU.gettransitionmatrix(wn,20))

for x=-2*pi:0.1:2*pi
	#println(x, "\t", WrappedUnivariateOU.logstat(wn, x))
	println(x, "\t", WrappedUnivariateOU.logtpd(wn, 1.0*pi, x))
end
#=
println(wn)
println(WrappedUnivariateOU.logstat(wn, 0.1))
println(WrappedUnivariateOU.logtpd(wn, 0.1,0.2))=#