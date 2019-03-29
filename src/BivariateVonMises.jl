module BivariateVonMises
	#= 
	Adapted from Mocapy++
 	vonmises2d.h --- Bivariate von Mises distribution, cosine variant
	Copyright (C) 2008, Wouter Boomsma, The Bioinformatics Centre, University of Copenhagen.
	=#
	using Distributions
	using SpecialFunctions
	using QuadGK
	using NLopt

	export BivariateVonMisesDist
	mutable struct BivariateVonMisesDist
		k::Array{Float64,1}
		mu::Array{Float64,1}
		logNormConst::Float64
		ess2COS1::Float64
		ess2SIN1::Float64
		ess2COS2::Float64
		ess2SIN2::Float64
		essCOS1MIN2::Float64
		essSIN1MIN2::Float64
		essCOS1PLUS2::Float64
		essSIN1PLUS2::Float64
		essC1::Float64
		essS1::Float64
		essC2::Float64
		essS2::Float64
		count::Float64
		data::Array{Array{Float64,1},1}
		kappa_prior::ContinuousUnivariateDistribution
		
		function BivariateVonMisesDist(k::Array{Float64,1}, mu::Array{Float64,1},kappa_prior_exp_rate::Float64=0.5)
			new(k,mu, computeLogNormConst(k),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, Array{Float64,1}[], Exponential(kappa_prior_exp_rate))
		end
	end
"""
	// Exact integrand for normalization
double logNormConstExactIntegrand(double *phi, double *extraArguments) {
     double k1 = extraArguments[0];
     double k2 = extraArguments[1];
     double k3 = extraArguments[2];
     double k13 = sqrt((k1*k1) + (k3*k3) - 2*k1*k3*cos(*phi));
     double res = (2*M_PI*exp(log(i0(k13)) + k2*cos(*phi)));
     return res;
}"""
	function logNormConstExactIntegrand(phi::Float64, k::Array{Float64,1})
		k13 = sqrt((k[1]*k[1]) + (k[3]*k[3]) - 2.0*k[1]*k[3]*cos(phi))
		#println(k13,"\t",besseli(0.0, k13))
		return 2.0*pi*exp(log(max(0.0, besseli(0.0, k13))) + k[2]*cos(phi))
	end


"""
// Integrand using approximation for Bessel function (works for large values of k13)
double logNormConstApproximateIntegrand(double *phi, double *extraArguments) {
     double k1 = extraArguments[0];
     double k2 = extraArguments[1];
     double k3 = extraArguments[2];
     double scaleTerm = extraArguments[3];
     double k13 = sqrt((k1*k1) + (k3*k3) - 2*k1*k3*cos(*phi));
     double e = exp(k13 + k2*cos(*phi) - scaleTerm);
     return sqrt((2*M_PI)/k13)*e;
}"""


	function logNormConstApproximateIntegrand(phi::Float64, k::Array{Float64,1}, scaleTerm::Float64)
	  k13 = sqrt((k[1]*k[1]) + (k[3]*k[3]) - 2.0*k[1]*k[3]*cos(phi))
      e = exp(k13 + k[2]*cos(phi) - scaleTerm)
      return sqrt((2.0*pi)/k13)*e
	end


#=
double k1 = k[0];
     double k2 = k[1];
     double k3 = k[2];
     double exactLimit = 400;
     if ((sqrt(k1*k1 + k2*k2 - 2*k1*k3) + k2) < exactLimit) {
          double extraArguments[3] = {k1, k2, k3};
          double resExact = integrateQuad(&logNormConstExactIntegrand, extraArguments, lower, upper);
          return (log(resExact));
     } else {
          double scaleTerm = sqrt(k1*k1 + k3*k3 - 2*k1*k3) + k2;
          double extraArguments[4] = {k1, k2, k3, scaleTerm};
          double resApproximate = integrateQuad(&logNormConstApproximateIntegrand, extraArguments, lower, upper);
          return (log(resApproximate) + scaleTerm);
     }
     return 0;
}=#

	function computeLogNormConst(k::Array{Float64,1}, lower::Float64=Float64(-pi), upper::Float64=Float64(pi))
	     exactLimit = 400.0
	     beta = 
	     if ((sqrt(k[1]*k[1] + k[2]*k[2] - 2.0*k[1]*k[3]) + k[2]) < exactLimit)
	     	 println("A")
	          #resExact = integrateQuad(logNormConstExactIntegrand, k, lower, upper);
	          resExact,err = quadgk(x -> logNormConstExactIntegrand(x,k), lower, upper)
	          return log(resExact)
	     else
	     	println("B")
	          scaleTerm = sqrt(k[1]*k[1] + k[3]*k[3] - 2.0*k[1]*k[3]) + k[2]
	          #resApproximate = integrateQuad(logNormConstExactIntegrand, k, lower, upper)
	          resApproximate,err = quadgk(x -> logNormConstApproximateIntegrand(x,k,scaleTerm), lower, upper)
	          
	          resExact,err = quadgk(x -> logNormConstExactIntegrand(x,k), lower, upper)
	          println("comp: ",log(resExact),"\t", log(resApproximate) + scaleTerm,"\t",scaleTerm)

	          return log(resApproximate) + scaleTerm
	     end
	end

	function logpdf(dist::BivariateVonMisesDist, x::Array{Float64,1})
		return dist.k[1]*cos(x[1]-dist.mu[1]) + dist.k[2]*cos(x[2]-dist.mu[2]) - dist.k[3]*cos(x[1]-dist.mu[1]-x[2]+dist.mu[2]) - dist.logNormConst
	end

	function pdf(dist::BivariateVonMisesDist, x::Array{Float64,1})
		return exp(logpdf(dist,x))
	end


	function addpoint(dist::BivariateVonMisesDist, x::Array{Float64,1})
		push!(dist.data, x)
		dist.ess2COS1 += cos(2*x[1])
		dist.ess2SIN1 += sin(2*x[1])
		dist.ess2COS2 += cos(2*x[2])
		dist.ess2SIN2 += sin(2*x[2])
		dist.essCOS1MIN2 += cos(x[1] - x[2])
		dist.essSIN1MIN2 += sin(x[1] - x[2])
		dist.essCOS1PLUS2 += cos(x[1] + x[2])
		dist.essSIN1PLUS2 += sin(x[1] + x[2])
		dist.essC1 += cos(x[1])
		dist.essS1 += sin(x[1])
		dist.essC2 += cos(x[2])
		dist.essS2 += sin(x[2])
		dist.count += 1.0
	end

	function reset(dist::BivariateVonMisesDist)
		dist.data = Array{Float64,1}[]
		dist.ess2COS1 = 0.0
		dist.ess2SIN1 = 0.0
		dist.ess2COS2 = 0.0
		dist.ess2SIN2 = 0.0
		dist.essCOS1MIN2 = 0.0
		dist.essSIN1MIN2 = 0.0
		dist.essCOS1PLUS2 = 0.0
		dist.essSIN1PLUS2 = 0.0
		dist.essC1 = 0.0
		dist.essS1 = 0.0
		dist.essC2 = 0.0
		dist.essS2 = 0.0
		dist.count = 0.0
	end

	function loglik_and_prior(bv::BivariateVonMisesDist, data::Array{Array{Float64,1},1}, k::Array{Float64,1})
		try
			temp = BivariateVonMisesDist(Float64[k[1], k[2], k[3]], Float64[k[4], k[5]])
			loglik = sum(Float64[logpdf(temp,x) for x in data]) 
			loglik += Distributions.logpdf(bv.kappa_prior, k[1]) 
			loglik += Distributions.logpdf(bv.kappa_prior, k[2])
			loglik += Distributions.logpdf(bv.kappa_prior, abs(k[3]))
			#loglik = sum(Float64[logpdf(temp,x) for x in data])
			println(loglik,"\t", k)
			return loglik
		catch e
			println("error", e)
			return -Inf
		end
	end

	function optimizeprior(bv::BivariateVonMisesDist)
		startstr = string("start value: ", bv.k, "\t", bv.mu)
		if bv.count == 0
			bv.k = Float64[1e-5,1e-5]
		else
			opt = Opt(:LN_COBYLA,5)
			data = filter(x -> x[1] > -100.0 && x[2] > -100.0, bv.data)
			localObjectiveFunction = ((params, grad) -> loglik_and_prior(bv, data, params))
			lower = ones(Float64, 5)*1e-5
			upper = ones(Float64, 5)*700.0
			lower[3] = -700.0
			lower[4] = Float64(-pi)
			lower[5] = Float64(-pi)
			upper[4] = Float64(pi)
			upper[5] = Float64(pi)
			lower_bounds!(opt, lower)
			upper_bounds!(opt, upper)
			xtol_rel!(opt,1e-5)
			maxeval!(opt, 1000)
			max_objective!(opt, localObjectiveFunction)
			(minf,minx,ret) = optimize(opt, Float64[min(bv.k[1], 700.0), min(bv.k[2], 700.0), max(-700.0, min(bv.k[3], 700.0)), bv.mu[1], bv.mu[2]])
			bv.k[1] = minx[1]
			bv.k[2] = minx[2]
			bv.k[3] = minx[3]
			bv.mu[1] = minx[4]
			bv.mu[2] = minx[5]
			println(ret)
		end		
		println(startstr)
		println("end value: ", bv.k, "\t", bv.mu)
	end

	#=
	function loglik_and_prior(bv::BivariateVonMisesDist, data::Array{Array{Float64,1},1}, k::Array{Float64,1})
		try
			temp = BivariateVonMisesDist(Float64[k[1],k[2], bv.k[3]], bv.mu)
			#loglik = sum(Float64[logpdf(temp,x) for x in data]) + Distributions.logpdf(bv.kappa_prior, k[1]) + Distributions.logpdf(bv.kappa_prior, k[2])
			loglik = sum(Float64[logpdf(temp,x) for x in data])
			println(loglik,"\t", k)
			return loglik
		catch e
			println("error", e)
			return -Inf
		end
	end

	function optimizeprior(bv::BivariateVonMisesDist)
		println("start value: ", bv.k)
		if bv.count == 0
			bv.k = Float64[1e-5,1e-5]
		else
			opt = Opt(:LN_COBYLA,2)
			data = filter(x -> x[1] > -100.0 && x[2] > -100.0, bv.data)
			localObjectiveFunction = ((params, grad) -> loglik_and_prior(bv, data, params))
			lower = ones(Float64, 2)*1e-5
			upper = ones(Float64, 2)*700.0
			lower_bounds!(opt, lower)
			upper_bounds!(opt, upper)
			xtol_rel!(opt,1e-5)
			maxeval!(opt, 1000)
			max_objective!(opt, localObjectiveFunction)
			(minf,minx,ret) = optimize(opt, Float64[min(bv.k[1], 700.0), min(bv.k[2], 700.0)])
			bv.k[1] = minx[1]
			bv.k[2] = minx[2]
			println(ret)
		end
		println("end value: ", bv.k)
	end=#

	export estimate
	function estimate(dist::BivariateVonMisesDist)
		mu1 = atan(dist.essS1, dist.essC1)
		mu2 = atan(dist.essS2, dist.essC2)

		S1 = 0.5 - ((1.0/(2.0*dist.count))*cos(2*mu1)*dist.ess2COS1 + (1.0/(2.0*dist.count))*sin(2*mu1)*dist.ess2SIN1)     	
     	S2 = 0.5 - ((1.0/(2.0*dist.count))*cos(2*mu2)*dist.ess2COS2 + (1.0/(2.0*dist.count))*sin(2*mu2)*dist.ess2SIN2)
     	S12 = ((1.0/(2.0*dist.count)) * cos(mu1-mu2) * dist.essCOS1MIN2 + (1.0/(2.0*dist.count)) * sin(mu1-mu2) * dist.essSIN1MIN2 - (1.0/(2.0*dist.count)) * cos(mu1+mu2) * dist.essCOS1PLUS2 - (1.0/(2.0*dist.count)) * sin(mu1+mu2) * dist.essSIN1PLUS2)

     	k1 = (S2-S12)/(S1*S2-(S12*S12))
    	k2 = (S1-S12)/(S1*S2-(S12*S12))
    	k3 = (-S12)/(S1*S2-(S12*S12))

		if k1 < 0.0
			mu1 -= pi
			k1 *= -1.0
		end

		if k2 < 0.0
			mu2 -= pi
			k2 *= -1.0
		end

		dist.k = Float64[k1,k2,k3]
		dist.mu = Float64[mu1,mu2]

		optimizeprior(dist)
		reset(dist)
 	end
end

bv = BivariateVonMises.BivariateVonMisesDist(Float64[2.0,4.0,0.0], Float64[0.0,0.0])
#println(BivariateVonMises.logpdf(bv,[0.0,0.0]))
for z=1:1000
	x = 0.5 + randn()*0.1
	y = 0.8 + randn()*0.1 - x*0.5
	# + x*0.02 +
	BivariateVonMises.addpoint(bv, Float64[x,y])
end
BivariateVonMises.estimate(bv)
println(bv.k,"\t",bv.mu)