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
	using Random

	function pimod(angle::Float64)
	  theta = mod2pi(angle)
	  if theta > pi
	    return theta -2.0*pi
	  else
	    return theta
	  end
	end

	export BivariateVonMisesNode
	mutable struct BivariateVonMisesNode
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
		
		function BivariateVonMisesNode(k::Array{Float64,1}=Float64[1e-5,1e-5,0.0], mu::Array{Float64,1}=Float64[0.0,0.0],kappa_prior_exp_rate::Float64=100.0)
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
	     if((sqrt(k[1]*k[1] + k[2]*k[2] - 2.0*k[1]*k[3]) + k[2]) < exactLimit)
	          resExact,err = quadgk(x -> logNormConstExactIntegrand(x,k), lower, upper)
	          return log(resExact)
	     else
	          scaleTerm = sqrt(k[1]*k[1] + k[3]*k[3] - 2.0*k[1]*k[3]) + k[2]
	          resApproximate,err = quadgk(x -> logNormConstApproximateIntegrand(x,k,scaleTerm), lower, upper)	          
	          #resExact,err = quadgk(x -> logNormConstExactIntegrand(x,k), lower, upper)
	          ##println("A\t",log(resExact),"\t", log(resApproximate) + scaleTerm)
	          return log(resApproximate) + scaleTerm
	     end
	end

	#export logpdf
	function logpdf(dist::BivariateVonMisesNode, x::Array{Float64,1})
		return dist.k[1]*cos(x[1]-dist.mu[1]) + dist.k[2]*cos(x[2]-dist.mu[2]) - dist.k[3]*cos(x[1]-dist.mu[1]-x[2]+dist.mu[2]) - dist.logNormConst
	end

	#export pdf
	function pdf(dist::BivariateVonMisesNode, x::Array{Float64,1})
		return exp(logpdf(dist,x))
	end


	export add_bvm_point
	function add_bvm_point(dist::BivariateVonMisesNode, x::Array{Float64,1})
		if x[1] > -100.0 && x[2] > -100.0
			push!(dist.data, x)
			dist.ess2COS1 += cos(2.0*x[1])
			dist.ess2SIN1 += sin(2.0*x[1])
			dist.ess2COS2 += cos(2.0*x[2])
			dist.ess2SIN2 += sin(2.0*x[2])
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
	end

	function reset(dist::BivariateVonMisesNode)
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

	function loglik_and_prior(bv::BivariateVonMisesNode, data::Array{Array{Float64,1},1}, k::Array{Float64,1})
		try
			temp = BivariateVonMisesNode(Float64[k[1], k[2], k[3]], Float64[k[4], k[5]])
			
			bvmll = sum(Float64[logpdf(temp,x) for x in data]) 
			loglik = bvmll
			if length(data) < 200.0
				loglik += Distributions.logpdf(bv.kappa_prior, k[1]) 
				loglik += Distributions.logpdf(bv.kappa_prior, k[2])
				loglik += Distributions.logpdf(bv.kappa_prior, abs(k[3]))
			end
			#loglik = sum(Float64[logpdf(temp,x) for x in data])
			#println(loglik,"\t",bvmll,"\t", k)
			return loglik
		catch e
			#println("error", e)
			return -Inf
		end
	end

	function optimizeprior(bv::BivariateVonMisesNode)
		startstr = string("start value: ", bv.k, "\t", bv.mu)
		if bv.count == 0
			bv.k = Float64[1e-5,1e-5, 0.0]
		else
			kappalimit = 5000.0
			opt = Opt(:LN_COBYLA,5)
			data = filter(x -> x[1] > -100.0 && x[2] > -100.0, bv.data)
			localObjectiveFunction = ((params, grad) -> loglik_and_prior(bv, data, params))
			lower = ones(Float64, 5)*1e-5
			upper = ones(Float64, 5)*kappalimit
			lower[3] = -kappalimit
			lower[4] = Float64(0.0)
			lower[5] = Float64(0.0)
			upper[4] = Float64(2.0*pi)
			upper[5] = Float64(2.0*pi)
			lower_bounds!(opt, lower)
			upper_bounds!(opt, upper)
			xtol_rel!(opt,1e-3)
			maxeval!(opt, 2000)
			max_objective!(opt, localObjectiveFunction)

			maxk3 = 0.5*abs((bv.k[1]*bv.k[1] + bv.k[2]*bv.k[2])/(2.0*bv.k[1]))
			mink3 = -maxk3
			(minf,minx,ret) = optimize(opt, Float64[min(bv.k[1], kappalimit), min(bv.k[2], kappalimit), max(mink3, min(bv.k[3], maxk3)), mod2pi(bv.mu[1]), mod2pi(bv.mu[2])])
			bv.k[1] = minx[1]
			bv.k[2] = minx[2]
			bv.k[3] = minx[3]
			bv.mu[1] = minx[4]
			bv.mu[2] = minx[5]
			k = bv.k
			if abs(2.0*k[1]*k[3]) >= (k[1]*k[1] + k[2]*k[2])
				bv.k[3] = 0.0
			end
			#println(bv.k)
			bv.logNormConst = computeLogNormConst(bv.k)
			#println(ret)
		end		
		#println(startstr)
		#println("end value: ", bv.k, "\t", bv.mu)
	end

	#=
	function loglik_and_prior(bv::BivariateVonMisesNode, data::Array{Array{Float64,1},1}, k::Array{Float64,1})
		try
			temp = BivariateVonMisesNode(Float64[k[1],k[2], bv.k[3]], bv.mu)
			#loglik = sum(Float64[logpdf(temp,x) for x in data]) + Distributions.logpdf(bv.kappa_prior, k[1]) + Distributions.logpdf(bv.kappa_prior, k[2])
			loglik = sum(Float64[logpdf(temp,x) for x in data])
			println(loglik,"\t", k)
			return loglik
		catch e
			println("error", e)
			return -Inf
		end
	end

	function optimizeprior(bv::BivariateVonMisesNode)
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

	export estimate_bvm
	function estimate_bvm(dist::BivariateVonMisesNode)
		mu1 = pimod(atan(dist.essS1, dist.essC1))
		mu2 = pimod(atan(dist.essS2, dist.essC2))

		if dist.count < 10
			dist.k = Float64[0.5,0.5,0.0]
			dist.mu = Float64[mu1,mu2]		
		else
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
		end

		dist.logNormConst = computeLogNormConst(dist.k)

		reset(dist)
 	end


	
	export vonmisesrand
	function vonmisesrand(rng::AbstractRNG, mu::Float64, kappa::Float64)
		if kappa  < 1e-6
			return rand(rng)*2.0*pi
		else
			a = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa)
			b = (a - sqrt(2.0 * a)) / (2.0 * kappa)
			r = (1.0 + b * b) / (2.0 * b)
			c = 0.0
			f = 0.0
			U2 = 0.0

			cont = true
			while cont || ( ( c * (2.0 - c) - U2 <= 0 ) && ( c * exp(1.0 - c) - U2 < 0 ) )
				cont = false
				U1 = rand(rng)
				z = cos(pi*U1)
				f = (1.0 + r * z) / (r + z)
				c = kappa * (r - f)
				U2 = rand(rng)
			end

			U3 = rand(rng)

			if U3 - 0.5 > 0
				theta = mu + acos(f)
			else
				theta = mu - acos(f)
			end

			return mod2pi(theta)
		end
	end

	function vonmisesrand(rng::AbstractRNG, vonmises::VonMises)
		return vonmisesrand(rng, vonmises.μ, vonmises.κ)
	end

	export sample
 	function sample(rng::AbstractRNG, bvm::BivariateVonMisesNode)
 		phivm = VonMises(bvm.mu[1], bvm.k[1])
 		psivm = VonMises(bvm.mu[2], bvm.k[2])
 		currx = Float64[vonmisesrand(rng,phivm),vonmisesrand(rng,psivm)]
 		propx = copy(currx)
 		currll = logpdf(bvm, currx)
 		propll = currll 
 		#fout = open("mcmc.log","w")		
 		#println(fout, "iter\tmu1\tmu2")
 		for i=1:10000
 			if i % 3 == 0
 				propx[1] = vonmisesrand(rng, phivm)
 				propratio = Distributions.logpdf(phivm, currx[1]) - Distributions.logpdf(phivm, propx[1])
 				propll = logpdf(bvm, propx)
 				if exp(propll-currll+propratio) > rand(rng)
 					currll = propll
 					currx[1] = propx[1]
 					#println("$(i) 1 ", currx)
 				else
 					propx[1] = currx[1]
 				end
 			elseif i % 3 == 1
 				propx[2] = vonmisesrand(rng, psivm)
 				propratio = Distributions.logpdf(psivm, currx[2]) - Distributions.logpdf(psivm, propx[2])
 				propll = logpdf(bvm, propx)
 				if exp(propll-currll+propratio) > rand(rng)
 					currll = propll
 					currx[2] = propx[2]
 					#println("$(i) 2 ", currx)
 				else
 					propx[2] = currx[2]
 				end
 			else 
 				index = rand(rng,1:2)
 				propx[index] =  vonmisesrand(rng, currx[index], bvm.k[index]*1.0)
 				propll = logpdf(bvm, propx)
 				if exp(propll-currll) > rand(rng)
 					currll = propll
 					currx[index] = propx[index]
 					#println("$(i) 3 ", currx)
 				else
 					propx[index] = currx[index]
 				end
 			end
 			#println(fout, "$(i-1)\t$(currx[1])\t$(currx[2])")
 		end
 		#close(fout)
 		return currx
 	end

#=
 	function vonMisesSampler(rng::AbstractRNG, k::Float64, mean::Float64)
 		res = 0.0
 		a = 1.0 + sqrt(1.0 + 4.0*k*k)
 		b = (a - sqrt(2.0*a)) / (2.0*k)
 		r = (1.0 + b*b)/(2.0*b)
 		f = 0.0
 		while true
 			U1 = rand(rng)
 			z = cos(pi*U1)
 			f = (1.0 + r*z)/(r + z)
 			c = k * (r - f)
			U2 = rand(rng)
		 	if (((c*(2.0-c) - U2) > 0.0) || ((log(c/U2) + 1.0 - c >= 0.0)))
		 		break
          	end
 		end
 		U3 = rand(rng)
 		if U3 > 0.5
 			res = mod(acos(f)+mean, 2.0*pi)
 		else
 			res = mod(-acos(f)+mean, 2.0*pi)
 		end
 		return res
 	end

 	function sample(rng::AbstractRNG, bvm::BivariateVonMisesNode)

	     double k13 = sqrt(k[1]*k[1] + k[3]*k[3] - 2.0*k[1]*k[3]*cos(psi));
	     double psi_mu = atan((-k3*sin(psi))/(k1 - k3*cos(psi)));
	     double phi = vonMisesSampler(k13, psi_mu, rg);

	     // Angles are in the interval [-pi, pi]. Add 3*pi to bring them
	     // to: [2*pi, 4*pi] which with mu values in [-pi, pi] brings it
	     // to [pi, 3pi], which can then readily be transformed back to [-pi, pi]
	     psi = fmod(psi + (3*M_PI) + mu2, 2*M_PI) - M_PI;
	     phi = fmod(phi + (3*M_PI) + mu1, 2*M_PI) - M_PI;

	     rVal.resize(2);
	     rVal[0] = phi;
	     rVal[1] = psi;
 	end

 	  double vonMisesSampler(double k, double mean, RandomGen* rg) {
     if (mean < -M_PI || mean > M_PI) {
          fprintf(stderr, "vonMises Error: mean must be in the interval (-pi,pi). Mean=%f\n", mean);
     }
     double res;
     double a = 1.0 + sqrt(1.0 + 4.0*k*k);
     double b = (a - sqrt(2*a)) / (2*k);
     double r = (1.0 + b*b)/(2*b);
     double f;
     while(true) {
          // double U1 = RandomLib::Random::Global.FixedN();
          double U1 = rg->get_rand();
          double z = cos(M_PI * U1);
          f = (1.0 + r*z)/(r + z);
          double c = k * (r - f);
          // double U2 = RandomLib::Random::Global.FixedN();
          double U2 = rg->get_rand();
          if (((c*(2.0-c) - U2) > 0.0) || ((log(c/U2) + 1.0 - c >= 0.0))){
               break;           // accept
          }
     }
     // double U3 = RandomLib::Random::Global.FixedN();
     double U3 = rg->get_rand();
     if (U3 > 0.5) {
          res = fmod(acos(f)+mean, 2*M_PI);
     } else {
          res = fmod(-acos(f)+mean, 2*M_PI);
     }
     return res;
}=#
end
#=
using Random
bvm = BivariateVonMises.BivariateVonMisesNode(Float64[400.0, 400.0, 200.0], Float64[pi,pi])
rng = MersenneTwister(1)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)
currx =  BivariateVonMises.sample(rng, bvm)
println(currx)=#