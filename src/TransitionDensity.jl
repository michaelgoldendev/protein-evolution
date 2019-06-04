using Distributions
using Random
push!(LOAD_PATH,@__DIR__)
using CommonUtils

function indextoangle(index::Int, numcats::Int=100)
	delta = 2.0*pi/numcats
	deltadiv2 = delta/2.0
	return index*delta - deltadiv2
end

function angletoindex(angle::Float64, numcats::Int=100)
	delta = 2.0*pi/numcats
	deltadiv2 = delta/2.0
	return Int(floor((mod2pi(angle) + deltadiv2)/delta)) + 1	
end

function sampleconvolvedvonmises(rng::AbstractRNG, v1::VonMises, v2::VonMises, maxiters::Int=5000)
	m = logpdf(v1, v1.μ) + logpdf(v2, v2.μ)
	for i=1:maxiters
		theta = 2.0*pi*rand(rng) - pi
		ll = logpdf(v1, theta) + logpdf(v2, theta)
		if ll > m + log(rand(rng))
			println(i)
			return theta
		end
	end
	println("X")

	len = 1000
	delta = 2.0*pi/len
	x = zeros(Float64,len)
	x[1] = delta/2.0
	for i=2:len
		x[i] = x[i-1]+delta
	end
	probs = logpdf.(v1, x) + logpdf.(v2, x)
	probs = exp.(probs .- maximum(probs))
	i1 = CommonUtils.sample(rng, probs)
	return x[i1] - delta/2.0 + delta*rand(rng)
end

function sampleconvolvedvonmises(rng::AbstractRNG, v1::VonMises, v2::VonMises, v3::VonMises, maxiters::Int=5000)
	m = logpdf(v1, v1.μ) + logpdf(v2, v2.μ)  + logpdf(v3, v3.μ)
	for i=1:maxiters
		theta = 2.0*pi*rand(rng) - pi
		ll = logpdf(v1, theta) + logpdf(v2, theta) + logpdf(v3, theta)
		if ll > m + log(rand(rng))
			println(i)
			return theta
		end
	end	
	println("X")

	len = 1000
	delta = 2.0*pi/len
	x = zeros(Float64,len)
	x[1] = delta/2.0
	for i=2:len
		x[i] = x[i-1]+delta
	end
	probs = logpdf.(v1, x) + logpdf.(v2, x) + logpdf.(v3, x)
	probs = exp.(probs .- maximum(probs))
	i1 = CommonUtils.sample(rng, probs)
	return x[i1] - delta/2.0 + delta*rand(rng)
end

function sampleinternalnode(rng::AbstractRNG, distphi1::VonMises, distpsi1::VonMises, distphi2::VonMises, distpsi2::VonMises)
	return sampleconvolvedvonmises(rng, distphi1, distphi2), sampleconvolvedvonmises(rng, distpsi1, distpsi2)
end

function sampleinternalnode(rng::AbstractRNG, distphi1::VonMises, distpsi1::VonMises, distphi2::VonMises, distpsi2::VonMises, distphi3::VonMises, distpsi3::VonMises)
	return sampleconvolvedvonmises(rng, distphi1, distphi2, distphi3), sampleconvolvedvonmises(rng, distpsi1, distpsi2, distpsi3)
end

rng = MersenneTwister(871903801543)
distphi1 = VonMises(2.0,100.0)
distpsi1 = VonMises(2.5,100.0)
distphi2 = VonMises(2.1,100.0)
distpsi2 = VonMises(3.0,100.0)
distphi3 = VonMises(2.5,100.0)
distpsi3 = VonMises(3.2,100.0)

phiarray = Float64[]
psiarray = Float64[]
for i=1:10000
	phi,psi = sampleinternalnode(rng, distphi1, distpsi1, distphi2, distpsi2, distphi3, distpsi3)
	push!(phiarray, phi)
	push!(psiarray, psi)
end
println(mean(phiarray))
println(mean(psiarray))