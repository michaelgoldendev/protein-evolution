using Random

push!(LOAD_PATH,string(@__DIR__,"/../../MolecularEvolution/src/"))
using MolecularEvolution

push!(LOAD_PATH,@__DIR__)
using LG
using CTMCs

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

println(gettindex(0.1))
println(gettindex(0.3))
println(gettindex(0.2))
println(gettindex(0.5))
println(gettindex(0.6))
println(gettindex(0.06))
println(gettindex(0.34))
println(gettindex(0.01))
println(gettindex(0.001))

aapathdict = Dict{Tuple{Int,Int,Int,Int},Array{Any,1}}()
function putaapath(h::Int, t::Float64, path::Array{Int,1}, time::Array{Float64,1})
	tindex = gettindex(t)
	key = (h,aa1,aa2,t)
	ls = get(aapathdict, key, [])
	push!(ls, (t,path,time))
	aapathdict[key] = ls
end

function getaapath(rng::AbstractRNG, h::Int, tin::Float64, R::Array{Float64,2}, aa1::Int, aa2::Int)
	tindex, t = gettindex(tin)
	key = (h,aa1,aa2,tindex)
	if haskey(aapathdict, key) && length(aapathdict[key]) > 0
		t,path,time = pop!(aapathdict[key])
		return t,path,time
	else
		cache = []
		path,time = modifiedrejectionsampling3(rng,aapathdict,key, R,t, aa1, aa2)
		return t,path,time
	end	
end

rng = MersenneTwister(1)
println(getaapath(rng, 1, 0.1, LGmatrix, 1, 2))
println(getaapath(rng, 1, 0.08, LGmatrix, 1, 2))