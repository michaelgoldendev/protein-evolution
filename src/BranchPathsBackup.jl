module BranchPaths
	export BranchState
	mutable struct BranchState
		h::Int
		aa::Int
		phi::Float64
		psi::Float64

		function BranchState(h::Int,aa::Int,phi::Float64,psi::Float64)
			new(h,aa,phi,psi)
		end
	end

	export BranchPath
	mutable struct BranchPath
		paths::Array{Array{BranchState,1},1}
		times::Array{Array{Float64,1},1}

		function BranchPath(len::Int)
			paths = Array{BranchState,1}[]
			times = Array{Float64,1}[]
			for i=1:len
				push!(paths, BranchState[])
				push!(times, Float64[])
			end
			return new(paths,times)
		end

		function BranchPath(paths::Array{Array{BranchState,1},1}, times::Array{Array{Float64,1},1})
			new(paths,times)
		end
	end

	export BranchPathIterator
	mutable struct BranchPathIterator
		branch::BranchPath
		indices::Array{Int,1}
		counters::Array{Int,1}
		prevstates::Array{BranchState,1}
		prevtime::Float64
		currstates::Array{BranchState,1}
		currtime::Float64
		mincol::Int
		minindex::Int
		mintime::Float64
		col::Int
		done::Bool

		function BranchPathIterator(branch::BranchPath, indices::Array{Int,1})
			counters = ones(Int,length(indices))*2
			prevstates = BranchState[branch.paths[index][1] for index in indices]
			currstates = BranchState[branch.paths[index][1] for index in indices]
			mincol = -1
			minindex = -1
			mintime = 0.0
			col = 1
			done = false
			return new(branch,indices,counters,prevstates,0.0,currstates,0.0,mincol,minindex,mintime,col,done)
		end
	end

	export iterate
	function Base.iterate(iter::BranchPathIterator, state=(nothing, 0))
       prevelement, count = deepcopy(state)

	   if prevelement == nothing
		   fill!(iter.counters, 2)
		   iter.done = false
		   iter.mincol = -1
		   prevelement = deepcopy((iter.prevstates, iter.prevtime, iter.currstates, iter.currtime, iter.mincol))
		   prevelement, count = deepcopy(next(iter,count))
	   end

       if iter.done && prevelement != nothing && prevelement[5] == -2
		    return nothing
       elseif iter.done
		   prevelement = deepcopy((prevelement[1],prevelement[2],prevelement[3],prevelement[4],-2))
		   return (prevelement, (prevelement,count+1))
	   end
	   nextelement, nextcount = deepcopy(next(iter,count))
	   ret = (prevelement, (nextelement,nextcount))
       return ret
   end

	function next(it::BranchPathIterator,state)
		if it.mincol != -1
			it.prevstates[it.mincol] = it.branch.paths[it.minindex][it.counters[it.mincol]-1]
			it.prevtime = it.mintime
		end
		it.mincol = -1
		it.minindex = -1
		it.mintime = 0.0
		it.col = 1
		for index in it.indices
			if it.counters[it.col] <= length(it.branch.paths[index])
				time = it.branch.times[index][it.counters[it.col]]
				if it.mincol == -1 || time < it.mintime
					it.mincol = it.col
					it.minindex = index
					it.mintime = time
				end
			end
			it.col += 1
		end

		prevcounters = copy(it.counters)
		if it.mincol == -1
			it.done = true
			it.currtime = 1.0
		else
			it.counters[it.mincol] += 1
		end

		if it.mincol != -1
			it.currstates[it.mincol] = it.branch.paths[it.minindex][it.counters[it.mincol]-1]
			it.currtime = it.mintime
		end
		return (it.prevstates, it.prevtime, it.currstates, it.currtime, it.mincol), state+1
	end
end

#=
paths = Array{BranchPaths.BranchState,1}[]
times = Array{Float64,1}[]
p1 = BranchPaths.BranchState[]
push!(p1, BranchPaths.BranchState(1,1,0.0,0.0))
push!(p1, BranchPaths.BranchState(2,2,0.0,0.0))
push!(p1, BranchPaths.BranchState(3,3,0.0,0.0))
push!(paths,p1)
push!(times, Float64[0.0,0.1,0.2])

p2 = BranchPaths.BranchState[]
push!(p2, BranchPaths.BranchState(4,4,0.0,0.0))
push!(p2, BranchPaths.BranchState(5,5,0.0,0.0))
push!(p2, BranchPaths.BranchState(6,6,0.0,0.0))
push!(paths,p2)
push!(times, Float64[0.0,0.11,0.21])

branchpath = BranchPaths.BranchPath(paths,times)

branchiterator = BranchPaths.BranchPathIterator(branchpath,Int[1,2])
for (prevstates, prevtime, currstates, currtime, mincol) in branchiterator
	println(prevstates,"\t", prevtime,"\t", currstates,"\t", currtime,"\t", mincol)
end=#
