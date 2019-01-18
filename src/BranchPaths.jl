module BranchPaths	

	export BranchPath
	mutable struct BranchPath
		paths::Array{Array{Int,1},1}
		times::Array{Array{Float64,1},1}

		Rmatrices::Array{Array{Float64,2},1}
		Pmatrices::Array{Array{Float64,2},1}
		vs::Array{Array{Float64,1},1}
		time::Array{Float64,1}

		function BranchPath(len::Int)
			paths = Array{Int,1}[]
			times = Array{Float64,1}[]
			for i=1:len
				push!(paths, Int[])
				push!(times, Float64[])
			end
			return new(paths,times, Array{Float64,2}[], Array{Float64,2}[], Array{Float64,1}[], Float64[])
		end

		function BranchPath(paths::Array{Array{Int,1},1}, times::Array{Array{Float64,1},1})
			new(paths,times, Array{Float64,2}[], Array{Float64,2}[], Array{Float64,1}[], Float64[])
		end
	end

	export BranchPathIterator
	mutable struct BranchPathIterator
		branch::BranchPath
		indices::Array{Int,1}
		counters::Array{Int,1}
		prevstates::Array{Int,1}
		prevtime::Float64
		currstates::Array{Int,1}
		currtime::Float64
		mincol::Int
		minindex::Int
		mintime::Float64
		col::Int
		done::Bool
		count::Int

		function BranchPathIterator(branch::BranchPath, indices::Array{Int,1})
			counters = ones(Int,length(indices))*2
			prevstates = Int[branch.paths[index][1] for index in indices]
			currstates = Int[branch.paths[index][1] for index in indices]
			mincol = -1
			minindex = -1
			mintime = 0.0
			col = 1
			done = false
			return new(branch,indices,counters,prevstates,0.0,currstates,0.0,mincol,minindex,mintime,col,done,0)
		end
	end

	export iterate
	function Base.iterate(iter::BranchPathIterator, state=(nothing, 0))
       prevelement, iter.count = deepcopy(state)

	   if prevelement == nothing
		   fill!(iter.counters, 2)
		   iter.done = false
		   iter.mincol = -1
		   prevelement = deepcopy((iter.prevstates, iter.prevtime, iter.currstates, iter.currtime, iter.mincol))
		   prevelement, iter.count = deepcopy(next(iter,iter.count))
	   end

       if iter.done && prevelement != nothing && prevelement[5] == -2
		    return nothing
       elseif iter.done
		   prevelement = deepcopy((prevelement[1],prevelement[2],prevelement[3],prevelement[4],-2))
		   return (prevelement, (prevelement,iter.count+1))
	   end
	   nextelement, nextcount = deepcopy(next(iter,iter.count))
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

	export peek
	function peek(it::BranchPathIterator)
		#=
		mincol = -1
		minindex = -1
		mintime = 0.0		
		col = 1
		for index in it.indices
			if it.counters[col] <= length(it.branch.paths[index])
				time = it.branch.times[index][it.counters[col]]
				if mincol == -1 || time < mintime
					mintime = time
					mincol = col
					minindex = index
				end
			end
			col += 1
		end
		return it.mintime,it.mincol=#
		mintime = Inf
		for (col,index) in enumerate(it.indices)
			if it.counters[col] <= length(it.branch.times[index])
				time = it.branch.times[index][it.counters[col]]
				if time < mintime
					mintime = time
				end
			end
		end
		return mintime
		#return it.prevtime,it.currtime,it.mintime,it.mincol
	end

	
	export MultiBranchPathIterator
	mutable struct MultiBranchPathIterator
		branchpathiterators::Array{BranchPathIterator,1}
		branchpathindex::Int
		mintime::Float64
		prevtime::Float64
		currtime::Float64
		done::Bool
		count::Int
		mincol::Int

		function MultiBranchPathIterator(branchpathiterators::Array{BranchPathIterator,1})
			return new(branchpathiterators,1,0.0,0.0,0.0,false,0,0)
		end


	end

	export iterate
	function Base.iterate(multi_iter::MultiBranchPathIterator,state=(nothing,0))
		prevelement, multi_iter.count = state
		multi_iter.prevtime = multi_iter.currtime

		miniterator = 1
		mintime = Inf
		for (itindex, it) in enumerate(multi_iter.branchpathiterators)
			peekedtime = peek(it)
			if peekedtime < mintime && !it.done
				mintime = peekedtime
				miniterator = itindex
			end
		end
		curriter = multi_iter.branchpathiterators[miniterator]
		multi_iter.branchpathindex = miniterator

		if multi_iter.count == 0
			for i=1:length(multi_iter.branchpathiterators)
				if i != miniterator
					next(multi_iter.branchpathiterators[i],0)
				end
			end
		end
		#=
		multi_iter.currtime = curriter.currtime
		println(multi_iter.prevtime,"\t",multi_iter.currtime)
		if curriter.mincol == -1
			#multi_iter.done = true
			#return nothing
			alldone = true
			for it in multi_iter.branchpathiterators
				if !it.done 
					alldone = false
				end
			end
			if alldone
				multi_iter.done = true
				return nothing
			end
		end=#

		
		#println(miniterator,"\t", next(curriter,curriter.count)[1])
		next(curriter,curriter.count)
		multi_iter.currtime = curriter.currtime
		multi_iter.mincol = curriter.mincol
		#println(multi_iter.prevtime,"\t",multi_iter.currtime)
		if multi_iter.count == 0 && curriter.mincol == -1	
			multi_iter.prevtime = 0.0
			multi_iter.currtime = 1.0
		end
		if multi_iter.count != 0 && curriter.mincol == -1			
			multi_iter.done = true
			return nothing
		end
		return (prevelement, (curriter,multi_iter.count+1))
	end

	#=
	export iterate
	function Base.iterate(iter::MultiBranchPathIterator, state=(nothing, 0))
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
   end=#

	
end

#=
paths = Array{Int,1}[]
times = Array{Float64,1}[]
push!(paths, Int[0,2,4,6,8])
push!(times, Float64[0.0,0.1,0.2,0.4,0.7])
push!(paths, Int[1,3,5,7,9])
push!(times, Float64[0.0,0.11,0.21,0.41,0.71])
branchpath = BranchPaths.BranchPath(paths,times)

branchiterator = BranchPaths.BranchPathIterator(branchpath,Int[1,2])
for (prevstates, prevtime, currstates, currtime, mincol) in branchiterator
	println(prevstates,"\t", prevtime,"\t", currstates,"\t", currtime,"\t", mincol)
end=#
