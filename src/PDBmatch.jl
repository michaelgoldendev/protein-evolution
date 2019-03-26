using FastaIO

function issubsequence(query::String, target::String)
	if length(query) > length(target)
		return false,0
	end
	b = 0
	len = length(target)
	numinserts = 0
	for q in query
		b += 1
		if b > len
			return false, numinserts
		end
		isinsert = false
		while q != target[b]
			isinsert = true
			b += 1
			if b > len
				return false, numinserts
			end
		end

		if isinsert
			numinserts += 1
		end
	end
	return true, numinserts
end


#fastafile = abspath("../data/influenza_a/HA/selection3.fasta")
#fastafile = abspath("../data/hiv/curated6.fasta")
#fastafile = abspath("../data/test_data/hiv_pol_selection.fasta")
#fastafile = abspath("../data/test_data/maise_streak_virus_coat_protein_selection.fasta")
fastafile = abspath("../data/test_data/westnile_dengue_selection.fasta")

sequences = AbstractString[]
names = AbstractString[]
FastaIO.FastaReader(fastafile) do fr
    for (desc, seq) in fr
        push!(names,desc)
        push!(sequences, replace(seq, "-" => ""))
    end
end

open("../data/pdbsequences.fasta") do file
	name = ""
	seq = ""
	count = 0
    for ln in eachline(file)
    	if startswith(ln,">")
    		name = ln
    		seq = ""
    	else
    		count += 1
    		if count % 10000 == 0
    			println(count)
    		end
    		seq = ln

    		for (query,queryname) in zip(sequences, names)
			 	if length(seq)/length(query) > 0.1 && issubsequence(seq, query)[1]
			 		subsequence, inserts = issubsequence(seq, query)
			 		if inserts <= 5
				 		println("X INSERTS ", inserts)
				 		println(length(seq)/length(query))
			        	println(name)
			        	println(seq)
			        	println()
			        end
		        end

		        if length(query)/length(seq) > 0.1 && issubsequence(query, seq)[1]
			 		subsequence, inserts = issubsequence(query, seq)
			 		if inserts <= 5
				 		println("Y INSERTS ", inserts)
				 		println(length(query)/length(seq))
			        	println(name)
			        	println(seq)
			        	println()
			        end
		        end
		    end
    	end
       
    end
end
