aminoacids = "ACDEFGHIKLMNPQRSTVWY"
aamapping = Dict()
alphabet = length(aminoacids)
for a=1:alphabet
  aamapping[aminoacids[a]] = a
end

fastaamapping = zeros(Int,26)
for (index,aa) in enumerate(aminoacids)
  fastaamapping[Int(aa)-64] = index
end


function kmertoindex(kmer::AbstractString)
  len = length(kmer)
  nuc = fastaamapping[Int(kmer[1])-64]
  if nuc == 0
    return -1
  end
  index = nuc-1
  for pos=2:len
    index *= alphabet
    nuc = fastaamapping[Int(kmer[pos])-64]
    if nuc == 0
      return -1
    end
    index += nuc-1
  end
  return index+1
end

function freqvector(sequence::AbstractString, k::Int, normalised::Bool=true)
  sequence2 = replace(sequence,"-" => "")
  f = zeros(Float64,alphabet^k)
  for startpos=1:length(sequence2)-k+1
    endpos = startpos+k-1
    kmer = sequence2[startpos:endpos]
    index = kmertoindex(kmer)
    if index >= 0
      f[index] += 1
    end
  end

  if normalised
    return f/sum(f)
  else
    return f
  end
end

function freqvectorint(sequence::AbstractString, k::Int)
  f = zeros(Int,alphabet^k)
  for startpos=1:length(sequence)-k+1
    endpos = startpos+k-1
    kmer = sequence[startpos:endpos]
    index = kmertoindex(kmer)
    if index >= 0
      f[index] += 1
    end
  end
  return f
end



function euclidean(f1::Array{Float64,1},f2::Array{Float64,1})
  return sqrt(mean([(v1-v2)^2.0 for (v1,v2) in zip(f1,f2)]))
end

function bintovec(bin::Array{Int, 1})
  sel = Int[]
  for i=1:length(bin)
    if bin[i] > 0
      push!(sel, i)
    end
  end
  return sel
end

function vectobin(sel::Array{Int,1}, numseqs::Int)
  bin = zeros(Int,numseqs)
  for s in sel
    bin[s] = 1
  end
  return bin
end