aminoacids = "ACDEFGHIKLMNPQRSTVWY"
aamapping = Dict()
alphabet = length(aminoacids)
for a=1:alphabet
  aamapping[aminoacids[a]] = a
end


function kmertoindex(kmer::AbstractString)
  len = length(kmer)
  nuc = get(aamapping,kmer[1], 0)
  if nuc == 0
    return -1
  end
  index = nuc-1
  for pos=2:len
    index *= alphabet
    nuc = get(aamapping,kmer[pos], 0)
    if nuc == 0
      return -1
    end
    index += nuc-1
  end
  return index+1
end

function freqvector(sequence::AbstractString, k::Int)
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
  return f/sum(f)
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