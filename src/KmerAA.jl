push!(LOAD_PATH,string(@__DIR__,"/../../MolecularEvolution/src/"))
using MolecularEvolution
using Nullables
using Base.Filesystem
using FastaIO
using Statistics
push!(LOAD_PATH,@__DIR__)
using Binaries
using CommonUtils

aminoacids = "ACDEFGHIKLMNPQRSTVWY"
nucmapping = Dict()
alphabet = length(aminoacids)
for a=1:alphabet
  nucmapping[aminoacids[a]] = a
end

function binarize!(tree::TreeNode)
    nodes = getnodelist(tree)
    counter = 0
    for n in nodes
        while length(n.children) > 2
            c1 = pop!(n.children)
            c2 = pop!(n.children)
            counter +=1 
            push!(n.children, TreeNode(0.0, "binarized_$counter"))
            n.children[end].children = [c1,c2]
            c1.parent = Nullable{TreeNode}(n.children[end])
            c2.parent = Nullable{TreeNode}(n.children[end])
            n.children[end].parent = Nullable{TreeNode}(n)
        end
    end
end

function nodedistances(nodelist::Array{TreeNode,1})
  for (index,node) in enumerate(nodelist)
    node.nodeindex = index
  end

  distmatrix = ones(Float64, length(nodelist), length(nodelist))*Inf
  for node in nodelist
    distmatrix[node.nodeindex,node.nodeindex] = 0.0
    for child in node
      distmatrix[node.nodeindex,child.nodeindex] = child.branchlength
      distmatrix[child.nodeindex,node.nodeindex] = child.branchlength
    end
  end

  mindistmatrix = copy(distmatrix)
  for n1 in nodelist
    for v1 in nodelist
      for v2 in nodelist
        mindistmatrix[n1.nodeindex,v2.nodeindex] = min(mindistmatrix[n1.nodeindex,v2.nodeindex], mindistmatrix[n1.nodeindex,v1.nodeindex]+mindistmatrix[v1.nodeindex,v2.nodeindex])
      end
    end
  end

  return mindistmatrix
end

#=
function root(node::TreeNode)
  parents = []
  current = node
  while !isnull(current.parent)
    current = get(current.parent)
    push!(parents, current)
  end
  println([p.nodeindex for p in parents])

  bottom = parents[end]
  while length(parents) > 0
    pop!(parents)
    println("A",[b.nodeindex for b in bottom])
    println("B",bottom.nodeindex,"\t",parents[end].nodeindex)   
    exit()
    leftchild = bottom.children[1]
    push!(leftchild.children,bottom)    
    bottom.children = bottom.children[2:end]


    bottom.parent = leftchild
    newroot = TreeNode(0.1, "dummy$(length(parents))")
    push!(newroot.children,leftchild.children[1])
    push!(newroot.children,leftchild)
    leftchild.children[1].parent = Nullable{TreeNode}(newroot)
    leftchild.parent = Nullable{TreeNode}(newroot)
    leftchild.children = leftchild.children[2:end]
    bottom = newroot
  end
  return bottom
end=#

function reorient(node, new_parent, new_branch_length)
    newchildren = TreeNode[]
    for c in node
      if c != new_parent
        push!(newchildren,c)
      end
    end
    node.children = newchildren
    # If this node has a parent, reorient the parent.
    if !isnull(node.parent)
        parent = node.parent.value
        reorient(parent, node, node.branchlength)
        # Add the parent as a child
        push!(node.children, parent)
    end
    # Set then new parent as the parent, with the new branch length
    node.parent = Nullable{TreeNode}(new_parent)
    node.branchlength = new_branch_length
end

function reroot(child_node, dist_above_child=(child_node.branchlength/2))
    if dist_above_child > child_node.branchlength
        print("This isn't going to work")
    end
        
    # Remembering stuff
    dist_below_parent = child_node.branchlength - dist_above_child
    old_parent = child_node.parent.value
        
    new_root = TreeNode(0.0,"root")
    child_node.branchlength = dist_above_child
    reorient(old_parent, child_node, dist_below_parent)
    new_root.children = [child_node, old_parent]
    return new_root
end


function midpoint_root(nodelistin::Array{TreeNode,1})
  nodelist = deepcopy(nodelistin)
  for node in nodelist
    if startswith(node.name,":") || node.name == ""
      node.name = string("node",node.nodeindex)
    end
  end

  mindistmatrix = nodedistances(nodelist)
  minnodetotipdistance = Inf
  minnodeindex = 0
  for nonleaf in nodelist
    if !isleafnode(nonleaf)     
      nodetotipdistance = 0.0
      vals = Float64[]
      for leaf in nodelist
        if isleafnode(leaf)
          nodetotipdistance += mindistmatrix[nonleaf.nodeindex, leaf.nodeindex]
          push!(vals, mindistmatrix[nonleaf.nodeindex, leaf.nodeindex])
        end
      end
      stdeviation = std(vals)
      if stdeviation < minnodetotipdistance
        minnodetotipdistance = stdeviation
        minnodeindex = nonleaf.nodeindex
      end
    end
  end

  newroot = reroot(nodelist[minnodeindex], nodelist[minnodeindex].branchlength/2.0)
  nodelist = getnodelist(newroot)
  for node in nodelist
    if length(node.children) == 1
      node.children = node.children[1].children
      for c in node.children
        c.parent = Nullable{TreeNode}(node)
      end
    end
  end

  return newroot
  #root(newroot)

end

function kmertoindex(kmer::AbstractString)
  len = length(kmer)
  nuc = get(nucmapping,kmer[1], 0)
  if nuc == 0
    return -1
  end
  index = nuc-1
  for pos=2:len
    index *= alphabet
    nuc = get(nucmapping,kmer[pos], 0)
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


function euclidean(f1::Array{Float64,1},f2::Array{Float64,1})
  return sqrt(mean([(v1-v2)^2.0 for (v1,v2) in zip(f1,f2)]))
end

function score(distmatrix::Array{Float64,2}, sel::Array{Int,1}, z::Int=0)

  len = length(sel)
  #=
  d = 0.0
  for i=1:len
    for j in sel
      d += distmatrix[i,j]
    end
  end

  if z != 0 && !(z in sel)
    for i=1:len
      d += distmatrix[i,z]
    end
  end
  return d=#
  dist = 0.0
  for i=1:len
    d = Inf
    for j in sel
      if i != j
        d = min(d, distmatrix[i,j])
      end
    end
    dist += d
  end
  return dist
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

function greedyselection(distmatrix::Array{Float64,2}, seln::Int)
  len = size(distmatrix,1)
  sel = Int[rand(1:len)]

  maxscore = 0.0
  while length(sel) < seln
    maxindex = 1
    maxscore = 0.0
    for i=1:len
      if !(i in sel) && score(distmatrix,sel,i) > maxscore
        maxindex = i
        maxscore = score(distmatrix,sel,i)
      end
    end
    push!(sel, maxindex)
  end

  return sel, maxscore
end

function swap(distmatrix::Array{Float64,2}, bin::Array{Int,1}, s::Int, z::Int)
  sel = bintovec(bin)
  #selset = IntSet(sel)
  comp = Int[]
  for i=1:length(bin)
    if bin[i] == 0
     push!(comp, i)
    end
  end
  #compset = IntSet(comp)

  maxsel = sel
  maxcomp = comp
  maxscore = score(distmatrix,sel,z)
  #maxscore = 0.0
  println(maxscore)
  for k=1:2000000
    for j=1:s
        i1 = rand(1:length(sel))
        i2 = rand(1:length(comp))
        v1 = sel[i1]
        v2 = comp[i2]
        sel[i1] = v2
        comp[i2] = v1
    end
    if score(distmatrix,sel,z) > maxscore
      maxscore = score(distmatrix,sel,z)
      maxsel = copy(sel)
      maxcomp = copy(comp)
      println(k,"\t",maxscore)
    else
      sel = copy(maxsel)
      comp = copy(maxcomp)
    end
  end

  return maxsel
end

mutable struct BranchSupportData <: NodeData
  branchsupport::Float64

  BranchSupportData(branchsupport::Float64) = new(branchsupport) 
end

function getfloatvalue(x::AbstractString, default::Float64=0.0)
  if x == ""
    return 1.0
  end
  v = tryparse(Float64, x)
  if v == nothing
    return default
  else
    return v
  end
end

function getcommonancestor(n1::TreeNode, n2::TreeNode)
  i1 = Int[n1.nodeindex]
  i2 = Int[n2.nodeindex]
  if n1.nodeindex == n2.nodeindex
    return n1
  end

  p1 = n1.parent
  while !isnull(p1)
    push!(i1, get(p1).nodeindex)
    p1 = get(p1).parent
  end

  p2 = n2.parent
  while !isnull(p2)
    p = get(p2)
    if p.nodeindex in i1
      return p
    end
    #push!(i2, p.nodeindex)
    p2 = p.parent
  end
 
end

function getcommonancestor2(n1::TreeNode, n2::TreeNode)
  i1 = Int[n1.nodeindex]
  i2 = Int[n2.nodeindex]
  if n1.nodeindex == n2.nodeindex
    return n1
  end

  p1 = n1.parent
  while !isnull(p1)
    push!(i1, get(p1).nodeindex)
    p1 = get(p1).parent
  end

  ret = nothing
  p2 = n2.parent
  while !isnull(p2)
    p = get(p2)
    if p.nodeindex in i1
      if ret == nothing
        ret = p
      end
    end
    push!(i2, p.nodeindex)
    p2 = p.parent
  end
  println(i1,"\t",i2,"\t",ret.nodeindex)
  return ret
end

function evaluatesupport(supportmatrix::Array{Float64,2}, binvec::Array{Int,1})
  minsupport = 1.0
  for i=1:length(binvec)
    for j=i+1:length(binvec)
      if binvec[i] == 1 && binvec[j] == 1
        minsupport = min(minsupport,supportmatrix[i,j])
      end
    end
  end
  return minsupport
end

function gettotalbranchlength(nodelist::Array{TreeNode,1}, binvec::Array{Int,1})
  branchlength = 0.0
  used = zeros(Int,length(nodelist))
  for node in nodelist
    if node.seqindex > 0 && binvec[node.seqindex] == 1
      
      if used[node.nodeindex] == 0
        branchlength += node.branchlength
        used[node.nodeindex] = 1
      end

      p1 = node.parent
      while !isnull(p1)
        p = get(p1)
        if used[p.nodeindex] == 0
          branchlength += p.branchlength
          used[p.nodeindex] = 1
        end
        p1 = p.parent
      end

    end
  end
  return branchlength
end

function removeleaves(nodelist::Array{TreeNode,1}, binvec::Array{Int,1})
  root = nodelist[1]
  for node in nodelist
    if node.seqindex > 0 && binvec[node.seqindex] == 0
      parent = get(node.parent)
      deleteat!(parent.children, findfirst(x -> x.nodeindex == node.nodeindex, parent.children))

      if !isnull(parent.parent) && length(parent.children) == 1
        parentparent = get(parent.parent)
        deletednode = deleteat!(parentparent.children, findfirst(x -> x.nodeindex == parent.nodeindex, parentparent.children))
        parent.children[1].branchlength += parent.branchlength
        push!(parentparent.children, parent.children[1])
        parent.children[1].parent = Nullable{TreeNode}(parentparent)
      else
        root = parent
      end
    end
  end

  if length(root.children) == 1 && !isleafnode(root.children[1])
    root.children = root.children[1].children
    for c in root.children
      c.parent = Nullable{TreeNode}(root)
    end
  end

  return root
end

function removeleaves(nodelist::Array{TreeNode,1}, nodenames::Array{String,1})
  root = nodelist[1]
  for node in nodelist
    if node.name in nodenames && isleafnode(node)
      parent = get(node.parent)
      deleteat!(parent.children, findfirst(x -> x.nodeindex == node.nodeindex, parent.children))

      if !isnull(parent.parent) && length(parent.children) == 1
        parentparent = get(parent.parent)
        deletednode = deleteat!(parentparent.children, findfirst(x -> x.nodeindex == parent.nodeindex, parentparent.children))
        parent.children[1].branchlength += parent.branchlength
        push!(parentparent.children, parent.children[1])
        parent.children[1].parent = Nullable{TreeNode}(parentparent)
      else
        root = parent
      end
    end
  end
  
  if length(root.children) == 1 && !isleafnode(root.children[1])
    root.children = root.children[1].children
    for c in root.children
      c.parent = Nullable{TreeNode}(root)
    end
  end

  return root
end

#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\westnile_virus_polyprotein\\ncbi.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\chikungunya_virus_nonstructural\\ncbi.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rabbit_hemorrhagic_disease_virus_VP60\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\hepatitis_c_polyprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\betacoronavirus_nucleoprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\lassavirus_nucleoprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\lassavirus_polymerase\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\enterovirus_a_2C\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\enterovirus_a_3D\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rhinovirus_a_VP1\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rhinovirus_a_VP2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rhinovirus_a_VP3\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\ebolavirus_GP\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\ebolavirus_NP\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\ebolavirus_L\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\myxoma_virus_m150R\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\betacoronavirus_rdrp\\viprbrc.fasta"

#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_herpesvirus3_glycoproteinE\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\puumala_orthohantavirus_N\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\puumala_orthohantavirus_S\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\dengue_virus_polyprotein\\viprbrc.fasta"

#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\orthohepevirus_capsid\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\huaiyangshan_banyangvirus_NS\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rift_valley_fever_phlebovirus_NS\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rift_valley_fever_phlebovirus_GP\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\crimean_congo_hemorrhagic_fever_orthonairovirus_E\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\california_encephalitis_orthobunyavirus_G2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\california_encephalitis_orthobunyavirus_glycoprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\lassa_mammarenavirus_nucleoprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\lassa_mammarenavirus_Z\\viprbrc.fasta"
#
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rotavirus_a_VP2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rotavirus_a_VP3\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rotavirus_a_NSP1\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\huaiyangshan_banyangvirus_GP\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\huaiyangshan_banyangvirus_GC\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\norovirus_capsid\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\rift_valley_fever_phlebovirus_NS\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\GB_virus_C_polyprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\pegivirus_a_polyprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\zika_virus_polyprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\tick_borne_encephalitis_virus_polyprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_poliovirus_1_VP1\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_poliovirus_1_VP2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_poliovirus_1_VP3\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_poliovirus_1_3D\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_poliovirus_1_2A\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\vaccinia_virus_HA\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\eastern_equine_encephalitis_virus_NSP1\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\eastern_equine_encephalitis_virus_NSP2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\eastern_equine_encephalitis_virus_NSP3\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_herpesvirus_5_UL97\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_herpesvirus_5_UL40\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_herpesvirus_5_UL135\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\human_herpesvirus_5_UL73\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\foot_and_mouth_disease_virus_polyprotein\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\avian_coronavirus_S\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\avian_coronavirus_NSP3\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\avian_coronavirus_NSP4\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\avian_coronavirus_NSP5\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bovine_viral_diarrhea_virus_1_E2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bovine_viral_diarrhea_virus_1_NS5A\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bovine_viral_diarrhea_virus_1_NS5B\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bovine_viral_diarrhea_virus_1_NTPase\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bovine_viral_diarrhea_virus_1_C\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_VP1\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_VP2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_VP3\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_VP4\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_VP5\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_NS1\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_NS2\\viprbrc.fasta"
#fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\bluetongue_virus_NS3\\viprbrc.fasta"
fastafile = "C:\\Google Drive\\GitHub\\protein-evolution\\data\\curated\\H1N1_upto2005_HA\\viprbrc.fasta"







#Binaries.mafft(fastafile)

prefix = joinpath(basename(fastafile), dirname(fastafile))
numsel = 200
k = 3

uniquedict = Dict{String,String}()
names = []
sequences = []
seqnametoindex  = Dict{AbstractString,Int}()
FastaIO.FastaReader(fastafile) do fr
 seqindex = 1
 for (desc, seq) in fr
   len = length(seq)
   if len >= 500 && count(c -> c == 'X', seq) == 0 
     canonicalseq = replace(uppercase(seq), "-" => "")
     h = CommonUtils.sha256base36(canonicalseq)
     if !haskey(uniquedict, h)       
       push!(names,desc)
       push!(sequences, canonicalseq)
       seqnametoindex[desc] = seqindex
       seqindex += 1
       uniquedict[h] = h
     end
   end
 end
end
numsel = min(numsel, length(sequences))
println(length(sequences))
sequences = sequences

numseqs = length(sequences)
distmatrix = zeros(Float64,numseqs,numseqs)

freqvectors = Array{Float64,1}[freqvector(sequence,k) for sequence in sequences]

for i=1:length(freqvectors)
  for j=i+1:length(freqvectors)
    distmatrix[i,j] = euclidean(freqvectors[i],freqvectors[j])
    distmatrix[j,i] = distmatrix[i,j]
  end
end

maxsel = collect(1:numsel)
if numsel != length(sequences)
  sel, maxscore = greedyselection(distmatrix, numsel)
  maxsel = swap(distmatrix, vectobin(sel,numseqs), 2, numsel)
  println(maxsel)
end

alignfile = string(fastafile,".norm", numsel)
outfile = open(alignfile, "w")
for s in maxsel
  #write(outfile, string(">seq", s, "\n"))
  write(outfile, ">$(names[s])\n")
  write(outfile, sequences[s], "\n")
end
close(outfile)

musclealignment, musclealignmentfile = Binaries.muscle(alignfile)
println(abspath(musclealignmentfile))
newick, dummy = Binaries.fasttreeaa(abspath(musclealignmentfile), branchsupport=true)
inputroot = gettreefromnewick(newick)

binarize!(inputroot)
nodelist = getnodelist(inputroot)
for (index,node) in enumerate(nodelist)
  node.nodeindex = index
end
root = gettreefromnewick(getnewick(midpoint_root(nodelist)))
nodelist = getnodelist(root)
for (index,node) in enumerate(nodelist)
  node.nodeindex = index
end

outfile = open(string(prefix,".rooted.nwk"),"w")
write(outfile, getnewick(root))
close(outfile)

names = []
sequences = []
seqnametoindex  = Dict{AbstractString,Int}()
FastaIO.FastaReader(musclealignmentfile) do fr
 seqindex = 1
 for (desc, seq) in fr       
     push!(names,desc)
     push!(sequences, seq)
     seqnametoindex[desc] = seqindex
     seqindex += 1
 end
end

for (nodeindex,node) in enumerate(nodelist)
  node.data = BranchSupportData(getfloatvalue(node.name))
  node.nodeindex = nodeindex
  if haskey(seqnametoindex,node.name)
    node.seqindex = seqnametoindex[node.name]
  end
end

supportmatrix = ones(Float64, numsel, numsel)
for n1 in nodelist
  if n1.seqindex > 0
    for n2 in nodelist
      if n2.seqindex > 0
        anc = getcommonancestor(n1,n2)
        supportmatrix[n1.seqindex,n2.seqindex] = anc.data.branchsupport
        println(n1.name,"\t",n2.name,"\t",anc.name,"\t",anc.data.branchsupport)
      end
    end
  end
end

best = zeros(Int,numsel)
bestcount = 0.0
for iter=1:2000
  global bestcount
  global best
  b1 = zeros(Int,numsel)
  bestbranchlength = 0.0
  for i=1:5000
    oldsupport = evaluatesupport(supportmatrix,b1)
    oldcount = sum(b1)
    oldbranchlength = gettotalbranchlength(nodelist, b1)
    r = rand(1:numsel)
    b1[r] = 1 - b1[r]
    newsupport = evaluatesupport(supportmatrix,b1)
    newcount = sum(b1)
    newbranchlength = gettotalbranchlength(nodelist, b1)
    threshold = 0.99
    if (newsupport > oldsupport && oldsupport < threshold) || ((newcount == 2 || newsupport >= threshold) && newbranchlength > oldbranchlength)
      #println(newsupport,"\t",newbranchlength)
      bestbranchlength = newbranchlength
    else
      b1[r] = 1 - b1[r]
    end
  end
  if bestbranchlength > bestcount
    bestcount = bestbranchlength
    best = copy(b1)
    println(iter,"\t",b1)
  end
end


println(prefix)
alignfile = string(prefix,".select.fasta")
outfile = open(alignfile, "w")
maxsel = bintovec(best)
for s in maxsel
  write(outfile, ">$(names[s])\n")
  write(outfile, sequences[s], "\n")
end
close(outfile)
musclealignment, musclealignmentfile = Binaries.muscle(alignfile)
#newick, dummy = Binaries.fasttreeaa(abspath(musclealignmentfile), branchsupport=true)
Filesystem.cp(musclealignmentfile, string(prefix,".select.fasta"), force=true)

outfile = open(string(prefix,".select.nwk"),"w")
nodelist = getnodelist(removeleaves(nodelist, best))
write(outfile, getnewick(removeleaves(nodelist, String["node1"])))
close(outfile)



exit()
