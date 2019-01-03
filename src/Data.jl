function createalignedprotein(seq1::String, phi1::Array{Float64,1}, psi1::Array{Float64,1}, ss1::String, align1::String)
    seq1align = Int[]
    phi1align = Float64[]
    psi1align = Float64[]
    ss1align = Int[]
    index = 0
    for i=1:length(align1)
        if align1[i] == '#'
            index += 1
            push!(seq1align, indexof(string(seq1[index]), aminoacids))
            push!(phi1align, phi1[index])
            push!(psi1align, psi1[index])
            push!(ss1align, indexof(string(ss1[index]), secondarystructure))
        else
            push!(seq1align, 0)
            push!(phi1align, -1000.0)
            push!(psi1align, -1000.0)
            push!(ss1align, 0)
        end
    end

    protein = BranchState[]
    for (aa,phi,psi) in zip(seq1align,phi1align,psi1align)
        push!(protein, BranchState(0,aa,phi,psi))
    end
    return protein
end

function loadpairs(infile)
    fin = open(infile,"r")

    seq1 = ""
    seq2 = ""
    phi1 = Float64[]
    phi2 = Float64[]
    psi1 = Float64[]
    psi2 = Float64[]
    ss1 = ""
    ss2 = ""
    align1 = ""
    align2 = ""
    lineno = 0
    pairs = Tuple[]
    for line in readlines(fin)
        if length(line) > 0 && line[1] == '>'
            lineno = 0
            seq1 = ""
            seq2 = ""
            phi1 = Float64[]
            phi2 = Float64[]
            psi1 = Float64[]
            psi2 = Float64[]
            ss1 = ""
            ss2 = ""
        elseif lineno == 1
            seq1 = line
        elseif lineno == 2
            seq2 = line
        elseif lineno == 3
            phi1 = Float64[parse(Float64,v) for v in split(line,",")]
        elseif lineno == 4
            psi1 = Float64[parse(Float64,v) for v in split(line,",")]
        elseif lineno == 5
            phi2 = Float64[parse(Float64,v) for v in split(line,",")]
        elseif lineno == 6
            psi2 = Float64[parse(Float64,v) for v in split(line,",")]
        elseif lineno == 7
            ss1 = line
        elseif lineno == 8
            ss2 = line
        elseif lineno == 9
            align1 = line
        elseif lineno == 10
            align2 = line
        end

        if lineno == 10
            protein1 = createalignedprotein(seq1, phi1, psi1, ss1, align1)
            protein2 = createalignedprotein(seq2, phi2, psi2, ss2, align2)
            push!(pairs, (protein1,protein2))
        end

        lineno += 1
    end
    return pairs
end

function createprotein(seq::String, phi_arr::Array{Float64,1}, psi_arr::Array{Float64,1})
    protein = BranchState[]
    for (aa,phi,psi) in zip(seq,phi_arr,psi_arr)
        push!(protein, BranchState(0,aa,phi,psi))
    end
    return protein
end
