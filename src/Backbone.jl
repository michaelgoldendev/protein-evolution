module Backbone
    using LinearAlgebra
    using Formatting
    using DataStructures
    using BioStructures

    push!(LOAD_PATH,@__DIR__)
    using Binaries

    export one_to_three
    one_to_three = DefaultDict{String,String}("XXX")
    one_to_three["A"] = "ALA"
    one_to_three["C"] = "CYS"
    one_to_three["D"] = "ASP"    
    one_to_three["E"] = "GLU"
    one_to_three["F"] = "PHE"
    one_to_three["G"] = "GLY"
    one_to_three["H"] = "HIS"    
    one_to_three["I"] = "ILE"
    one_to_three["L"] = "LEU"    
    one_to_three["K"] = "LYS"  
    one_to_three["M"] = "MET"
    one_to_three["N"] = "ASN"
    one_to_three["P"] = "PRO"
    one_to_three["Q"] = "GLN"
    one_to_three["R"] = "ARG"
    one_to_three["S"] = "SER"    
    one_to_three["T"] = "THR"
    one_to_three["W"] = "TRP"
    one_to_three["V"] = "VAL"
    one_to_three["Y"] = "TYR"

    export three_to_one
    three_to_one = DefaultDict{String,String}("X")
    for oneletter in keys(one_to_three)
        threeletter = one_to_three[oneletter]
        three_to_one[threeletter] = oneletter
    end

    export get_ideal_bond_angle
    function get_ideal_bond_angle(atoms::String, residue::String)
        bond_angle = 0.0
        bond_angle_std = 0.0
        if atoms == "C-N-CA" || atoms == "CA-N-C"
            if residue == "GLY"
                bond_angle = 120.6 / 180 * pi
                bond_angle_std = 1.7  / 180 * pi
            elseif residue == "PRO"
                bond_angle = 122.6 / 180 * pi
                bond_angle_std = 5.0  / 180 * pi
            else
                bond_angle = 121.7 / 180 * pi
                bond_angle_std = 1.8  / 180 * pi
            end
        elseif atoms == "N-CA-C" || atoms == "C-CA-N"
            if residue == "GLY"
                bond_angle = 112.5 / 180 * pi
                bond_angle_std = 2.9  / 180 * pi
            elseif residue == "PRO"
                bond_angle = 111.8 / 180 * pi
                bond_angle_std = 2.5  / 180 * pi
            else
                bond_angle = 111.2 / 180 * pi
                bond_angle_std = 2.8  / 180 * pi
            end
        elseif atoms == "CA-C-N" || atoms == "N-C-CA"
            if residue == "GLY"
                bond_angle = 116.4/ 180 * pi
                bond_angle_std = 2.1  / 180 * pi
            elseif residue == "PRO"
                bond_angle = 116.9 / 180 * pi
                bond_angle_std = 1.5  / 180 * pi
            else
                bond_angle = 116.2 / 180 * pi
                bond_angle_std = 2.0  / 180 * pi
            end
        end
        return bond_angle+pi, bond_angle_std
    end

    export get_ideal_bond_length
    function get_ideal_bond_length(atoms::String, residue::String)
        bond_length = 0.0
        bond_length_std = 0.0
        if atoms == "N-CA" || atoms == "CA-N"
            if residue == "GLY"
                bond_length = 1.451
                bond_length_std = 0.016
            elseif residue == "PRO"
                bond_length = 1.466
                bond_length_std = 0.015
            else
                bond_length = 1.458
                bond_length_std = 0.019
            end
        elseif atoms == "CA-C" || atoms == "C-CA"
            if residue == "GLY"
                bond_length = 1.516
                bond_length_std = 0.018
            else
                bond_length = 1.525
                bond_length_std = 0.021
            end
        elseif atoms == "C-N" || atoms == "N-C"
            if residue == "PRO"
                bond_length = 1.341
                bond_length_std = 0.016
            else
                bond_length = 1.329
                bond_length_std = 0.014
            end
        end
        return bond_length, bond_length_std
    end

    export Atom
    mutable struct Atom
        name::String
        coord::Array{Float64,1}    
        bfactor::Float64
        occupancy::Float64
        element::String
        
        function Atom(name::String, coord::Array{Float64,1}; bfactor::Float64=1.0, occupancy::Float64=1.0, element::String="")
            new(name,coord,bfactor,occupancy,element)
        end
    end

    export Residue
    mutable struct Residue
        resname::String
        atoms::Array{Atom,1}
        
        function Residue(resname::String)
            atoms = Atom[]
            new(resname, atoms)
        end
    end

    export Chain
    mutable struct Chain
        id::String
        residues::Array{Residue,1}
        
        function Chain(id::String)
            residues = Residue[]
            new(id,residues)
        end
    end

    export add_residue
    function add_residue(chain::Chain, residue::Residue)
        push!(chain.residues, residue)
        return residue
    end

    export add_atom
    function add_atom(residue::Residue, atom::Atom)
        push!(residue.atoms, atom)
        return atom
    end

    export writepdb
    function writepdb(io, chain::Chain)
        atom_index = 1
        for (residue_index, residue) in enumerate(chain.residues)
            for atom in residue.atoms
                println(io, "ATOM $(lpad(atom_index,6))  $(rpad(atom.name,4))$(rpad(residue.resname,3)) $(chain.id)$(lpad(residue_index,4))    $(lpad(format(atom.coord[1],precision=3),8))$(lpad(format(atom.coord[2],precision=3),8))$(lpad(format(atom.coord[3],precision=3),8))$(lpad(format(atom.occupancy,precision=2),6))$(lpad(format(atom.bfactor,precision=2),6))  $(lpad(atom.element,10))")
                atom_index += 1
            end
        end
        println(io, "TER  $(lpad(atom_index,6))  $(rpad("",4))$(rpad(chain.residues[end].resname,3)) $(chain.id)$(lpad(length(chain.residues),4))")
        println(io, "END")
    end

    export build_structure_from_angles
    function build_structure_from_angles(polypeptide; use_input_bond_angles::Bool=false, use_input_bond_lengths::Bool=false)
        sequence = polypeptide["sequence"]
        phi_psi = polypeptide["phi_psi"]        
        omega = polypeptide["omega"]
        bond_angles = polypeptide["bond_angles"]
        bond_lengths = polypeptide["bond_lengths"]
        return build_structure_from_angles(sequence, phi_psi, omega, bond_angles, bond_lengths, use_input_bond_angles, use_input_bond_lengths)
    end

    export build_structure_from_angles
    function build_structure_from_angles(sequence, phi_psi, omega, bond_angles, bond_lengths; use_input_bond_angles::Bool=false, use_input_bond_lengths::Bool=false)
        chain = Chain("A")

        residue = add_residue(chain, Residue(one_to_three[string(sequence[1])]))
        N = add_atom(residue, Atom("N", Float64[-4.308, -7.299, 1.075], element="N"))
        CA = add_atom(residue, Atom("CA", Float64[-3.695, -7.053, 2.378], element="C"))
        C = add_atom(residue, Atom("C", Float64[-4.351, -5.868, 3.097], element="C"))

        for pos=2:length(sequence)
            residue = add_residue(chain, Residue(one_to_three[string(sequence[pos])]))

            torsion_angle = phi_psi[pos-1][2] + pi
            bond_angle, bond_angle_std = get_ideal_bond_angle("CA-C-N", residue.resname)
            if use_input_bond_angles                
                bond_angle = bond_angles[pos][1] + pi
            end
            ab = CA.coord - N.coord
            bc = C.coord - CA.coord
            R, R_std = get_ideal_bond_length("C-N", residue.resname)
            if use_input_bond_lengths
                R = bond_lengths[pos][1]
            end
            D = Float64[R*cos(bond_angle), R*cos(torsion_angle)*sin(bond_angle), R*sin(torsion_angle)*sin(bond_angle)]    
            bchat = bc/LinearAlgebra.norm(bc)
            n = cross(ab, bchat)/norm(cross(ab, bchat))
            M = hcat(bchat,cross(n,bchat),n)
            N = add_atom(residue, Atom("N", M*D + C.coord, element="N"))

            torsion_angle = omega[pos] + pi           
            bond_angle, bond_angle_std = get_ideal_bond_angle("C-N-CA", residue.resname)
            if use_input_bond_angles                
                bond_angle = bond_angles[pos][2] + pi
            end
            ab = C.coord - CA.coord
            bc = N.coord - C.coord
            R, R_std = get_ideal_bond_length("N-CA", residue.resname)
            if use_input_bond_lengths
                R = bond_lengths[pos][2]
            end
            D = Float64[R*cos(bond_angle), R*cos(torsion_angle)*sin(bond_angle), R*sin(torsion_angle)*sin(bond_angle)]    
            bchat = bc/LinearAlgebra.norm(bc)
            n = cross(ab, bchat)/norm(cross(ab, bchat))
            M = hcat(bchat,cross(n,bchat),n)
            CA = add_atom(residue, Atom("CA", M*D + N.coord, element="C"))

            torsion_angle = phi_psi[pos][1] + pi            
            bond_angle, bond_angle_std = get_ideal_bond_angle("N-CA-C", residue.resname)
            if use_input_bond_angles                
                bond_angle = bond_angles[pos][3] + pi
            end
            ab = N.coord - C.coord
            bc = CA.coord - N.coord
            R, R_std = get_ideal_bond_length("CA-C", residue.resname)
            if use_input_bond_lengths
                R = bond_lengths[pos][3]
            end
            D = Float64[R*cos(bond_angle), R*cos(torsion_angle)*sin(bond_angle), R*sin(torsion_angle)*sin(bond_angle)]    
            bchat = bc/LinearAlgebra.norm(bc)
            n = cross(ab, bchat)/norm(cross(ab, bchat))
            M = hcat(bchat,cross(n,bchat),n)
            C = add_atom(residue, Atom("C", M*D + CA.coord, element="C"))
        end
        return chain
    end

    export backbone_angles_and_bond_lengths_from_pdb
    function backbone_angles_and_bond_lengths_from_pdb(chain)
        sequence = ""
        phi_psi_list = Tuple{Float64,Float64}[]
        omega_list = Float64[]
        bond_angles_list = Tuple{Float64,Float64,Float64}[]
        bond_lengths_list = Tuple{Float64,Float64,Float64}[]
        Ntempfactor_list = Float64[]        
        CAtempfactor_list = Float64[]
        Ctempfactor_list = Float64[]

        residues = []
        for resid in resids(chain)
            res = chain[resid]
            if !ishetero(res)
                push!(residues, res)
            end
        end

        if length(residues) > 0
            currresnumber = min(1, resnumber(residues[1]))
            pos = 1
            while pos <= length(residues)
                thisresnumber = resnumber(residues[pos])

                aa = "X"
                omega = -1000.0
                phi = -1000.0
                psi = -1000.0
                a1 = -1000.0
                a2 = -1000.0
                a3 = -1000.0
                b1 = -1000.0
                b2 = -1000.0
                b3 = -1000.0

                while currresnumber < thisresnumber && resid(residues[pos]) == string(thisresnumber)
                    sequence = string(sequence, "X")
                    push!(phi_psi_list, (-1000.0,-1000.0))
                    push!(omega_list, -1000.0)
                    push!(bond_angles_list, (-1000.0,-1000.0,-1000.0))
                    push!(bond_lengths_list, (-1000.0,-1000.0,-1000.0))
                    push!(Ntempfactor_list, 0.0)
                    push!(CAtempfactor_list, 0.0)
                    push!(Ctempfactor_list, 0.0)
                    currresnumber += 1
                end

                if resid(residues[pos]) == string(thisresnumber)
                    prevresi = nothing
                    currresi = residues[pos]
                    nextresi = nothing
                    if pos > 1 && resnumber(residues[pos-1])+1 == resnumber(residues[pos])
                        prevresi = residues[pos-1]
                    end
                    if pos < length(residues) && resnumber(residues[pos])+1 == resnumber(residues[pos+1])
                        nextresi = residues[pos+1]
                    end

                    if in("N", atomnames(currresi))
                        push!(Ntempfactor_list, tempfactor(currresi["N"]))
                    else
                        push!(Ntempfactor_list, 0.0)
                    end
                    if in("CA", atomnames(currresi))
                        push!(CAtempfactor_list, tempfactor(currresi["CA"]))
                    else
                        push!(CAtempfactor_list, 0.0)
                    end
                    if in("C", atomnames(currresi))
                        push!(Ctempfactor_list, tempfactor(currresi["C"]))
                    else
                        push!(Ctempfactor_list, 0.0)
                    end

                    if prevresi != nothing
                        if in("CA", atomnames(prevresi)) && in("C", atomnames(prevresi)) && in("N", atomnames(currresi)) && in("CA", atomnames(currresi))
                            omega = dihedralangle(prevresi["CA"], prevresi["C"], currresi["N"], currresi["CA"])
                        end
                        if in("C", atomnames(prevresi)) && in("N", atomnames(currresi)) && in("CA", atomnames(currresi)) && in("C", atomnames(currresi))
                            phi = dihedralangle(prevresi["C"], currresi["N"], currresi["CA"], currresi["C"])
                        end
                        if in("CA", atomnames(prevresi)) && in("C", atomnames(prevresi)) && in("N", atomnames(currresi))
                            a1 = bondangle(prevresi["CA"], prevresi["C"], currresi["N"])
                        end
                        if in("C", atomnames(prevresi)) && in("N", atomnames(currresi)) && in("CA", atomnames(currresi))
                            a2 = bondangle(prevresi["C"], currresi["N"], currresi["CA"])
                        end
                        if in("N", atomnames(currresi)) && in("C", atomnames(prevresi))
                            b1 = sqrt(sqdistance(currresi["N"], prevresi["C"]))
                        end
                    end
                    if in("N", atomnames(currresi)) && in("CA", atomnames(currresi)) && in("C", atomnames(currresi))
                        a3 = bondangle(currresi["N"], currresi["CA"], currresi["C"])
                    end
                    if in("CA", atomnames(currresi)) && in("N", atomnames(currresi))
                        b2 = sqrt(sqdistance(currresi["CA"], currresi["N"]))
                    end
                    if in("C", atomnames(currresi)) && in("CA", atomnames(currresi))
                        b3 = sqrt(sqdistance(currresi["C"], currresi["CA"]))
                    end
                    if nextresi != nothing 
                        if in("N", atomnames(currresi)) && in("CA", atomnames(currresi)) && in("C", atomnames(currresi)) && in("N", atomnames(nextresi))
                            psi  = dihedralangle(currresi["N"], currresi["CA"], currresi["C"], nextresi["N"])
                        end
                    end
                    threeletter = resname(currresi)
                    aa = three_to_one[threeletter]

                    sequence = string(sequence, aa)
                    push!(phi_psi_list, (phi,psi))
                    push!(omega_list, omega)
                    push!(bond_angles_list, (a1,a2,a3))
                    push!(bond_lengths_list, (b1,b2,b3))

                    #println(aa,"\t",currresnumber,"\t", omega, "\t", phi, "\t", psi,"\t",a1,"\t",a2,"\t",a3,"\t",b1,"\t",b2,"\t",b3)V
                    currresnumber += 1                
                end
                pos += 1
            end
        end
        dict = Dict{String,Any}()
        dict["sequence"] = sequence
        dict["phi_psi"] = phi_psi_list
        dict["omega"] = omega_list
        dict["bond_angles"] = bond_angles_list
        dict["bond_lengths"] = bond_lengths_list
        dict["Ntempfactor"] = Ntempfactor_list
        dict["CAtempfactor"] = CAtempfactor_list
        dict["Ctempfactor"] = Ctempfactor_list
        return dict
    end
end
#=
#using Backbone
using BioStructures
#=
using JSON

polypeptide = JSON.parse(open("../data/families/5_3_exonuclease.ali_2.fam", "r"))["proteins"][1]
chain = PDBBuilder.build_structure_from_angles(polypeptide, use_input_bond_angles=false, use_input_bond_lengths=true)
fout = open("ideal-bond-angles.pdb", "w")
PDBBuilder.writepdb(fout, chain)
close(fout)

chain = PDBBuilder.build_structure_from_angles(polypeptide, use_input_bond_angles=true, use_input_bond_lengths=true)
fout = open("actual-bond-angles.pdb", "w")
PDBBuilder.writepdb(fout, chain)
close(fout)=#

pdbdir = abspath("../data/pdbs/")
pdbname = "6n41"
pdbfile = downloadpdb(pdbname, pdb_dir=pdbdir)
structure = read(pdbfile, PDB)
model = structure[1]
chainids = collect(keys(model.chains))
sort!(chainids)
firstchain = model[chainids[1]]

#=
structure = read("../data/pdbs/pdb1a0i.ent", PDB)
model = structure[1]
chainids = collect(keys(model.chains))
sort!(chainids)
firstchain = model[chainids[1]]
println(firstchain)=#
dict = Backbone.backbone_angles_and_bond_lengths_from_pdb(firstchain)
println(dict)=#