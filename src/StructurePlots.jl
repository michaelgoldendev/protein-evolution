#module StructurePlots
	include("Main.jl")

	push!(LOAD_PATH,@__DIR__)
	#using BivariateVonMises

	using PyPlot
	using AngleUtils
	using StatsBase

	#=
	function cluster()
		rng = MersenneTwister(840184110481014)
		family_directories = ["../data/nonhomologous_singles_xlarge/","../data/homstrad_curated_highquality/", "../data/curated_rna_virus_structures/"]
		#family_directories = []
		family_names = String[]
		traininghash = zero(UInt)
		count = 1
		philist = Float64[]
		psilist = Float64[]
		maxsize = 1000
		k = 30
		for family_dir in family_directories
			family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))
			for family_file in family_files
				full_path = abspath(joinpath(family_dir, family_file))
				json_family = JSON.parse(open(full_path, "r"))
				for z=1:length(json_family["proteins"])
					sequence = json_family["proteins"][z]["sequence"]
					for (pos,(phi,psi)) in enumerate(json_family["proteins"][z]["phi_psi"])
						if phi > -100.0 && psi > -100.0
							if length(philist) < maxsize
								push!(philist,phi)
								push!(psilist,psi)
							else
								r = rand(rng,1:maxsize)
								philist[r] = phi
								psilist[r] = psi
							end
						end
					end
					count += 1
				end
				if count > 500
					break
				end
			end
		end	

		len = length(philist)
		mediods = Int[rand(rng,1:len) for i=1:k]
		assignments = Int[rand(rng,1:k) for r=1:len]
		for iter=1:100
			for m=1:k
				for a=1:len

				end
			end
		end
	end=#

	function index(i::Int, j::Int, dim1::Int, dim2::Int)
		reti = i
		while reti < 1
			reti += dim1
		end
		while reti > dim1
			reti -= dim1
		end

		retj = j
		while retj < 1
			retj += dim1
		end
		while retj > dim2
			retj -= dim2
		end

		return reti,retj
	end

	function createfilter(dim::Int,std::Float64)
		#'dist = Normal(0,std)
		dist = Cauchy(0.0,std)
		filt = zeros(Float64,dim,dim)
		d = div(dim,2)
		for i=1:dim
			for j=1:dim
				a = i - d - 1
				b = j - d - 1
				filt[i,j] = pdf(dist,a)*pdf(dist,b)
			end
		end
		return filt./sum(filt)
	end

	function smoothmat(mat::Array{Float64,2},dim=13,std=2.0)
		filt =  createfilter(dim,std)
		d = div(size(filt,1),2)
		out = zeros(Float64, size(mat,1), size(mat,2))
		for i=1:size(mat,1)
			for j=1:size(mat,2)
				for x=-d:d
					for y=-d:d
						a,b = index(i+x,j+y,size(mat,1),size(mat,2))
						out[i,j] += filt[x+d+1,y+d+1]*mat[a,b]
					end
				end
			end
		end
		return out
	end

	function plot_empirical()
		#json_family["proteins"][seqnametoindex[name]]["aligned_phi_psi"][col]
		N = 500
		fontsize = 24

		conditions = []
		#push!(conditions, ("Glycine","glycine", function(seq,pos) seq[pos] == 'G' end))
		push!(conditions, ("Alanine","alanine", function(seq,pos) seq[pos] == 'A' end))
		push!(conditions, ("Alanine before Proline", "alaninebeforeproline", function(seq,pos) seq[pos] == 'A' && seq[pos+1] == 'P' end))

		matrices = []
		maxval = 0.0
		for (title,name,cond) in conditions
			mat = zeros(Float64, N, N)
			family_directories = ["../data/nonhomologous_singles_xlarge/","../data/homstrad_curated_highquality/", "../data/curated_rna_virus_structures/"]
			#family_directories = []
			family_names = String[]
			traininghash = zero(UInt)
			count = 1
			for family_dir in family_directories
				family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))
				for family_file in family_files
					full_path = abspath(joinpath(family_dir, family_file))
					json_family = JSON.parse(open(full_path, "r"))
					for z=1:length(json_family["proteins"])
						sequence = json_family["proteins"][z]["sequence"]
						for (pos,(phi,psi)) in enumerate(json_family["proteins"][z]["phi_psi"])
							if phi > -100.0 && psi > -100.0 && pos < length(sequence) && cond(sequence,pos)
								x = round(Int, 1 + (mod2pi(phi+pi) / pi / 2.0)*(N-1))
								y = round(Int, 1 + (mod2pi(psi+pi) / pi / 2.0)*(N-1))
								mat[y,x] += 1.0
							end
						end
						count += 1
					end
					if count > 1000
						#break
					end
				end
			end		

			mat = smoothmat(mat)
			mat ./= sum(mat)
			push!(matrices,mat)

			m =  maximum(mat)
			maxval = max(maxval,m)
			println("MAXIMUM ", m)
		end

		for ((title,name,cond),mat) in zip(conditions,matrices)
			mat = mat./maxval
			fig = plt.figure(figsize=(8,8))
			plt.rc("text", usetex=true)
			plt.rc("font", family="serif")
			plt.clf()		

			ax = plt.imshow(mat, vmin=0.0, vmax=1.0)

			angle_tick_positions = [0, div(N-1,4), div(N-1,2), div((N-1)*3,4), N-1]
			angle_labels = ["\$-\\pi\$","\$-\\pi/2\$", "0", "\$\\pi/2\$", "\$\\pi\$"]
			#ax[:set_xticks](angle_tick_positions)
			#ax[:set_yticks](angle_tick_positions)
			plt.xticks(angle_tick_positions, angle_labels, fontsize=fontsize)
			plt.yticks(angle_tick_positions, reverse(angle_labels), fontsize=fontsize)
			plt.title(title, fontsize=28)
			plt.xlabel("Phi (\$\\phi\$)", fontsize=fontsize)
			plt.ylabel("Psi (\$\\psi\$)", fontsize=fontsize, labelpad=-22)
			cb = plt.colorbar(shrink=0.805)
			cb.set_label(label="Normalised density",size=20)
			cb.ax.tick_params(labelsize=20)

			plt.xlim(0.0,N)
			plt.ylim(0.0,N)

			plt.savefig("plots/ramachandran_$(name).png", transparent=true)
			plt.close()
		end
	end

	function plot_nodes(modelfile)

		fin = open(modelfile, "r")
		modelparams = Serialization.deserialize(fin)
		close(fin)

		println("PRIOR ", modelparams.hiddennodes[1].phi_node.kappa_prior)

		toprowonly = true	
		N = 600
		nrows = 2
		figsize = (10,6)
		if toprowonly
			nrows = 1
			figsize = (11,3.6)
			#figsize = (11*1.5,3.6*1.5)
		end
		ncols = 3

		hiddenfreqs = (modelparams.transitionprobs^100)[1,:]

		for h=1:modelparams.numhiddenstates
			#for aa=1:20
				fig = plt[:figure](figsize=figsize)
				plt[:clf]
				plt[:rc]("text", usetex=true)
				plt[:rc]("font", family="serif")

				plt[:suptitle](string(" Hidden state $(h) (", @sprintf("%0.2f", hiddenfreqs[h]*100.0), "\\%)"), fontsize=16,x=0.515)

				mat = zeros(Float64, N, N)
				if modelparams.hidden_conditional_on_aa
					for aa=1:20
						for (i, x) in enumerate(range(-pi,stop=pi,length=N))
							for (j, y) in enumerate(range(-pi,stop=pi,length=N))
								mat[N-j+1,i] += pdf(modelparams.hiddennodes[h].phi_nodes[aa].dist, x)*pdf(modelparams.hiddennodes[h].psi_nodes[aa].dist, y)*modelparams.hiddennodes[h].aa_node.probs[aa]
							end
						end
					end
				else
					for (i, x) in enumerate(range(-pi,stop=pi,length=N))
						for (j, y) in enumerate(range(-pi,stop=pi,length=N))
							mat[N-j+1,i] = pdf(modelparams.hiddennodes[h].phi_node.dist, x)
							mat[N-j+1,i] *= pdf(modelparams.hiddennodes[h].psi_node.dist, y)
						end
					end
				end
				ax = plt[:subplot](nrows, ncols, 1, aspect="auto")
				ax[:imshow](mat)
				angle_tick_positions = [0, div(N-1,4), div(N-1,2), div((N-1)*3,4), N-1]
				angle_labels = ["\$-\\pi\$","\$-\\pi/2\$", "0", "\$\\pi/2\$", "\$\\pi\$"]
				ax[:set_xticks](angle_tick_positions)
				ax[:set_yticks](angle_tick_positions)
				ax[:set_xticklabels](angle_labels)
				ax[:set_yticklabels](reverse(angle_labels))
				plt[:xlabel]("Phi (\$\\phi\$)", fontsize=13)
				plt[:ylabel]("Psi (\$\\psi\$)", fontsize=13)
				
				ax = plt[:subplot](nrows, ncols, 2, aspect="auto")
				x = range(-pi,stop=pi,length=N)
				y = zeros(Float64, N)
				if modelparams.hidden_conditional_on_aa
					for aa=1:20
						y = y .+ pdf.(modelparams.hiddennodes[h].omega_nodes[aa].dist, x)*modelparams.hiddennodes[h].aa_node.probs[aa]
					end
				else
					y = pdf.(modelparams.hiddennodes[h].omega_node.dist, x)
				end
				y = pdf.(modelparams.hiddennodes[h].omega_node.dist, x)
				ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
				ax[:set_xlim](-pi, pi)
				ax[:set_ylim](0.0, 10.0)
				ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
				ax[:set_xticklabels](angle_labels)
				plt[:xlabel]("Omega (\$\\omega\$)", fontsize=13)

				aminoacidstext = ["Ala","Cys","Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr"]
			    barcolors = ["#777775", "#ffd415", "#ef3340", "#ef3340", "#000000", "#ffd415", "#3f48cc", "#333334", "#3f48cc", "#333334", "#333334", "#9e2fe3", "#ffd415", "#9e2fe3", "#3f48cc", "#3f48cc", "#333334", "#777775", "#000000", "#000000"]    
			    barwidth = 0.5
			    ax = plt[:subplot](nrows, ncols, 3, aspect="auto")
				ax[:bar](1:20 + barwidth, modelparams.hiddennodes[h].aa_node.probs, barwidth, color=barcolors, label="Men")
				ax[:set_ylim](0.0, 1.0)
				ax[:set_xticks](1:20 + barwidth + barwidth/2.0)    
			    ax[:set_xticklabels](aminoacidstext, rotation="vertical", fontsize=8)

			    if !toprowonly
				    ax = plt[:subplot](nrows, ncols, 4, aspect="auto")
					x = range(-pi,stop=pi,length=N)
					y = pdf.(modelparams.hiddennodes[h].bond_angle1_node.dist, x)
					ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
					ax[:set_xlim](-pi, pi)
					ax[:set_ylim](0.0, 15.0)
					ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
					ax[:set_xticklabels](angle_labels)
					plt[:xlabel]("\${\\angle}C_{i-1}-N_{i}-CA_{i}\$")

					ax = plt[:subplot](nrows, ncols, 5, aspect="auto")
					x = range(-pi,stop=pi,length=N)
					y = pdf.(modelparams.hiddennodes[h].bond_angle2_node.dist, x)
					ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
					ax[:set_xlim](-pi, pi)
					ax[:set_ylim](0.0, 15.0)
					ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
					ax[:set_xticklabels](angle_labels)
					plt[:xlabel]("\${\\angle}N_{i}-CA_{i}-C_{i}\$")

					ax = plt[:subplot](nrows, ncols, 6, aspect="auto")
					x = range(-pi,stop=pi,length=N)
					y = pdf.(modelparams.hiddennodes[h].bond_angle3_node.dist, x)
					ax[:plot](x, y, color="blue", linewidth=2.0, linestyle="-")
					ax[:set_xlim](-pi, pi)
					ax[:set_ylim](0.0, 15.0)
					ax[:set_xticks]([-pi, -pi/2.0, 0.0, pi/2.0, pi])
					ax[:set_xticklabels](angle_labels)
					plt[:xlabel]("\${\\angle}CA_{i}-C_{i}-N_{i+1}\$")
				end

				plt[:subplots_adjust](bottom=0.2, hspace=0.5)
				#plt[:savefig]("plots/hidden/plot$(h)_nocond.png")
				plt[:savefig]("plots/hidden/plot$(h).png")
				plt[:close]
			#end
		end
	end


	using JSON
	function plotratenetwork(modelfile)
		fin = open(modelfile, "r")
		modelparams = Serialization.deserialize(fin)
		close(fin)

		jsondict = Dict()

		jsondict["numhiddenstates"] = modelparams.numhiddenstates
		jsondict["hiddenfreqs"] = (modelparams.transitionprobs^100)[1,:]

		jsondict["aa_exchangeablities"] = modelparams.aa_exchangeablities
		jsondict["transitionrates"] = modelparams.transitionrates

		for h=1:modelparams.numhiddenstates
			jsondict["aafreqs_h$h"] = modelparams.hiddennodes[h].aa_node.probs
		end

		fout = open("ratenetwork.json", "w")
		JSON.print(fout, jsondict)
		close(fout)
	end

	function angulardist_to_angle(r::Float64)
		return acos((4.0-r*r)/4.0)
	end

	#=
	pdb5kon_A       5.0     0.10497765462660003
	pdb5kon_A       10.0    0.15504705879597758
	pdb5kon_A       15.0    0.1952006839547977
	pdb5kon_A       20.0    0.2363761857881891
	pdb5kon_A       25.0    0.2753585030174491
	pdb5kon_A       30.0    0.31742743975553606
	pdb5kon_A       35.0    0.35733717157303746
	pdb5kon_A       40.0    0.39962262036628277
	pdb5kon_A       45.0    0.4436988090579543
	pdb5kon_A       50.0    0.4930929801444881
	pdb5kon_A       55.0    0.5447908568027151
	pdb5kon_A       60.0    0.6012837939020248
	pdb5kon_A       65.0    0.6681579154888537
	pdb5kon_A       70.0    0.7522127707019409
	pdb5kon_A       75.0    0.8512081083302513
	pdb5kon_A       80.0    0.9853479454377839
	pdb5kon_A       85.0    1.210492764013601
	pdb5kon_A       90.0    1.7017839460994724
	pdb5kon_A       95.0    2.0153764641273937
	pdb5kon_A       100.0   2.827588897883591

	pdb5kon_A       5.0     0.1201021426065796
	pdb5kon_A       10.0    0.17502817016176758
	pdb5kon_A       15.0    0.23282797240820557
	pdb5kon_A       20.0    0.2895681089134446
	pdb5kon_A       25.0    0.35064183163686635
	pdb5kon_A       30.0    0.41814965652203573
	pdb5kon_A       35.0    0.4912150109371237
	pdb5kon_A       40.0    0.5755203605876927
	pdb5kon_A       45.0    0.6703321677459775
	pdb5kon_A       50.0    0.7804351132135359
	pdb5kon_A       55.0    0.9116972833656863
	pdb5kon_A       60.0    1.0829697011592112
	pdb5kon_A       65.0    1.3287474362318772
	pdb5kon_A       70.0    1.6627402575284427
	pdb5kon_A       75.0    1.8936097122000308
	pdb5kon_A       80.0    1.9779699868588898
	pdb5kon_A       85.0    2.0144371434617256
	pdb5kon_A       90.0    2.1132150019257687
	pdb5kon_A       95.0    2.3317382836531944
	pdb5kon_A       100.0   2.8280129048920477
	=#

	function positionofvalue(arr::Array{Float64,1}, v::Float64)
		sort!(arr)
		for (i,a) in enumerate(arr)
			if v < a
				return (i-1.0)/(length(arr)-1.0)
			end
		end
		return 1.0
	end

	function plotradius(plt, colour, r::Float64; linestyle="--", linewidth::Float64=1.0)
		N = 100
		coordinates = []
		xcoordinates = Float64[]
		ycoordinates = Float64[]
		delta = 0.0
		for x in delta:0.001:2.0*pi-delta
			phix = x
			term = (r*r - 4.0 + 2.0*cos(phix-pi))/-2.0
			if abs(term) <= 1.0
				psiy = mod2pi(acos(term))
				if !isnan(phix) && !isnan(psiy)
					if delta <= phix <= 2.0*pi-delta && delta <= psiy <= 2.0*pi-delta
						x1 = phix/2.0/pi*N
						y1 = mod2pi(psiy+pi)/2.0/pi*N
						push!(coordinates, (x1,y1))
						push!(xcoordinates, x1)
						push!(ycoordinates, y1)
					end
				end
			end
		end
		for x in delta:0.001:2.0*pi-delta
			phix = 2.0*pi - x
			term = (r*r - 4.0 + 2.0*cos(phix-pi))/-2.0
			if abs(term) <= 1.0
				psiy = mod2pi(acos(term))
				if !isnan(phix) && !isnan(psiy)
					if delta <= phix <= 2.0*pi-delta && delta <= psiy <= 2.0*pi-delta
						x1 = phix/2.0/pi*N
						y1 = mod2pi(-psiy-pi)/2.0/pi*N
						push!(coordinates, (x1,y1))
						push!(xcoordinates, x1)
						push!(ycoordinates, y1)
					end
				end
			end
		end
		push!(xcoordinates,xcoordinates[1])
		push!(ycoordinates,ycoordinates[1])
		#println(xcoordinates)
		#println(ycoordinates)
		plt.plot(xcoordinates, ycoordinates, linestyle, color=colour, linewidth=linewidth)
	end

	function plotaccuracy(samplefile="output/2019-05-08.14h02m51s.408.2lae8pznqv3xj/output.samples")
		outdir = dirname(samplefile)
		N = 100
		startpos = 1
		fontsize = 24

		fin = open(samplefile, "r")	
		samples = Serialization.deserialize(fin)
		close(fin)
		for name in keys(samples)
			if !startswith(name,"metadata:")			
				proteinsample = samples[name]
				json_family = proteinsample.json_family
				seqnametoindex = Dict{String,Int}()
				for (index,protein) in enumerate(proteinsample.json_family["proteins"])
					seqnametoindex[proteinsample.json_family["proteins"][index]["name"]] = index
				end
				if length(filter(x -> x[1] > -100.0 || x[2] > -100.0, proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_phi_psi"])) > 0
					modelparams = proteinsample.modelparams
					hiddenfreqs = (modelparams.transitionprobs^100)[1,:]

					numcols = length(proteinsample.aasamples[1])
					fig = plt.figure(figsize=(8,8))
					plt.rc("text", usetex=true)
					plt.rc("font", family="serif")

					numsamples = length(proteinsample.aasamples)
					startiter = max(1, div(numsamples,3))
					enditer = numsamples

					distances = Float64[]
					for iter=startiter:numsamples
						for col=1:numcols
							truephipsi = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_phi_psi"][col]
							sampledphipsi = proteinsample.phipsisamples[iter][col]
							if truephipsi[1] > -100.0 && truephipsi[2] > -100.0 && sampledphipsi[1] > -100.0 && sampledphipsi[2] > -100.0
								push!(distances, angular_rmsd(truephipsi[1],sampledphipsi[1],truephipsi[2],sampledphipsi[2]))
							end
						end			
					end
					meandistance = median(distances)
					#println("Median: $(meandistance)")

					cpos = [0.0, 0.75, 1.25, 1.75]
					rs = Float64[0.5, 1.0, 1.5, 2.0]
					rcolours = ["wo", "wo", "wo", "wo"]

					cpos = [0.0, 0.380-0.025, 0.75-0.025, 1.25-0.025, 1.75-0.025]
					rs = Float64[0.25, 0.5, 1.0, 1.5, 2.0]
					distlabelr = rs .+ 0.065
					distlabelr[1] = 0.17
					percfontsizes = Int[21,22,24,24,24]
					distfontsizes = Int[20,24,24,24,24]
					#colours = ["#65ca44", "#84d432", "#abdb20", "#d4e115", "#fce51e"]
					colours = reverse(["#5599ff", "#2a7fff", "#0066ff", "#0055d4", "#0044aa", "#002255"])


					#=
					if length(distances) > 0
						for perc=5.0:5.0:100.0
							println(name,"\t",perc,"\t",StatsBase.percentile(distances,perc))
						end
						for r in rs
							println(r,"\t",positionofvalue(distances,r))
						end
					end=#
					mat = zeros(Float64, N, N)

					plt.clf()
					
					#ax = plt.imshow(mat)

					ax = plt.gca()
					ax.set_facecolor("white")
					ax.set_aspect("equal", "box")

					#perc =	@sprintf("%d", positionofvalue(distances,0.25)*100.0)
					#plt.text(0.5*N,0.5*N, string("$(perc)\\%"), fontsize=percfontsize, color="white",  horizontalalignment="center", verticalalignment="center")
					#plt.text((0.5+(0.15/pi))*N,0.5*N, string("0.25"), fontsize=10, color="white",  horizontalalignment="left", verticalalignment="center")

					#plotradius(plt, "#999999", meandistance; linestyle="--")

					for (colour,percfontsize,distfontsize,c,d,r) in zip(colours,percfontsizes,distfontsizes,cpos,distlabelr,rs)
						perc =	@sprintf("%d", positionofvalue(distances,r)*100.0)
						theta = acos(1.0 - c*c/2.0)/2.0
						plt.text(0.5*N,(0.5+(theta/pi))*N, string("$(perc){\\normalsize{\\%}}"), fontsize=percfontsize, color=colour,  horizontalalignment="center", verticalalignment="center")
						theta = angulardist_to_angle(d)/2.0
						plt.text((0.5+(theta/pi))*N,(0.5+(theta/pi))*N, string(r), fontsize=distfontsize, color=colour,  horizontalalignment="left", verticalalignment="center")
					end

					for (colour,r) in zip(colours,rs)
						plotradius(plt, colour, r)
					end

					#=
					for (r, pointcol) in zip(rs,rcolours)
						for x in 0.0:0.0025:2.0*pi
							phix = x
							term = (r*r - 4.0 + 2.0*cos(phix))/-2.0
							if abs(term) <= 1.0
								psiy = mod2pi(acos(term))
								if !isnan(phix) && !isnan(psiy)
									if 0.0 <= phix <= 2.0*pi && 0.0 <= psiy <= 2.0*pi
										plt.plot(mod2pi(phix-pi)/2.0/pi*N, mod2pi(psiy-pi)/2.0/pi*N, pointcol, markersize=0.25)
									end
								end
								if !isnan(phix) && !isnan(psiy)
									if 0.0 <= phix <= 2.0*pi && 0.0 <= psiy <= 2.0*pi
										plt.plot(mod2pi(phix-pi)/2.0/pi*N, mod2pi(-psiy-pi)/2.0/pi*N, pointcol, markersize=0.25)
									end
								end
							end
						end
					end=#

				

					angle_tick_positions = [0, div(N-1,4), div(N-1,2), div((N-1)*3,4), N-1]
					angle_labels = ["\$-\\pi\$","\$-\\pi/2\$", "0", "\$\\pi/2\$", "\$\\pi\$"]
					#ax[:set_xticks](angle_tick_positions)
					#ax[:set_yticks](angle_tick_positions)
					plt.xticks(angle_tick_positions, angle_labels, fontsize=fontsize)
					plt.yticks(angle_tick_positions, reverse(angle_labels), fontsize=fontsize)
					#plt.title("Site $(col)", fontsize=15)
					plt.xlabel("Phi (\$\\phi\$)", fontsize=fontsize)
					plt.ylabel("Psi (\$\\psi\$)", fontsize=fontsize, labelpad=-22)

					plt.xlim(0.0,N)
					plt.ylim(0.0,N)

					outplotdir = joinpath(outdir, "plots")
					if !isdir(outplotdir)
						mkdir(outplotdir)	
					end
					plt.savefig(joinpath(outplotdir, "target_$(name)"), transparent=false, dpi=600)
					plt.close()
				end
			end
		end
	end

	function plotstructuresamples(samplefile="output/2019-05-08.14h35m19s.66.2lae8pznqv3xj/output.samples", othernames::Array{AbstractString,1}=["pdb4rpd_A"])
		outdir = dirname(samplefile)
		N = 100
		startpos = 1

		fin = open(samplefile, "r")	
		samples = Serialization.deserialize(fin)
		close(fin)
		for name in keys(samples)
			if !startswith(name,"metadata:")
				#=
				if name != plotname
					continue
				end=#
				proteinsample = samples[name]

				json_family = proteinsample.json_family

				seqnametoindex = Dict{String,Int}()
				for (index,protein) in enumerate(proteinsample.json_family["proteins"])
					seqnametoindex[proteinsample.json_family["proteins"][index]["name"]] = index
				end

				modelparams = proteinsample.modelparams
				hiddenfreqs = (modelparams.transitionprobs^100)[1,:]

				numcols = length(proteinsample.aasamples[1])
				fig = plt.figure(figsize=(7,5))
				plt.rc("text", usetex=true)
				plt.rc("font", family="serif")

				numsamples = length(proteinsample.aasamples)
				startiter = max(1, div(numsamples,3))
				enditer = numsamples

				distances = Float64[]
				for iter=startiter:numsamples
					for col=1:numcols
						truephipsi = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_phi_psi"][col]
						sampledphipsi = proteinsample.phipsisamples[iter][col]
						if truephipsi[1] > -100.0 && truephipsi[2] > -100.0 && sampledphipsi[1] > -100.0 && sampledphipsi[2] > -100.0
							push!(distances, angular_rmsd(truephipsi[1],sampledphipsi[1],truephipsi[2],sampledphipsi[2]))
						end
					end			
				end

				for col=startpos:numcols
					phi,psi = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_phi_psi"][col]
					aa = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_sequence"][col]

					if phi > -100.0 && psi > -100.0
						mat = zeros(Float64, N, N)				
						for iter=startiter:numsamples
							h = proteinsample.hiddensamples[iter][col]
							aa = proteinsample.aasamples[iter][col]

							#k = Float64[modelparams.hiddennodes[h].phi_nodes[aa].kappa, modelparams.hiddennodes[h].psi_nodes[aa].kappa, 0.0]
							#mu = Float64[modelparams.hiddennodes[h].phi_nodes[aa].mu, modelparams.hiddennodes[h].psi_nodes[aa].mu]

							for (i, x) in enumerate(range(-pi,stop=pi,length=N))
								for (j, y) in enumerate(range(-pi,stop=pi,length=N))
									#mat[N-j+1,i] += pdf(modelparams.hiddennodes[h].phi_nodes[aa].dist, x)*pdf(modelparams.hiddennodes[h].psi_nodes[aa].dist, y)*modelparams.hiddennodes[h].aa_node.probs[aa]
									if modelparams.hidden_conditional_on_aa
										#mat[N-j+1,i] += BivariateVonMises.pdf(modelparams.hiddennodes[h].phipsi_nodes[aa], Float64[x,y])
										mat[j,i] += BivariateVonMises.pdf(modelparams.hiddennodes[h].phipsi_nodes[aa], Float64[x,y])
									else
										#mat[N-j+1,i] += BivariateVonMises.pdf(modelparams.hiddennodes[h].phipsi_node, Float64[x,y])
										mat[j,i] += BivariateVonMises.pdf(modelparams.hiddennodes[h].phipsi_node, Float64[x,y])
									end
								end
							end
						end
						mat /= (enditer-startiter+1.0)

						plt.clf()
						#plt[:suptitle](string(" Hidden state $(h) (", @sprintf("%0.1f", hiddenfreqs[h]*100.0), "\\%)"), fontsize=16,x=0.515)
						#ax = plt[:gca]

						
						ax = plt.imshow(mat)
						phi,psi = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_phi_psi"][col]
						aa = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_sequence"][col]
						x = ((phi+pi)/(2.0*pi))*N
						#y = (1.0-((psi+pi)/(2.0*pi)))*N
						y = (((psi+pi)/(2.0*pi)))*N
						#plt.plot(x, y, "o", "None", markersize=12, markeredgecolor="k")
						plt.text(x+1*(100.0/N),y+1*(100.0/N), string("\\textsf{\\textbf{",aa,"}}"), color="white",  horizontalalignment="center", verticalalignment="center")

						for othername in othernames
							if haskey(seqnametoindex, othername)
								phihomolog,psihomolog = proteinsample.json_family["proteins"][seqnametoindex[othername]]["aligned_phi_psi"][col]
								aa = proteinsample.json_family["proteins"][seqnametoindex[othername]]["aligned_sequence"][col]
								if phihomolog > -100.0 && psihomolog > -100.0
									x = ((phihomolog+pi)/(2.0*pi))*N
									#y = (1.0-((psihomolog+pi)/(2.0*pi)))*N
									y = (((psihomolog+pi)/(2.0*pi)))*N
									#plt.plot(x, y, "o", "None", markersize=10, markeredgecolor="r")
									plt.text(x+1*(100.0/N),y+1*(100.0/N), string("\\textsf{\\textbf{",aa,"}}"), color="red",  horizontalalignment="center", verticalalignment="center")
								end
							end
						end

					

						angle_tick_positions = [0, div(N-1,4), div(N-1,2), div((N-1)*3,4), N-1]
						angle_labels = ["\$-\\pi\$","\$-\\pi/2\$", "0", "\$\\pi/2\$", "\$\\pi\$"]
						#ax[:set_xticks](angle_tick_positions)
						#ax[:set_yticks](angle_tick_positions)
						plt.xticks(angle_tick_positions, angle_labels)
						plt.yticks(angle_tick_positions, reverse(angle_labels))
						plt.title("Site $(col)", fontsize=15)
						plt.xlabel("Phi (\$\\phi\$)", fontsize=13)
						plt.ylabel("Psi (\$\\psi\$)", fontsize=13)

						plt.xlim(0.0,N)
						plt.ylim(0.0,N)

						outplotdir = joinpath(outdir, "samples")
						if !isdir(outplotdir)
							mkdir(outplotdir)	
						end
						plt.savefig(joinpath(outplotdir, "$(name)_site$(col).png"), transparent=true)
						plt.close()
					end
				end
			end
		end
	end

	function parse_plotting_commandline()
	    settings = ArgParseSettings()
	    settings.prog = prog
	    settings.version = version
	    settings.add_version = true

	    add_arg_group(settings, "plotting")
	    @add_arg_table settings begin
	        "plot_model_file"
	        	help = "specify model file to plot hidden nodes"
	          	arg_type = String
	          	required = true
	    end

	    return parse_args(settings)
	end

	#parsed_args = parse_plotting_commandline()
	#plot_nodes(parsed_args["plot_model_file"])
	#plotratenetwork(parsed_args["plot_model_file"])
	#plotratenetwork("models/model_h.1.thresh3.rerun.hiddenaascaling.ratemode1.model")
	#plot_nodes("models/model_h.20.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.16.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.15.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.21.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.14.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.17.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.10.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.11.thresh2.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.10.thresh2.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.12.thresh2.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.15.thresh2.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.20.thresh3.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.10.thresh3.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
	#plot_nodes("models/model_h.30.thresh2.rerun.hiddenaascaling.nobranchsampling.ratemode1.model")
	#plot_nodes("models/model_h.30.thresh2.rerun.hiddenaascaling.nobranchsampling.anglescondaa15.ratemode1.model")
	#plotratenetwork("models/model_h.10.thresh3.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
	#plotaccuracy()
	#plotstructuresamples()
	#plot_empirical()
	#cluster()
#end

#using StructurePlots

plot_nodes("models/model.h40.2glmjug4seflu.thresh2.hiddenaascaling.anglescondaa25.precluster.ratemode1.model")
