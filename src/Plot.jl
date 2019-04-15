include("Main.jl")

push!(LOAD_PATH,@__DIR__)
#using BivariateVonMises

using PyPlot
using AngleUtils
using StatsBase

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
		    barcolors = ["#777775", "#fedd00", "#ef3340", "#ef3340", "#000000", "#fedd00", "#0087c7", "#333334", "#0087c7", "#333334", "#333334", "#65428a", "#fedd00", "#65428a", "#0087c7", "#0087c7", "#333334", "#777775", "#000000", "#000000"]    
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
			plt[:savefig]("plots/hidden/plot$(h).svg")
			#plt[:savefig]("plots/hidden/plot$(h)_$(aminoacids[aa]).png")
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

function plotstructuresamples()
	othername = "pdb4rpd_A"
	#othername = ""
	samplefile = string("structure.samples")

	N = 100
	startpos = 1

	fin = open(samplefile, "r")	
	samples = Serialization.deserialize(fin)
	close(fin)
	for name in keys(samples)
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
		if length(distances) > 0
			for perc=5.0:5.0:100.0
				println(name,"\t",perc,"\t",StatsBase.percentile(distances,perc))
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
								mat[N-j+1,i] += BivariateVonMises.pdf(modelparams.hiddennodes[h].phipsi_nodes[aa], Float64[x,y])
							else
								mat[N-j+1,i] += BivariateVonMises.pdf(modelparams.hiddennodes[h].phipsi_node, Float64[x,y])
							end
						end
					end
				end
				mat /= (enditer-startiter+1.0)

				plt.clf()
				#plt[:suptitle](string(" Hidden state $(h) (", @sprintf("%0.1f", hiddenfreqs[h]*100.0), "\\%)"), fontsize=16,x=0.515)
				#ax = plt[:gca]

				
				ax = plt.imshow(mat)

				#=
				#plt[:plot](0.5*N, 0.5*N, "wo", markersize=1.0*N)
				#plt[:scatter](0.5*N, 0.5*N, s=2.0*(N)^2.0, facecolors="w", edgecolors="r")
				#plt[:gcf]()[:gca]()[:add_artist](plt[:Circle]((0.5*N, 0.5*N), 0.5*angulardist_to_angle(1.25)/pi*N, fill=false, edgecolor="w"))
				#plt[:gcf]()[:gca]()[:add_artist](plt[:Circle]((0.5*N, 0.5*N), 0.5*angulardist_to_angle(1.5)/pi*N, fill=false, edgecolor="w"))
				#plt[:gcf]()[:gca]()[:add_artist](plt[:Circle]((0.5*N, 0.5*N), 0.5*angulardist_to_angle(sqrt(8)-0.001)/pi*N, fill=false, edgecolor="w"))
				#rs = Float64[0.911, 1.465, 2.03]
				#rs = Float64[0.394, 0.797, 1.859, 2.106]
				#rs = Float64[0.327, 0.648, 1.493, 2.031] # far
				#rs = Float64[0.237, 0.424, 0.733, 1.366] # close
				#rs = Float64[0.663, 1.806, 2.06, 2.317] # seqonly
				#rs = Float64[0.679, 1.707, 2.042, 2.284] # random
				#rs = Float64[0.294, 0.556, 1.129, 1.993] # far, rate variation
				rs = Float64[0.232, 0.416, 0.72, 1.269] # close, rate variation
				rs = Float64[0.318, 0.555, 0.913, 1.578] # far, rate variation, no angle cond
				rcolours = ["wo", "wo", "wo", "wo"]
				for x in 0.0:0.005:2.0*pi
					for (r, pointcol) in zip(rs,rcolours)
						phix = x
						term = (r*r - 4.0 + 2.0*cos(phix))/-2.0
						if abs(term) <= 1.0
							psiy = mod2pi(acos(term))
							if !isnan(phix) && !isnan(psiy)
								if 0.0 <= phix <= 2.0*pi && 0.0 <= psiy <= 2.0*pi
									plt[:plot](mod2pi(phix-pi)/2.0/pi*N, mod2pi(psiy-pi)/2.0/pi*N, pointcol, markersize=0.25)
								end
							end
							if !isnan(phix) && !isnan(psiy)
								if 0.0 <= phix <= 2.0*pi && 0.0 <= psiy <= 2.0*pi
									plt[:plot](mod2pi(phix-pi)/2.0/pi*N, mod2pi(-psiy-pi)/2.0/pi*N, pointcol, markersize=0.25)
								end
							end
						end
					end
				end=#
				phi,psi = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_phi_psi"][col]
				aa = proteinsample.json_family["proteins"][seqnametoindex[name]]["aligned_sequence"][col]
				x = ((phi+pi)/(2.0*pi))*N
				y = (1.0-((psi+pi)/(2.0*pi)))*N
				#plt.plot(x, y, "o", "None", markersize=12, markeredgecolor="k")
				plt.text(x+1*(100.0/N),y+1*(100.0/N), string("\\textsf{\\textbf{",aa,"}}"), color="white",  horizontalalignment="center", verticalalignment="center")

				if haskey(seqnametoindex, othername)
					phihomolog,psihomolog = proteinsample.json_family["proteins"][seqnametoindex[othername]]["aligned_phi_psi"][col]
					aa = proteinsample.json_family["proteins"][seqnametoindex[othername]]["aligned_sequence"][col]
					if phihomolog > -100.0 && psihomolog > -100.0
						x = ((phihomolog+pi)/(2.0*pi))*N
						y = (1.0-((psihomolog+pi)/(2.0*pi)))*N
						#plt.plot(x, y, "o", "None", markersize=10, markeredgecolor="r")
						plt.text(x+1*(100.0/N),y+1*(100.0/N), string("\\textsf{\\textbf{",aa,"}}"), color="red",  horizontalalignment="center", verticalalignment="center")
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

				plt.savefig("plots/$(name)_site$(col).png", transparent=true)
				plt.close()
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
#plot_nodes("models/model_h.30.thresh2.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
#plotratenetwork("models/model_h.10.thresh3.rerun.hiddenaascaling.anglescondaa.ratemode1.model")
plotstructuresamples()