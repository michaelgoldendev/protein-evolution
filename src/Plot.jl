include("Main.jl")


using PyPlot
function plot_nodes(modelfile)

	fin = open(modelfile, "r")
	modelparams = Serialization.deserialize(fin)
	close(fin)

	toprowonly = true	
	N = 300
	nrows = 2
	figsize = (10,6)
	if toprowonly
		nrows = 1
		figsize = (11,3.6)
	end
	ncols = 3

	hiddenfreqs = (modelparams.transitionprobs^100)[1,:]

	for h=1:modelparams.numhiddenstates
		fig = plt[:figure](figsize=figsize)
		plt[:clf]
		plt[:rc]("text", usetex=true)
		plt[:rc]("font", family="serif")

		plt[:suptitle](string(" Hidden state $(h) (", @sprintf("%0.1f", hiddenfreqs[h]*100.0), "\\%)"), fontsize=16,x=0.515)

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
		plt[:savefig]("plot$(h).png")
		plt[:close]
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

parsed_args = parse_plotting_commandline()
#plot_nodes(parsed_args["plot_model_file"])
plotratenetwork(parsed_args["plot_model_file"])