include("Main.jl")

function testerrors(parsed_args=Dict{String,Any}())
	rng = MersenneTwister(85649871544631)
	Random.seed!(rand(rng,UInt))
	modelfile = parsed_args["inference_model_file"]

	modelparams = Serialization.deserialize(open(modelfile, "r"))
	family_directories = ["../data/homstrad_curated_highquality/", "../data/curated_rna_virus_structures/"]
	family_names = String[]
	traininghash = zero(UInt)
	trainingexamples = Tuple[]
	for family_dir in family_directories
		family_files = filter(f -> endswith(f,".fam"), readdir(family_dir))
		for family_file in family_files
			full_path = abspath(joinpath(family_dir, family_file))
			json_family = JSON.parse(open(full_path, "r"))
			traininghash = hash(json_family, traininghash)
			if 1 <= length(json_family["proteins"]) <= 1
				training_example = training_example_from_json_family(rng, modelparams, json_family)		
				println(json_family["newick_tree"])
				push!(trainingexamples, training_example)
				push!(family_names, family_file)
			end

			if length(trainingexamples) > 2
				break
			end
		end
	end

	data = Dict{Tuple{Int,Int},Any}()
	for (proteins,nodelist,json_family,sequences) in trainingexamples
		node = nodelist[1]
		backwardsamplesingle(rng, node, modelparams)

		numcols = length(node.data.protein.sites)
		
		for col=1:numcols			
			h = node.data.branchpath.paths[col][end]
			aa = node.data.aabranchpath.paths[col][end]
			bvmphi  = modelparams.hiddennodes[h].phipsi_nodes[aa].mu[1]
			bvmpsi  = modelparams.hiddennodes[h].phipsi_nodes[aa].mu[2]

			aa = 1
			bvmphi  = modelparams.hiddennodes[h].phipsi_node.mu[1]
			bvmpsi  = modelparams.hiddennodes[h].phipsi_node.mu[2]

			h = 1
			key = (h,aa)
			ls = get(data, key, [])
			phipsi = json_family["proteins"][1]["aligned_phi_psi"][col]
			if phipsi[1] > -100.0 && phipsi[2] > -100.0
				CAtempfactor = json_family["proteins"][1]["aligned_CAtempfactor"][col]
				#if haskey(json_family["proteins"][1], "resolution")
					rvalue = json_family["proteins"][1]["resolution"]
					dist = angular_rmsd(bvmphi, phipsi[1], bvmpsi, phipsi[2])
					push!(ls, (dist, CAtempfactor))
					data[key] = ls
				#end
			end
		end
	end

	for key in keys(data)
		ls = data[key]
		for (dist,rvalue) in ls
			println(rvalue,"\t",dist)
		end
		break
	end

end

function parse_errordistributions_commandline()
    settings = ArgParseSettings()
    settings.prog = prog
    settings.version = version
    settings.add_version = true

    add_arg_group(settings, "inference (prediction)")
    @add_arg_table settings begin
        "inference_model_file"
        	help = "specify model file to use for inference"
          	arg_type = String
          	required = true
        "--dataset"
        	help = "alignment in fasta format or .fam file"
        	arg_type = String
        "--tree"
        	help = "tree in newick format"
        	arg_type = String        
        "--blindproteins"
        	help = "comma-seperated list of protein names, whose sequences and structures will be blinded"
        	arg_type = String
    	"--blindstructures"
        	help = "comma-seperated list of protein names, whose structures only  will be blinded"
        	arg_type = String           
        "--samplebranchlengths"
        	help = ""
      	    action = :store_true
        "--samplescalingfactor"
        	help = "keep relative branch lengths fixed, but sample tree scaling factor"
          	action = :store_true
        "--samplesiterates"
        	help = ""
          	action = :store_true
    end


    return parse_args(settings)
end

parsed_args = parse_errordistributions_commandline()
testerrors(parsed_args)