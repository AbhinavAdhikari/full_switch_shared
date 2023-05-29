
##
include("Include.jl")
function objective_function(parameter_guess_array,time_start,time_step_size,time_stop,model_data_dictionary, exp_data_dictionary)  

    # what is the host_type?
    host_type = :cell_free

    # Phase 1: parameter update =========================================================================== #
    # update the paramaters in the model data dictionary -
    # for now - lets only search over dG's -
    R = model_data_dictionary["R"]
    T_K = model_data_dictionary["T_K"] 

    # compute W -
    tmp_W_array = Float64[]
    for index = 1:8
        parameter_guess = parameter_guess_array[index]
        value = exp(-1*(parameter_guess))
        push!(tmp_W_array,value)
    end

    # update the control W's -
    control_parameter_dictionary = model_data_dictionary["control_parameter_dictionary"]
	control_parameter_dictionary["W_RNAP_P70"] = tmp_W_array[1]
    control_parameter_dictionary["W_RNAP_P28"] = tmp_W_array[2]
    control_parameter_dictionary["W_RNAP_mP70"] = tmp_W_array[3]
    control_parameter_dictionary["W_S70_RNAP_P70"] = tmp_W_array[4]
	control_parameter_dictionary["W_S70_RNAP_mP70"] = tmp_W_array[5]
	control_parameter_dictionary["W_S28_RNAP_P28"] = tmp_W_array[6]
	control_parameter_dictionary["W_GntR_mP70"] = tmp_W_array[7]
	control_parameter_dictionary["W_AS28_S28"] = tmp_W_array[8] # for now just keep this but tighten the bounds.
    model_data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	#print(control_parameter_dictionary)

    binding_parameter_dictionary = model_data_dictionary["binding_parameter_dictionary"]
    binding_parameter_dictionary["n_S70_RNAP_P70"]=parameter_guess_array[9]
	binding_parameter_dictionary["K_S70_RNAP_P70"]=parameter_guess_array[10]
	binding_parameter_dictionary["n_S70_RNAP_mP70"]=parameter_guess_array[11]
	binding_parameter_dictionary["K_S70_RNAP_mP70"]=parameter_guess_array[12]
	binding_parameter_dictionary["n_GntR_mP70"]=parameter_guess_array[13]
    binding_parameter_dictionary["K_GntR_mP70"]=parameter_guess_array[14]
	binding_parameter_dictionary["n_S28_RNAP_P28"]=parameter_guess_array[15]
    binding_parameter_dictionary["K_S28_RNAP_P28"]=parameter_guess_array[16]
	binding_parameter_dictionary["n_AS28_S28"]=parameter_guess_array[17]
	binding_parameter_dictionary["K_AS28_S28"]=parameter_guess_array[18]
    model_data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

	#print(binding_parameter_dictionary)
    # time constant modifier -

	time_constant_modifier_array = [
		0.0							;	# 1	GntR
		0.0							;	# 2	S28
		0.0							;	# 3	AS28
        0.0                         ;   # 4 Venus 
        0.0                         ;   # 5 BFP
		parameter_guess_array[19]	;	# 6	mRNA_GntR 
		parameter_guess_array[20]	;	# 7	mRNA_S28
		parameter_guess_array[21]	;	# 8	mRNA_AS28
        parameter_guess_array[22]	;	# 9	mRNA_Venus							;
        parameter_guess_array[23]   ;   # 10 mRNA_BFP
		parameter_guess_array[24]	;	# 11	protein_GntR
		parameter_guess_array[25]	;	# 12	protein_S28
        parameter_guess_array[26]	;	# 13	protein_AS28
        parameter_guess_array[27]	;	# 14	protein_Venus
        parameter_guess_array[28]	;	# 15	protein_BFP
	]
    model_data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    # setup degradation_modifier_array -

	degradation_modifier_array = [
		0.0	;	# 1	GntR
		0.0	;	# 2	Venus
		0.0	;	# 3	sigma_70
        0.0 ;   # 4 Venus
        0.0 ;   # 5 BFP
		parameter_guess_array[29]	;	# 6	mRNA_GntR
		parameter_guess_array[30]	;	# 7	mRNA_S28
		parameter_guess_array[31]   ;   # 8 mRNA_AS28
        parameter_guess_array[32]   ;   # 9 mRNA_Venus
        parameter_guess_array[33]   ;   # 10 mRNA_BFP
		parameter_guess_array[34]	;	# 11	protein_GntR
		parameter_guess_array[35]	;	# 12	protein_Venus
		parameter_guess_array[36]	;	# 13	protein_sigma_70
        parameter_guess_array[37]	;	# 14	protein_Venus
        parameter_guess_array[38]	;	# 13	protein_BFP
	]
	model_data_dictionary["degradation_modifier_array"] = degradation_modifier_array

    # update the translation time -
    model_data_dictionary["half_life_translation_capacity"] = parameter_guess_array[39]

    # lastly, update KL and KX -
    biophysical_constants_dictionary = model_data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["translation_saturation_constant"] = parameter_guess_array[40]
    biophysical_constants_dictionary["transcription_saturation_constant"] = parameter_guess_array[41]

    model_data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

	# gluconate GntR binding parameters
	gluconate_parameter_dictionary = model_data_dictionary["gluconate_parameter_dictionary"]
	gluconate_parameter_dictionary["n_gluconate_GntR"] = parameter_guess_array[42]
	gluconate_parameter_dictionary["K_gluconate_GntR"] = parameter_guess_array[43]
	model_data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary

    # sigma 70 half_life
    model_data_dictionary["half_life_sigma_70"] = parameter_guess_array[44]

    # grab defaults -
    species_symbol_type_array = model_data_dictionary["species_symbol_type_array"]
    protein_coding_length_array = model_data_dictionary["protein_coding_length_array"]
    gene_coding_length_array = model_data_dictionary["gene_coding_length_array"]
    time_constant_modifier_array = model_data_dictionary["time_constant_modifier_array"]
    initial_condition_array = model_data_dictionary["initial_condition_array"]

    # # get gene IC -
    idx_gene = findall(x->x==:gene,species_symbol_type_array)
    gene_abundance_array = initial_condition_array[idx_gene]

    # Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
    model_data_dictionary["translation_parameter_array"] = translation_parameter_array

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
    model_data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

    # Dilution degrdation matrix -
    dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
    model_data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
    # ===================================================================================================== #
	#print(model_data_dictionary)
    # Phase 2:  solve model equations ===================================================================== #
    # solve the balance equations -
     
    (TSIM,XSIM)= SolveBalances(time_start, time_stop, time_step_size, model_data_dictionary)

    # ===================================================================================================== #
	#print(TSIM, XSIM)

    # Phase 3:  compute simulation error ================================================================== #
    # # compute the error - we need to do a bunch of interpolation -

	 tsim_exp_protein = exp_data_dictionary["prot_data_array"][:,1]
	 tsim_exp_mRNA = exp_data_dictionary["mRNA_data_array"][:,1]


    #Venus Protein

    itp_Venus_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,14]);
    protein_Venus_sim = itp_Venus_protein[tsim_exp_protein]

    # # BFP Protein

    itp_BFP_protein =  Interpolations.LinearInterpolation(TSIM, XSIM[:,15]);
    protein_BFP_sim = itp_BFP_protein[tsim_exp_protein]


    #only for BFP and Venus do we have experimental data and mRNA for venus and GntR

    protein_Venus_exp = exp_data_dictionary["prot_data_array"][:,2]       # 
	protein_Venus_std_exp = exp_data_dictionary["prot_data_array"][:,3]   # 
 
    protein_BFP_exp = exp_data_dictionary["prot_data_array"][:,5]       # mean 
	protein_BFP_std_exp = exp_data_dictionary["prot_data_array"][:,6]   # stdev 

    #mRNA data

    itp_Venus_mRNA =  Interpolations.LinearInterpolation(TSIM, (1000)*XSIM[:,9]);
    # mRNA_Venus_sim = itp_Venus_mRNA[tsim_exp_protein]
    mRNA_Venus_sim = itp_Venus_mRNA[tsim_exp_mRNA]

    itp_BFP_mRNA =  Interpolations.LinearInterpolation(TSIM, (1000)*XSIM[:,10]);
    # mRNA_BFP_sim = itp_BFP_mRNA[tsim_exp_protein]
    mRNA_BFP_sim = itp_BFP_mRNA[tsim_exp_mRNA]



    mRNA_Venus_exp = exp_data_dictionary["mRNA_data_array"][:,2]  #nM
    mRNA_Venus_exp_stderr=exp_data_dictionary["mRNA_data_array"][:,3]  
    mRNA_BFP_exp= exp_data_dictionary["mRNA_data_array"][:,4] #nM
    mRNA_BFP_exp_stderr=exp_data_dictionary["mRNA_data_array"][:,5]    


    itp_Venus_mRNA =  Interpolations.LinearInterpolation(TSIM, 1000*XSIM[:,9]);
    B = DataInterpolations.LinearInterpolation(mRNA_Venus_exp,tsim_exp_mRNA)
    # mRNA_Venus_exp = B.(tsim_exp_protein)
    mRNA_Venus_exp = B.(tsim_exp_mRNA)


    itp_BFP_mRNA =  Interpolations.LinearInterpolation(TSIM,1000*XSIM[:,10]); #uM to nM
    C = DataInterpolations.LinearInterpolation(mRNA_BFP_exp,tsim_exp_mRNA)
    # mRNA_BFP_exp = C.(tsim_exp_protein)
    mRNA_BFP_exp = C.(tsim_exp_mRNA)

    #  compute error terms -
    # error_term_array = zeros(4,1) 
    error_term_array = zeros(4,1) # 2 mRNA and protein each. 3 for penalizing GntR, S28 and AS28. 

    # # mRNA Venus -
	
    error_vector_1 = mRNA_Venus_exp .- mRNA_Venus_sim

    # protein Venus -
    # protein_Venus_std_exp[1] = 1.0
    # tmp_arr1 = 1.0./((protein_Venus_std_exp).^2)
    # W_prot1 = diagm(tmp_arr1)


    error_vector_2 = protein_Venus_exp .- protein_Venus_sim 


    #protein BFP- not having units..
    # protein_BFP_std_exp[1] = 1.0
    # tmp_arr2 = 1.0./((protein_BFP_std_exp).^2)
    # W_prot2 = diagm(tmp_arr2)

    error_vector_3 = protein_BFP_exp .- protein_BFP_sim 

    #mRNA BFP 

    error_vector_4 = mRNA_BFP_exp .- mRNA_BFP_sim

  #  # upper bounds on GntR, S28 and AS28 protein concentrations
    # GntR_prot_for_bound = XSIM[:,11]
    # S28_prot_for_bound = XSIM[:,12]
    # AS28_prot_for_bound = XSIM[:,13]
    # UB = 2.0 # set a realistic upper bound
    # MV_GntR = maximum(GntR_prot_for_bound) # grab max from simulation
    # MV_S28 = maximum(S28_prot_for_bound) # grab max from simulation
    # MV_AS28 = maximum(AS28_prot_for_bound) # grab max from simulation
    
    # error_term_GntR = 1*max(0,(MV_GntR - UB))
    # error_term_S28 = 1*max(0,(MV_S28 - UB))
    # error_term_AS28 = 1*max(0,(MV_AS28 - UB))


    # error_term_array[2] = transpose(error_vector_2a)*error_vector_2a + transpose(error_vector_2b)*error_vector_2b this includes mRNA also- for now im not doin it- the product takes care of the square 
    error_term_array[1] = 1.0*transpose(error_vector_1)*error_vector_1 
    error_term_array[2] = 1.0*transpose(error_vector_2)*error_vector_2
    error_term_array[3] = 1.0*transpose(error_vector_3)*error_vector_3
    error_term_array[4] = 1.0*transpose(error_vector_4)*error_vector_4
    # error_term_array[5] = error_term_GntR
    # error_term_array[6] = error_term_S28
    # error_term_array[7] = error_term_AS28

    # # ===================================================================================================== #
	#error_total = sum(error_term_array)

    # # return -
    return error_term_array
end

# # Evaluates the objective function values -
function local_refinement_step(path_to_data_dir, parameter_array; sigma=0.05, iteration_max=100)


    #inner functions - computes a value total for each error term and squishes it one single value 
    function _compute_error_total(objective_array,W)
        value = transpose(objective_array)*W*(objective_array) 
        return value[1] 
    end

    # initialize -
    
    l= length(parameter_array) 
    number_of_parameters=l 
    BIG = 1e10

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 10.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

#     # wght array -  
    W = diagm(ones(4)) #this is needed for making it one by one 


#     # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

#     # Set some values for the weights, gene lengths etc -
#     # model_data_dictionary = customize_data_dictionary(default_data_dictionary,host_type)
	model_data_dictionary = default_data_dictionary

#     # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

#     # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)

    # calculate the starting error -
    parameter_array_best = parameter_array
    error_array = BIG*ones(4) #
    error_array[1]=_compute_error_total(OF(parameter_array_best), W)
    
    # main refinement loop -
    iteration_counter = 1
    while (iteration_counter<iteration_max)

#         # take a step up -
        parameter_up = parameter_array_best.*(1 .+ sigma*rand(number_of_parameters))  
        parameter_up = check_parameter_bounds(parameter_up)

        # take a step down -
        parameter_down = parameter_array_best.*(1 .- sigma*rand(number_of_parameters))
        parameter_down = check_parameter_bounds(parameter_down)

        # Evaluate the obj function -
        error_array[2] = _compute_error_total(OF(parameter_up),W)
        error_array[3] = _compute_error_total(OF(parameter_down),W)

#         # Calculate a correction factor -
        a = error_array[2]+error_array[3] -2.0* error_array[1] 
        parameter_corrected = parameter_array_best
        if (a>0.0)
            amda = -0.5*(error_array[3]- error_array[2])/a 
            parameter_corrected = parameter_array_best .+ amda*rand(number_of_parameters)
            parameter_corrected = check_parameter_bounds(parameter_corrected)
            error_array[4] = _compute_error_total(OF(parameter_corrected), W) #so far okay
        end

        # Which step has the min error?
        min_index = argmin(error_array)
        if (min_index == 1)
            parameter_array_best = parameter_array_best
        elseif (min_index == 2)
            parameter_array_best = parameter_up
        elseif (min_index == 3)
            parameter_array_best = parameter_down
        elseif (min_index == 4)
            parameter_array_best = parameter_corrected
        end

        # Update the local error
        error_array[1] = error_array[min_index]

        @show iteration_counter,error_array[min_index]

        # update local counter -
        iteration_counter = iteration_counter + 1
    end

    return parameter_array_best
end

function check_parameter_bounds(parameter_array)

    # setup paramter bounds - User is supposed to give Parameter_bounds array
    pvec_bounds = [

        # dG's - log values (the exp of these values will be the W used in)

        0.1 5.0  	;   # 1     W_RNAP_P70
        0.1 5.0    	;   # 2     W_RNAP_P28
        0.1 5.0   	;   # 3     W_RNAP_mP70 was initially 5 for limit
        -5.0  -0.1       ;   # 4     W_S70_RNAP_P70
        -5.0  -0.1       ;   # 5     W_S70_RNAP_mP70
        -5.0  -0.1      ;   # 6     W_S28_RNAP_P28
        -5.0   -0.1       ;   # 7     W_GntR_mP70
        1.0 1.0000001       ;   # 8     W_AS28_S28_P28- MAKING SURE IT bounds within 0.01 changed lb from -3.5 to -5
       
       

        # binding parameters -
        0.5 10.0            ;   # 9     n_S70_RNAP_P70
        0.01 100         ;   # 10     K_S70_RNAP_P70
		0.5 10.0            ;   # 11     n_S70_RNAP_mP70
        0.01 100         ;   # 12     K_S70_RNAP_mP70
        0.5 10.0            ;   # 13     n_GntR_mP70
        1e-6 1e-3        ;   # 14     K_GntR_mP70
        0.5 10.0            ;   # 15     n_S28_RNAP_P28
        1e-3 100.0         ;   # 16     K_S28_RNAP_P28        
		0.5 10.0            ;   # 17    n_AS28_S28 
        1e-3 100.0         ;   # 18    K_AS28_S28
 

        # time constants -
		0.1 10.0         ;	# 19	    mRNA_GntR
		0.05 10.0         ;	# 20	    mRNA_S28
        0.1 10.0         ;	# 21	    mRNA_AS28
		0.01 10.0         ;	# 22	    mRNA_Venus
        0.1 10.0         ;	# 23	    mRNA_BFP
		0.1 10.0         ;	# 24	    protein_GntR
		0.1 10.0         ;	# 25	    protein_S28
        0.5 100.0         ;	# 26	    protein_AS28
		0.01 100.0         ;	# 27	    protein_Venus
        0.01 100.0         ;	# 28	    protein_BFP
       
        # degradation mods -

        1e-2 100.0         ;	# 29	    mRNA_GntR
		1e-2 100.0         ;	# 30	    mRNA_S28
        1e-2 100.0         ;	# 31	    mRNA_AS28
		1e-2 100.0         ;	# 32	    mRNA_Venus
        1e-2 100.0         ;	# 33	    mRNA_BFP
		1e-2 100.0         ;	# 34	    protein_GntR
		1e-2 100.0         ;	# 35	    protein_S28
        1e-2 100.0         ;	# 36	    protein_AS28
		1e-2 100.0         ;	# 37	    protein_Venus
        1e-2 100.0         ;	# 38	    protein_BFP
	

         # w -
        2.0 12.0           ;   # 39    translation capacity half-life

        # KL and KX values -
        10.0 500.0         ;   # 40    KL in muM
        0.05 0.50001         ;   # 41    KX in muM. bound for now.

		# n_gluconate_GntR-
		1.0 5.0					; # 42

		# K_gluconate_GntR-
		0.01 5.0				; # 43 mM

        # sigma 70 half life
        10.0 30.0            ;   # 44    sigma 70 half life
    ];
	
    # tmp -
    pvec_initial = parameter_array

    # check bounds -
  
    number_of_parameters = length(pvec_initial) 

    for parameter_index = 1:number_of_parameters

        # what is the parameter value?
        p_i = pvec_initial[parameter_index]

        # is p_i outside of the bounds?
        lb_value = pvec_bounds[parameter_index,1]
        ub_value = pvec_bounds[parameter_index,2]

        if (p_i<lb_value)
            pvec_initial[parameter_index,1] = lb_value
        end

        if (p_i>ub_value)
            pvec_initial[parameter_index,1] = ub_value
        end
    end

    # return -
    return pvec_initial
end

function neighbor_function(parameter_array; sigma=0.05)

    # setup -
    number_of_parameters = length(parameter_array)

    # calculate new parameter array -
    new_parameter_array = parameter_array.*(1 .+ sigma*randn(number_of_parameters))

    # check the bounds and return -
    return check_parameter_bounds(new_parameter_array)
end

function cooling_function(temperature)

  # define my new temperature -
  alpha = 0.9
  return alpha*temperature
end

function acceptance_probability_function(rank_array,temperature)
    return (exp(-rank_array[end]/temperature))
end

function main(path_to_data_dir::String, initial_parameter_array::Array{Float64,1}; rank_cutoff::Int64=4, maximum_number_of_iterations::Int64=100) 

    # load the default data_dictionary -
    time_start = 0.0
    time_stop = 10.0
    time_step_size = 0.01

    # what is the host_type?
    host_type = :cell_free
    number_of_parameters = length(initial_parameter_array)

    # path to parameters -
    path_to_biophysical_constants_file = "./CellFree.json"

    # Load the data dictionary (uses the default biophysical_constants file)
    default_data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)

    # Set some values for the weights, gene lengths etc -
   
	model_data_dictionary = deepcopy(default_data_dictionary)
    # load the experimental data -
    exp_data_dictionary = load_experimental_data_dictionary(path_to_data_dir)

    # setup the functions -
    OF(P) = objective_function(P,time_start,time_step_size,time_stop,model_data_dictionary,exp_data_dictionary)
    NF(P) = neighbor_function(P;sigma=0.01)

    # make call to POETs -
    (EC,PC,RA) = estimate_ensemble(OF,NF,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=rank_cutoff,maximum_number_of_iterations=maximum_number_of_iterations)

    # return -
    return (EC,PC,RA)
end

##
# import initial param guess Array
initial_param_import = readdlm("initial_param_NEW_01.csv")
pvec_initial = initial_param_import[:,1]
# pvec_initial = [

# 	# dG's - Kj/mol
# 	5.0    	;   # 1     W_RNAP_P70- From Abhi Work of GntR Model- were initially 2.5 1.5 and 1 respectively, because for the model, it will be exp(-W)
# 	5.0   	;   # 2     W_RNAP_P28 
#     5.0    ;   # 3     W_RNAP_mP70 
#     -5.0    ;   # 4     W_S70_RNAP_P70 
#     -2.5    ;   # 5     W_S70_RNAP_mP70, s70 should have more trouble binding to operator region
#     -5.0   ;   # 6     W_S28_RNAP_P28
#     -2.5    ;   # 7     W_GntR_mP70 - from Abhi Work of GntR model
# 	1.0    ;   # 8     W_AS28_S28_P28 

# 	# binding parameters -
# 	1.0            ;   # 9     n_S70_RNAP_P70 all n were 1 and K were 30
# 	0.001       ;   # 10         K_S70_RNAP_P70, was 100
# 	1.0            ;   # 11     n_S70_RNAP_mP70
# 	0.001         ;   # 12     K_S70_RNAP_mP70 was 100
#     1.0            ;   # 13     n_GntR_mP70
# 	0.01         ;   # 14     K_GntR_mP70
# 	1.0            ;   # 15     n_S28_RNAP_P28
# 	0.00001         ;   # 16     K_S28_RNAP_P28
# 	2.0            ;   # 17     n_AS28_S28
# 	1         ;   # 18    K_AS28_S28

# 	# time constants - 
# 	1.0         ;	#19	        #mRNA_GntR
# 	1.0         ;	# 20	    mRNA_S28
#     1.0         ;	# 21	    mRNA_AS28
# 	1.0         ;	# 22	    mRNA_Venus
#     1.0         ;   # 23        mRNA_BFP

# 	1.0          ;	# 24	    protein_GntR
# 	1.0         ;	# 25	    protein_S28
#     1.0         ;   # 26	    protein_AS28
#     1.0         ;	# 27	    protein_Venus changed from 
#     1.0         ;   # 28	    protein_BFP
    

# 	# degradation mods -
# 	1.0         ;	# 29	    mRNA_GntR
# 	1.0         ;	# 30	    mRNA_S28
#     1.0         ;	# 31	    mRNA_AS28
# 	1.0         ;	# 32	    mRNA_Venus
#     1.0         ;   # 33        mRNA_BFP

#     1.0          ;	# 34	    protein_GntR 
# 	1.0         ;	# 35	    protein_S28
#     1.0         ;   # 36	   protein_AS28-changed from 1 to 1.2 and 0.6 to 0.8
#     1.0         ;	# 37	    protein_Venus Changed from 8,7 to 9.,9
#     1.0         ;   # 38	    protein_BFP


# 	 # w -
# 	6.0           ;   # 39    translation capacity half-life

# 	# KL value -
# 	100.0         ;   # 40    KL in muM
	
#     # KX value -
# 	0.05         ;   # 41    KL in muM
# 	# n_gluconate_GntR-
# 	1.0					; # 42

# 	# K_gluconate_GntR-
# 	1.0					; #43 mM

#     # transcription capacity terms-
#     #0.0; # 49 decay (hours)-
#     #0.6; # 50 slope

#     # translation capacity terms-
#     #0.0; #51 decay (hours)
#     #0.1; #52 slope


# ];

# setup -
path_to_data_dir = "$(pwd())"
# pV = neighbor_function(pvec_initial; sigma=0.05)
pV=pvec_initial
EC = 0
PC = 0
RA = 0

##

# execute -
number_of_trials = 10
for trial_index = 1:number_of_trials

    global pV
    global EC
    global PC
    global RA
 

#    # do a local step - 
#     if (mod(trial_index,2) == 0)

#         # find the lowest score pV -
#         sum_error_array = sum(EC,dims=1)
#         best_p_index = argmin(vec(sum_error_array))
#         pV_best = PC[:,best_p_index]

#         # local refine -
#         pV = local_refinement_step(path_to_data_dir, pV_best; iteration_max=20)
#     end

    # main -
    (EC,PC,RA) = main(path_to_data_dir, vec(pV); rank_cutoff=4,maximum_number_of_iterations=20)

    # dump results to disk -
    fname = "./simulated/poets_ensemble/RA_T$(trial_index).dat"
    writedlm(fname,RA)
    fname = "./simulated/poets_ensemble/EC_T$(trial_index).dat"
    writedlm(fname,EC)
    fname = "./simulated/poets_ensemble/PC_T$(trial_index).dat"
    writedlm(fname,PC)

    @show trial_index
end


