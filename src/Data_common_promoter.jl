# ----------------------------------------------------------------------------------- #
# Copyright (c) 2023 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2021-04-29T20:40:53.886
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{String,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function build_data_dictionary(time_span::Tuple{Float64,Float64,Float64}, path_to_biophysical_constants_file::String = "./CellFree.json", host_type::Symbol = :cell_free)::Dict{String,Any}

	# load the biophysical_constants dictionary
	biophysical_constants_dictionary = build_biophysical_dictionary(path_to_biophysical_constants_file, host_type)

	# stoichiometric_matrix and dilution_matrix -
	stoichiometric_matrix = readdlm("./Network.dat")

	# number of states, and rates -
	(number_of_states,number_of_rates) = size(stoichiometric_matrix)

	# array of species types -
	species_symbol_type_array = [
		:gene	;	# 1	GntR
		:gene	;	# 2	S28
		:gene	;	# 3	AS28 
		:gene	;	# 4	Venus 
		:gene	;	# 5	BFP 
		:mrna	;	# 6	mRNA_GntR
		:mrna	;	# 7	mRNA_S28		
		:mrna	;	# 8	mRNA_AS28
		:mrna	;	# 9	mRNA_Venus
		:mrna	;	# 10 mRNA_BFP
		:protein	;	# 11	protein_GntR
		:protein	;	# 12	protein_S28
		:protein	;	# 13	protein_AS28
		:protein	;	# 14	protein_Venus
		:protein	;	# 15	protein_BFP
	]

	# we need to store the species symbol array for later -
	biophysical_constants_dictionary["species_symbol_type_array"] = species_symbol_type_array

	# array of gene lengths -
	gene_coding_length_array = [
		996.0	;	# 1	GntR
		720.0	;	# 2	S28 
		294.0	;	# 3	AS28 
		717.0   ;   # 4 Venus
		714.0   ;   # 5 BFP
	]

	# array of mRNA coding lengths -
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 6	1	mRNA_GntR
		gene_coding_length_array[2]	;	# 7	2	mRNA_S28
		gene_coding_length_array[3]	;	# 8	3	mRNA_AS28
		gene_coding_length_array[4]	;	# 9	4	mRNA_Venus
		gene_coding_length_array[5]	;	# 10 5	mRNA_BFP
	]

	# array of mRNA coding lengths -
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 11	1	protein_GntR
		round((0.33)*mRNA_coding_length_array[2])	;	# 12	2	protein_S28
		round((0.33)*mRNA_coding_length_array[3])	;	# 13	3	protein_AS28
		round((0.33)*mRNA_coding_length_array[4])	;	# 14	4	protein_Venus
		round((0.33)*mRNA_coding_length_array[5])	;	# 15	5	protein_BFP
	]

	# array of intracellular gene copy numbers - in cell free- this is just the concentration in uM
	gene_abundance_array = [
		0.010	;	#1 GntR
		0.007	;	#2 S28
		0.007   ;   #3 AS28
		0.007   ;   #4 Venus
		0.010	;	#5 BFP
	]

	# initial condition array - #concentration should be in uM
	initial_condition_array = [
		gene_abundance_array[1]	;	# 1	GntR
		gene_abundance_array[2]	;	# 2	S28
		gene_abundance_array[3]	;	# 3	AS28
		gene_abundance_array[4]	;	# 4	Venus
		gene_abundance_array[5]	;	# 5	BFP

		0.0	;	# 6	 mRNA_GntR
		0.0	;	# 7	 mRNA_S28
		0.0	;	# 8	 mRNA_AS28
		0.0	;	# 9	 mRNA_Venus
		0.0	;	# 10 mRNA_BFP
		0.0;    # 11 protein_GntR
		0.0	;	# 12 protein_S28 
		0.0	;	# 13 protein AS28
		0.0;    # 14 protein Venus
		0.0;    #  15 protein BFP
		# translation capacity -
		100.0	;	# 16 translation capacity
		35E-3 ; # 17 initial sigma 70 concentration in uM.
	]

	binding_parameter_dictionary = Dict{String,Float64}() #set an initial value. Optimization later will adjust it
	binding_parameter_dictionary["n_S70_RNAP_P70"]=1
	binding_parameter_dictionary["K_S70_RNAP_P70"]=0.005#I used Abhi's work for setting initial guess, only this is in uM
	binding_parameter_dictionary["n_S70_RNAP_mP70"]=1
    binding_parameter_dictionary["K_S70_RNAP_mP70"]=0.05
	binding_parameter_dictionary["n_GntR_mP70"]=1 #binds to mP70, WAS 2
	binding_parameter_dictionary["K_GntR_mP70"]=0.05 #say 
	binding_parameter_dictionary["n_S28_RNAP_P28"]=1 
	binding_parameter_dictionary["K_S28_RNAP_P28"]=2
	binding_parameter_dictionary["n_AS28_S28"]=1 #WAS 1 AND 50- place holder for Km
	binding_parameter_dictionary["K_AS28_S28"]=0.209

	
	# Alias the control function parameters -
	control_parameter_dictionary = Dict{String,Float64}()
	control_parameter_dictionary["W_RNAP_P70"] = 0.005 #while the RNAP can bind to any DNA sequence, the energy with which it binds will depend on DNA sequence 
	control_parameter_dictionary["W_RNAP_P28"] = 0.005 
	control_parameter_dictionary["W_RNAP_mP70"] = 0.005
	control_parameter_dictionary["W_S70_RNAP_P70"] = 3 #based on parameter fit, was 1 initially, try changing it to 2.7
	control_parameter_dictionary["W_S70_RNAP_mP70"] = 3 #based on p fit, was 0.25
	control_parameter_dictionary["W_S28_RNAP_P28"] = 3 #
	control_parameter_dictionary["W_GntR_mP70"] = 5 #
	# control_parameter_dictionary["W_AS28_S28_P28"] = 1 #WAS 1 
	
	# D- Gluconate parameter values
	gluconate_parameter_dictionary = Dict{String,Float64}()
	gluconate_parameter_dictionary["gluconate_concentration"]=1e-5 # units mM- should be 1e-5 mM for fittng the model (The range works from 1e2mM to 1e-5mM)
	gluconate_parameter_dictionary["n_gluconate_GntR"] = 1 #units are in mM
	gluconate_parameter_dictionary["K_gluconate_GntR"] = 1 # mM
	gluconate_parameter_dictionary["Protein_sigma_70"] = 35E-3 #uM units

# degradation modifiers - #for all zero degradation modifier, we get a michelis curve.- this is for θ 
	degradation_modifier_array = [
		0.0	;	# 1	GntR
		0.0	;	# 2	S28
		0.0	;	# 3	AS28
		0.0	;	# 4	Venus
		0.0	;	# 5	BFP
		1.0	;	# 6	mRNA_GntR
		1.0	;	# 7	mRNA_S28
		1.0	;	# 8	mRNA_AS28
		1.0 ;   # 9 mRNA_Venus
		1.0 ;   #10 mRNA_BFP
		10.0	;	#11	protein_GntR
		12.0	;	#12	protein_S28
		8.0	;	#13	protein_AS28 INCREASING THIS VALUE REDUCES THE VALUE/RANGE OF ACTOR3/Simulation val
		20.0 ;   #14 protein_Venus
		20.0 ;   #15 protein_BFP #Was 1.5 ish
	]


	# time constant modifiers - how do i decide how to choose
	time_constant_modifier_array = [
		0.0	;	# 1	GntR
		0.0	;	# 2	S28
		0.0	;	# 3	AS28
		0.0 ;   # 4 Venus
		0.0	;	# 5 BFP
		1.0	;	# 6	 mRNA_GntR
		1.0	;	# 7	 mRNA_S28
		1.0	;	# 8	 mRNA_AS28
		1.0	;	# 9	 mRNA_Venus
		1.0	;	# 10 mRNA_BFP
		4.0	;	# 11 protein_GntR
		5.0	;	# 12 protein_S28
		5.0	;	# 13 protein_AS28 # Decreasing this value increases the actor3 range of values and broadens it. - basically increases the sim value
		2.0	;	# 14 protein_Venus
		2.0	;	# 15 protein_BFP
	]

	

	#Dilution degrdation matrix -
	dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary,species_symbol_type_array,degradation_modifier_array)

	# Precompute the translation parameters -
	translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array, host_type)

	# Precompute the kinetic limit of transcription -
	transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)

	# Parameter name index array - #all these parameters are used to check for the morris sensitivity- 
	parameter_name_mapping_array = [
		"W_RNAP_P70"; #1
		"W_RNAP_P28" ;#2
		"W_RNAP_mP70" ;#3
		"W_S70_RNAP_P70"; #4
		"W_S70_RNAP_mP70"; #5
		"W_S70_RNAP_P28" ;#6
		"W_GntR_mP70" ;#7
		"W_AS28_S28_P28"; #9
		"n_S70_RNAP_P70"; #10
		"K_S70_RNAP_P70";#11
		"n_S70_RNAP_mP70" ;#14
		"K_S70_RNAP_mP70" ;#15
		"n_GntR_mP70" ;#18
		"K_GntR_mP70" ;#19
		"n_S28_RNAP_P28" ;#22
		"K_S28_RNAP_P28" ;#23
		"n_AS28_S28" ;#24
		"K_AS28_S28" ;#25
		"n_gluconate_GntR" ;#26
		"K_gluconate_GntR" ;#27
		"RNAPII_concentration"; #28		
		"ribosome_concentration"	; #29	
		"degradation_constant_mRNA"	;	# 30 
		"degradation_constant_protein"	;	#31 
		"kcat_transcription"	;	# 32
		"kcat_translation"	;	# 33
		"saturation_constant_transcription"	;	#34 
		"saturation_constant_translation"	;	# 35
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{String,Any}()
	data_dictionary["number_of_states"] = number_of_states
	data_dictionary["species_symbol_type_array"] = species_symbol_type_array
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix
	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array
	data_dictionary["translation_parameter_array"] = translation_parameter_array
	data_dictionary["degradation_modifier_array"] = degradation_modifier_array
	data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array
	data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary
	data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #

	data_dictionary["R"] = 8.314 			# J mol^-1 K^-1
	data_dictionary["T_K"] = 273.15 + 30.0 	# K
	data_dictionary["half_life_translation_capacity"] = 8.0	 # hr
	data_dictionary["base_line_weight"] = (1.0/1.0)
	data_dictionary["sigma_70_halflife"] = 20.0 #hr


	return data_dictionary
end
