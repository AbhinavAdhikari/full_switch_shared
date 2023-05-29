# ----------------------------------------------------------------------------------- #
# Copyright (c) 2019 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

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

# include -
include("Include.jl")

# Script to solve the balance equations -
time_start = 0.0
time_stop = 16.0
time_step_size = 0.01

# what is the host_type?
host_type = :cell_free

# path to parameters -
path_to_biophysical_constants_file = "./CellFree.json"

PC = readdlm("./poets_ensemble_W_test/PC_T5.dat")

# plot experimental data protein- Venis
prot_data = CSV.read("Protein_data_Venus.csv",DataFrame)
time = prot_data[!,"time(h)"]
Venus = prot_data[!,"mean_10mM"]
stdev_prot = prot_data[!,"sd_10mM"]


p2 = Plots.scatter(time,Venus,marker = "o",markercolor = "black", markersize = 3,yerror=stdev_prot, label = "", xlabel="Time (hr)",ylabel = "Concentration (nM)")

Plots.plot(p2)

for i in 1:size(PC)[2]

    parameter_array = PC[:,i]

    # Load the data dictionary (uses the default biophysical_constants file)
    data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
    #print(data_dictionary)

    R = data_dictionary["R"]
    T_K = data_dictionary["T_K"]
    # compute W -


    # Update the data dictionary
    control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]

    control_parameter_dictionary["W_RNAP_P70"] = exp(-1*(parameter_array[1]/100)/(R*T_K))
    control_parameter_dictionary["W_RNAP_P28"] = exp(-1*(parameter_array[2]/100)/(R*T_K))
    control_parameter_dictionary["W_RNAP_mP70"] = exp(-1*(parameter_array[3]/100)/(R*T_K))
    control_parameter_dictionary["W_S70_RNAP_P70"] = exp(-1*(parameter_array[4]/100)/(R*T_K))
	control_parameter_dictionary["W_S70_RNAP_mP70"] = exp(-1*(parameter_array[5]/100)/(R*T_K))
	control_parameter_dictionary["W_S28_RNAP_P28"] = exp(-1*(parameter_array[6]/100)/(R*T_K))
	control_parameter_dictionary["W_GntR_mP70"] = exp(-1*(parameter_array[7]/100)/(R*T_K))
    control_parameter_dictionary["W_GntR_gluconate_protein"] = exp(-1*(parameter_array[8]/100)/(R*T_K))
	control_parameter_dictionary["W_AS28_S28"] = exp(-1*(parameter_array[9]/100)/(R*T_K))
	

     data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

    binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
    binding_parameter_dictionary["n_S70_RNAP_GntR"]=parameter_array[10]
	binding_parameter_dictionary["K_S70_RNAP_GntR"]=parameter_array[11]
	binding_parameter_dictionary["n_S70_RNAP_S28"]=parameter_array[12]
	binding_parameter_dictionary["K_S70_RNAP_S28"]=parameter_array[13]
	binding_parameter_dictionary["n_S70_RNAP_AS28"]=parameter_array[14]
    binding_parameter_dictionary["K_S70_RNAP_AS28"]=parameter_array[15]
	binding_parameter_dictionary["n_S70_RNAP_Venus"]=parameter_array[16]
    binding_parameter_dictionary["K_S70_RNAP_Venus"]=parameter_array[17]
	binding_parameter_dictionary["n_GntR_mP70_AS28"]=parameter_array[18]
	binding_parameter_dictionary["K_GntR_mP70_AS28"]=parameter_array[19]
	binding_parameter_dictionary["n_GntR_mP70_Venus"]=parameter_array[20]
	binding_parameter_dictionary["K_GntR_mP70_Venus"]=parameter_array[21]
	binding_parameter_dictionary["n_S28_RNAP_BFP"]=parameter_array[22]
	binding_parameter_dictionary["K_S28_RNAP_BFP"]=parameter_array[23]
	binding_parameter_dictionary["n_AS28_S28_BFP"]=parameter_array[24]
	binding_parameter_dictionary["K_AS28_S28_BFP"]=parameter_array[25]
    data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

    data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

    time_constant_modifier_array = [
        0.0							;	# 1	GntR
        0.0							;	# 2	S28
        0.0							;	# 3	AS28
        0.0                         ;   # 4 Venus
        0.0                         ;   # 5 BFP
        parameter_array[26]            ;  #  6 mRNA GntR
        parameter_array[27]            ;  #  7 mRNA S28
        parameter_array[28]            ;  #  8 mRNA AS28
        parameter_array[29]            ;  #9 mRNA Venus
        parameter_array[30]            ; #10 mRNA BFP
        parameter_array[31]            ;  #11 Protein GntR
        parameter_array[32]            ; #12 Protein S28
        parameter_array[33]            ; #13 Protein AS28
        parameter_array[34]            ; #14 Protein Venus
        parameter_array[35]            ; #15 Protein BFP
    ]

    data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

    degradation_modifier_array = [
        0.0							;	# 1	GntR
        0.0							;	# 2	S28
        0.0							;	# 3	AS28
        0.0                         ;   # 4 Venus
        0.0                         ;   # 5 BFP
        parameter_array[36]            ;  #  6 mRNA GntR
        parameter_array[37]            ;  #  7 mRNA S28
        parameter_array[38]            ;  #  8 mRNA AS28
        parameter_array[39]            ;  #9 mRNA Venus
        parameter_array[40]            ; #10 mRNA BFP
        parameter_array[41]            ;  #11 Protein GntR
        parameter_array[42]            ; #12 Protein S28
        parameter_array[43]            ; #13 Protein AS28
        parameter_array[44]        ;     #14 Protein Venus
        parameter_array[45]            ; #15 Protein BFP
    ]

    data_dictionary["degradation_modifier_array"] = degradation_modifier_array

    # update the translation time -
    data_dictionary["half_life_translation_capacity"] = parameter_array[46]

    biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["translation_saturation_constant"] = parameter_array[47]
    data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

    # update the transcription time -
    biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["transcription_elongation_rate"] = parameter_array[48]
    data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary
    # lastly, update KX -
    biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
    biophysical_constants_dictionary["transcription_saturation_constant"] = parameter_array[49]
    data_dictionary["biophysical_constants_dictionary"] = biophysical_constants_dictionary

    species_symbol_type_array = data_dictionary["species_symbol_type_array"]
    protein_coding_length_array = data_dictionary["protein_coding_length_array"]
    gene_coding_length_array = data_dictionary["gene_coding_length_array"]
    time_constant_modifier_array = data_dictionary["time_constant_modifier_array"]
    initial_condition_array = data_dictionary["initial_condition_array"]

    # # get gene IC -
    idx_gene = findall(x->x==:gene,species_symbol_type_array)
    gene_abundance_array = initial_condition_array[idx_gene]

    # Precompute the translation parameters -
    translation_parameter_array = precompute_translation_parameter_array(biophysical_constants_dictionary, protein_coding_length_array, time_constant_modifier_array,host_type)
    data_dictionary["translation_parameter_array"] = translation_parameter_array

    # Precompute the kinetic limit of transcription -
    transcription_kinetic_limit_array = precompute_transcription_kinetic_limit_array(biophysical_constants_dictionary, gene_coding_length_array, gene_abundance_array, time_constant_modifier_array, host_type)
    data_dictionary["transcription_kinetic_limit_array"] = transcription_kinetic_limit_array

    # Dilution degrdation matrix -
    dilution_degradation_matrix = build_dilution_degradation_matrix(biophysical_constants_dictionary, species_symbol_type_array, degradation_modifier_array)
    data_dictionary["dilution_degradation_matrix"] = dilution_degradation_matrix

    # Solve the model equations -
    (T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)

    p2 = Plots.plot!(T,X[:,14],linewidth=2, label = "")

end

Plots.plot(p2)
Plots.savefig("./Ensemble_venus.pdf")