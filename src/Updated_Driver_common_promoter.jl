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
#include("Include.jl")- for now, no need
include("Include.jl")

# Script to solve the balance equations -
time_start = 0.0
time_stop = 10.0
time_step_size = 0.01

# what is the host_type?
host_type = :cell_free

# path to parameters -
path_to_biophysical_constants_file = "./CellFree.json"

# Load the data dictionary (uses the default biophysical_constants file)
data_dictionary = build_data_dictionary((time_start,time_stop,time_step_size), path_to_biophysical_constants_file, host_type)
#print(data_dictionary)

R = data_dictionary["R"]
T_K = data_dictionary["T_K"]
# compute W -

PC = readdlm("./simulated/poets_ensemble/PC_T6.dat")
average_params = mean(PC,dims=2)
poets_params=average_params
#override
# poets_params = readdlm("initial_param_NEW_01.csv")


# gluconate GntR binding parameters
gluconate_parameter_dictionary = data_dictionary["gluconate_parameter_dictionary"]
gluconate_parameter_dictionary["n_gluconate_GntR"] = poets_params[42]
gluconate_parameter_dictionary["K_gluconate_GntR"] = poets_params[43]




gluconate_parameter_dictionary["gluconate_concentration"] = 1e-3


data_dictionary["gluconate_parameter_dictionary"] = gluconate_parameter_dictionary

# update sigma 70 half_life
data_dictionary["sigma_70_half_life"] = poets_params[44]

# Update the data dictionary
control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
control_parameter_dictionary["W_RNAP_P70"] = exp(-1*(poets_params[1]))
control_parameter_dictionary["W_RNAP_P28"] = exp(-1*(poets_params[2]))*1
control_parameter_dictionary["W_RNAP_mP70"] = exp(-1*(poets_params[3]))
control_parameter_dictionary["W_S70_RNAP_P70"] = exp(-1*(poets_params[4]))
control_parameter_dictionary["W_S70_RNAP_mP70"] = exp(-1*(poets_params[5]))
control_parameter_dictionary["W_S28_RNAP_P28"] = exp(-1*(poets_params[6]))*1
control_parameter_dictionary["W_GntR_mP70"] = exp(-1*(poets_params[7]))
control_parameter_dictionary["W_AS28_S28_P28"] = exp(-1*(poets_params[8]))*1

data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary

    binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]
    binding_parameter_dictionary["n_S70_RNAP_P70"]=poets_params[9]
	binding_parameter_dictionary["K_S70_RNAP_P70"]=poets_params[10]
	binding_parameter_dictionary["n_S70_RNAP_mP70"]=poets_params[11]
	binding_parameter_dictionary["K_S70_RNAP_mP70"]=poets_params[12]
	binding_parameter_dictionary["n_GntR_mP70"]=poets_params[13]
    binding_parameter_dictionary["K_GntR_mP70"]=poets_params[14]  
	binding_parameter_dictionary["n_S28_RNAP_P28"]=poets_params[15] 
    binding_parameter_dictionary["K_S28_RNAP_P28"]=poets_params[16]
	binding_parameter_dictionary["n_AS28_S28"]=poets_params[17]
	binding_parameter_dictionary["K_AS28_S28"]=poets_params[18]

data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary

time_constant_modifier_array = [
    0.0							;	# 1	GntR
    0.0							;	# 2	S28
    0.0							;	# 3	AS28
    0.0                         ;   # 4 Venus
    0.0                         ;   # 5 BFP
    poets_params[19]            ;  #  6 mRNA GntR
    poets_params[20]            ;  #  7 mRNA S28
    poets_params[21]            ;  #  8 mRNA AS28
    poets_params[22]          ;  #9 mRNA Venus
    poets_params[23]             ; #10 mRNA BFP
    poets_params[24]            ;  #11 Protein GntR
    poets_params[25]            ; #12 Protein S28
    poets_params[26]           ; #13 Protein AS28
    poets_params[27]             ; #14 Protein Venus
    poets_params[28]            ; #15 Protein BFP
]

data_dictionary["time_constant_modifier_array"] = time_constant_modifier_array

degradation_modifier_array = [
    0.0							;	# 1	GntR
    0.0							;	# 2	S28
    0.0							;	# 3	AS28
    0.0                         ;   # 4 Venus
    0.0                         ;   # 5 BFP
    poets_params[29]            ;  #  6 mRNA GntR
    poets_params[30]            ;  #  7 mRNA S28
    poets_params[31]            ;  #  8 mRNA AS28
    poets_params[32]           ;  #9 mRNA Venus
    poets_params[33]           ; #10 mRNA BFP
    poets_params[34]            ;  #11 Protein GntR
    poets_params[35]            ; #12 Protein S28
    poets_params[36]            ; #13 Protein AS28
    poets_params[37]        ;     #14 Protein Venus
    poets_params[38]           ; #15 Protein BFP
    
]

data_dictionary["degradation_modifier_array"] = degradation_modifier_array

# update the translation time -
data_dictionary["half_life_translation_capacity"] = poets_params[39]

biophysical_constants_dictionary = data_dictionary["biophysical_constants_dictionary"]
biophysical_constants_dictionary["translation_saturation_constant"] = poets_params[40]
biophysical_constants_dictionary["trancription_saturation_constant"] = poets_params[41]
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

# Plot protein
# plot sim data
p1 = Plots.plot(T,X[:,14], label = "Venus_Model_Prot",xlabel="Time (hr)",ylabel = "Concentration (μM)",linewidth=3)

# plot experimental data Venus
prot_data = CSV.read("./data/FINAL_FULL_MODEL_FIT_00001mM.csv",DataFrame)
time = prot_data[!,"Time(h)"]
Venus = prot_data[!,"Venus_uM"]
stdev_prot = prot_data[!,"Venus_uM_STDERR"]

p2 = Plots.scatter!(time,Venus,label = "Venus_Exp_Prot",legend = :topleft,markercolor = "black", markersize = 3,yerror=stdev_prot)

# # # plot experimental data BFP 

prot_data2 = CSV.read("./data/FINAL_FULL_MODEL_FIT_00001mM.csv",DataFrame)
time = prot_data2[!,"Time(h)"]
BFP = prot_data2[!,"BFP_uM"]
stdev_prot1 = prot_data[!,"BFP_uM_STDERR"]

p3 = Plots.plot!(T,X[:,15], label = "BFP_Model_Prot",xlabel="Time (hr)",ylabel = "Concentration (μM)",linewidth=3)

p4 = Plots.scatter!(time,BFP,label = "BFP_Exp_Prot",legend = :topleft,markercolor = "blue", markersize = 3,yerror=stdev_prot1)
Plots.savefig("Fig_PC__commonPromoter_protein.pdf")


# load mRNA data 
mRNA_data_table = CSV.read("./data/LATEST_CSV_Venus_BFP_mRNA.csv",DataFrame)
T_mRNA = mRNA_data_table[!,"Time(h)"]
  mean_mRNA_Venus = mRNA_data_table[!,"Venus(nM)"]
  stdev_mRNA_Venus = mRNA_data_table[!,"Venus_STDERR"]
  mean_mRNA_BFP = mRNA_data_table[!,"BFP(nM)"]
  stdev_mRNA_BFP = mRNA_data_table[!,"BFP_STDERR"]
  mRNA_data_array = [T_mRNA mean_mRNA_Venus stdev_mRNA_Venus mean_mRNA_BFP stdev_mRNA_BFP]

  mRNA_Venus_exp = mRNA_data_array[:,2]  #nM
  mRNA_BFP_exp= mRNA_data_array[:,4] #nM

  p11 = Plots.plot(T,1000*X[:,9], label = "Venus_Model_mRNA",xlabel="Time (hr)",ylabel = "Concentration (nM)",linewidth=3)
  p12 = Plots.scatter!(T_mRNA,mean_mRNA_Venus,label = "VenusP_Exp_Prot",legend = :topleft,markercolor = "blue", markersize = 3,yerror=stdev_mRNA_Venus)
  p13 = Plots.plot!(T,1000*X[:,10], label = "BFP_Model_mRNA",xlabel="Time (hr)",ylabel = "Concentration (nM)",linewidth=3)
  p14 = Plots.scatter!(T_mRNA,mean_mRNA_BFP,label = "BFP_Exp_Prot",legend = :topleft,markercolor = "green", markersize = 3,yerror=stdev_mRNA_BFP)
  Plots.savefig("Fig_PC_commonPromoter_mRNA.pdf")

p111= Plots.plot(T,X[:,11],label = "GntR",color="black")
p112= Plots.plot!(T,X[:,12],label = "S28",color="blue")
p113= Plots.plot!(T,X[:,13],label = "AS28",color="green")
Plots.savefig("Fig_PC_commonPromoter_GntR_S28_AS28.pdf")


@show (X[:,17])