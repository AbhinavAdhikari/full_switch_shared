# ----------------------------------------------------------------------------------- # #Note for myself- The control array values are robust 
using PyCall
np = pyimport("numpy")
# Copyright (c) 2021 Varnerlab
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
# Function: calculate_transcription_control_array
# Description: Calculate the transcriptional control array at time t
# Generated on: 2021-04-29T20:40:54.294
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => 
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Transcriptional control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_transcription_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control 
	control_array = zeros(5) #5 genes, 5 CFS

	# Alias the species - 
	GntR = x[1]
	S28 = x[2]
	AS28 = x[3]
    Venus=  x[4]
	BFP= x[5]
	mRNA_GntR = x[6]
	mRNA_S28 = x[7]
	mRNA_AS28 = x[8]
	mRNA_Venus = x[9]
	mRNA_BFP = x[10]
	protein_GntR = x[11]
	protein_S28 = x[12]
	protein_AS28 = x[13]
	protein_Venus = x[14]
	protein_BFP = x[15]	
	actual_protein_sigma_70 = x[17]

	# Alias the binding parameters -
	binding_parameter_dictionary = data_dictionary["binding_parameter_dictionary"]

	# n_S70_RNAP_GntR=binding_parameter_dictionary["n_S70_RNAP_GntR"]
	# K_S70_RNAP_GntR=binding_parameter_dictionary["K_S70_RNAP_GntR"]
	n_S70_RNAP_P70=binding_parameter_dictionary["n_S70_RNAP_P70"]
	K_S70_RNAP_P70=binding_parameter_dictionary["K_S70_RNAP_P70"]
	# n_S70_RNAP_S28=binding_parameter_dictionary["n_S70_RNAP_S28"]
	# K_S70_RNAP_S28=binding_parameter_dictionary["K_S70_RNAP_S28"]
	# n_S70_RNAP_AS28=binding_parameter_dictionary["n_S70_RNAP_AS28"]
	# K_S70_RNAP_AS28=binding_parameter_dictionary["K_S70_RNAP_AS28"]
	n_S70_RNAP_mP70=binding_parameter_dictionary["n_S70_RNAP_mP70"]
	K_S70_RNAP_mP70=binding_parameter_dictionary["K_S70_RNAP_mP70"]
	# n_S70_RNAP_Venus=binding_parameter_dictionary["n_S70_RNAP_Venus"]
	# K_S70_RNAP_Venus=binding_parameter_dictionary["K_S70_RNAP_Venus"]
	# n_GntR_mP70_AS28=binding_parameter_dictionary["n_GntR_mP70_AS28"]
	# K_GntR_mP70_AS28=binding_parameter_dictionary["K_GntR_mP70_AS28"]
	# n_GntR_mP70_Venus=binding_parameter_dictionary["n_GntR_mP70_Venus"]
	# K_GntR_mP70_Venus=binding_parameter_dictionary["K_GntR_mP70_Venus"]
	n_GntR_mP70=binding_parameter_dictionary["n_GntR_mP70"]
	K_GntR_mP70=binding_parameter_dictionary["n_GntR_mP70"]
	n_S28_RNAP_P28=binding_parameter_dictionary["n_S28_RNAP_P28"] 
	K_S28_RNAP_P28=binding_parameter_dictionary["K_S28_RNAP_P28"]
	n_AS28_S28=binding_parameter_dictionary["n_AS28_S28"]
	K_AS28_S28=binding_parameter_dictionary["K_AS28_S28"] 




	# Alias the control function parameters -
	control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
	W_RNAP_P70 = control_parameter_dictionary["W_RNAP_P70"] 
	W_RNAP_P28=control_parameter_dictionary["W_RNAP_P28"]
	W_RNAP_mP70=control_parameter_dictionary["W_RNAP_mP70"] 
	W_S70_RNAP_P70=control_parameter_dictionary["W_S70_RNAP_P70"] 
	W_S70_RNAP_mP70=control_parameter_dictionary["W_S70_RNAP_mP70"] 
	W_S28_RNAP_P28=control_parameter_dictionary["W_S28_RNAP_P28"] 
	W_GntR_mP70=control_parameter_dictionary["W_GntR_mP70"]  #binds to mP70	
	# W_AS28_S28=control_parameter_dictionary["W_AS28_S28"] # not used in the current model
	

	#Since there are two GntR operators in which binding takes place:
	gluconate_parameter_dictionary = data_dictionary["gluconate_parameter_dictionary"]
	Gluc_conc = gluconate_parameter_dictionary["gluconate_concentration"]
	n_gluconate_GntR = gluconate_parameter_dictionary["n_gluconate_GntR"]
	K_gluconate_GntR = gluconate_parameter_dictionary["K_gluconate_GntR"]
	# protein_S70= gluconate_parameter_dictionary["Protein_sigma_70"]

	A1= Gluc_conc ^ n_gluconate_GntR 
	term1A=real(Gluc_conc^abs(n_gluconate_GntR))
	term1B=real(abs(K_gluconate_GntR)^abs(n_gluconate_GntR))
	B1= term1A+ term1B

	f_gluc_GntR = real(A1/B1)  #equivalent to f bound

	
	protein_GntR= (1-f_gluc_GntR)*protein_GntR #Remaining GntR unbound fraction (free from gluconate) #this protein changes but not actor set 2

	

	actor_set2= [protein_GntR]'

	actor=prod(protein_GntR)

 #################################################################################################################################################
	#Control Terms

	#1- For GntR

	# actors= protein_S70 #Since this is assumed to be a time invariant source, no need to take prod of it ie 3.5 e-3, because k for this is in uM 

	actors_s70 = [actual_protein_sigma_70]'
	actors = prod(actors_s70)
	
	term_num= real(actors ^real(n_S70_RNAP_P70))
	term_den_1= real(actors ^real(n_S70_RNAP_P70))
	term_den_2= (real(abs(K_S70_RNAP_P70))^(real(abs(n_S70_RNAP_P70))))

	term_den= term_den_1+ term_den_2

	f_S70_RNAP_GntR=  real(term_num/term_den)   


	

	num1= W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_GntR

	den1= 1+ W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_GntR


	control_array[1]= num1/den1  


	######################################

	#2 for S28 production- 

	# f_S70_RNAP_S28= real(actors ^n_S70_RNAP_S28) / ( real(actors ^n_S70_RNAP_S28) + real(abs(K_S70_RNAP_S28) ^abs(n_S70_RNAP_S28)) ) 

	# num2= W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_S28 #because S28 gene is under the control of P70 promoter
	# den2= 1+ W_RNAP_P70 + W_S70_RNAP_P70*f_S70_RNAP_S28


	control_array[2]= num1/den1 # identical to GntR (both P70 promoter)

	######################################

	#3 for AS28 production- 

    A=abs(actor)^abs(n_GntR_mP70) #actor is GntR as defined above
	B2=real(abs(K_GntR_mP70)^abs(n_GntR_mP70))
	B=A+B2

	f_GntR_mP70_AS28= real(abs(A/B)) 

	f_S70_RNAP_AS28= (actors ^(n_S70_RNAP_mP70)) / ( actors ^abs(n_S70_RNAP_mP70) + abs(K_S70_RNAP_mP70) ^abs(n_S70_RNAP_mP70))

	a1=W_RNAP_mP70 + W_S70_RNAP_mP70*f_S70_RNAP_AS28
	# a2=W_GntR_mP70*(f_gluc_GntR)*(f_GntR_mP70_AS28)

	# num3= a1 + a2

    # den3=1+a1+a2

	#control_array[3]=num3/den3
	control_array[3]=a1/(1+a1+W_GntR_mP70*f_GntR_mP70_AS28)


	######################################
	
	##4 Control function for Venus Production

	
	A4= abs(actor) ^ abs(n_GntR_mP70)
	B4= real(abs(actor) ^ abs(n_GntR_mP70) + abs(K_GntR_mP70) ^ abs(n_GntR_mP70))

    f_GntR_mP70_Venus= A4/B4 

	f_S70_RNAP_Venus= (actors ^abs(n_S70_RNAP_mP70)) / ( actors ^ abs(n_S70_RNAP_mP70) + abs(K_S70_RNAP_mP70) ^ abs(n_S70_RNAP_mP70)) 

	

	# num4= W_RNAP_mP70 + W_S70_RNAP_mP70*f_S70_RNAP_Venus + f_gluc_GntR* W_GntR_mP70*f_GntR_mP70_Venus #Here, both S70 AND GntR bind to mP70 promoter sequence 
    # den4= 1 + num4



	# control_array[4]= num4/den4 
	control_array[4]= (W_RNAP_mP70 + W_S70_RNAP_mP70*f_S70_RNAP_Venus)/(1+W_RNAP_mP70 + W_S70_RNAP_mP70*f_S70_RNAP_Venus+W_GntR_mP70*f_GntR_mP70_Venus) 

	#######################################

	#5 for BFP 

	#AS28 is like an inhibitor molecule for S28- like how Gluconate is to GntR, but in reverse .. as more and more as28 binds to s28, inhibition starts

	actor_set3= [protein_AS28]'  
		
	actor3= prod(actor_set3)


	#INHIBITED S28 PARTICIPATES IN THIS 

	f_AS28_S28_BFP = abs(actor3)^abs(n_AS28_S28) / ( abs(actor3)^abs(n_AS28_S28) + abs(K_AS28_S28)^abs(n_AS28_S28) )
	
	protein_S28= (1-f_AS28_S28_BFP)*protein_S28

	actor_set5= [protein_S28]'

	actor5= prod(actor_set5)



	# a= (1 + (actor3/K_AS28_S28_BFP))
	# Km=Km_AS28_S28
	# f_den= actor3_4 /(a*Km + actor3_4) 
	# f_AS28_S28_BFP= f_den 

	# protein_S28= (1-f_AS28_S28_BFP)*protein_AS28 #free S28
	
	# actor_set4= [protein_S28]'

	# actor4= prod(actor_set4)

	term_num4= real(abs(actor5)^ abs(n_S28_RNAP_P28))
	term_den4= real(abs(actor5))^abs(n_S28_RNAP_P28) + abs(K_S28_RNAP_P28)^abs(n_S28_RNAP_P28)

	f_S28_RNAP_BFP =term_num4/term_den4 


	num5= W_RNAP_P28 + W_S28_RNAP_P28 *f_S28_RNAP_BFP 

	# Factor=10

	# den5= 1 + num5 + (Factor*W_AS28_S28_P28*(1-f_den))  #Only the free fraction after being bound by as28 is considered for repression/derepresion
	den5= 1 + num5

	control_array[5]= (num5/den5 )
	
	#print("The control_array term is $(control_array[5])","\n")
	correction_term = (x[16]/100.0) 
    control_array = control_array*correction_term
	
	# return -
	return control_array
end 

#
# ----------------------------------------------------------------------------------- #
# Function: calculate_translation_control_array
# Description: Calculate the translation control array at time t
# Generated on: 2021-04-29T20:40:54.447
#
# Input arguments:
# t::Float64 => Current time value (scalar)
# x::Array{Float64,1} => State array (number_of_species x 1) (Basically its a unit vector with 9 entries of 1.0 each)
# data_dictionary::Dict{String,Any} => Dictionary holding model parameters
#
# Output arguments:
# control_array::Array{Float64,1} => Translation control array (number_of_genes x 1) at time t
# ----------------------------------------------------------------------------------- #
function calculate_translation_control_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{String,Any})

	# initialize the control -
	control_array = ones(5)

	# #manually add txtl resources decay term for now (logistic decay)
	# f(y)=(1+exp((-6*pi)/(0.6*sqrt(3))))/(1+exp(((y-6)*pi)/(0.6*sqrt(3))))
	# correction_term_manual=f(t)
	# control_array=control_array*correction_term_manual

	# estimate transcription_capacity_delay and transcription_capacity_slope parameters #these are what he added
	#translation_capacity_delay = data_dictionary["translation_capacity_delay"] 
	#translation_capacity_slope = data_dictionary["translation_capacity_slope"] 
	# build control array with these two logistic function parameters
	#f(y)= (1+exp((-translation_capacity_delay*pi)/(translation_capacity_slope*sqrt(3))))/(1+exp(((y-translation_capacity_delay)*pi)/(translation_capacity_slope*sqrt(3))))
	#correction_term_manual=f(t)
	#control_array=control_array*correction_term_manual #THIS WAS ORIGINALLY COMMENTED OUT
	
	
	correction_term = (x[16]/100.0) 
    control_array = control_array*correction_term
	# # # return -
	return control_array
end
