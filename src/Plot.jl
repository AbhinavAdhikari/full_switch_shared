using DataFrames
using Plots
using CSV

#not including the 2mM or the standard deviations as of now. this file will be updated as and when time comes. 

DF1=DataFrame(CSV.File("Protein_data_BFP.csv"))
DF2=DataFrame(CSV.File("Protein_data_Venus.csv"))
a=DF1[!,"time(h)"]
b=DF1[!,"mean_0mM"]
c=DF1[!,"mean_10mM"]
plot(a,b, label="No Gluconate BFP") 
plot!(a,c, label="Gluconate 10mM- BFP") 
plot!(xlabel= "Time in h", ylabel= "RFU readings" )

d=DF2[!, "mean_0mM"]
e=DF2[!,"mean_10mM"]
plot!(a,d, label ="No Gluconate  Venus")
plot!(a,e, label= "Gluconate 10mM Venus")
savefig("Experimental_Run1.png")