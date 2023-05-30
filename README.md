This is a Systems Biology ML model which predicts the response of a genetic full switch circuit to various gluconate concentrations . This mechanistic model predicts the transcriptional and translational dynamics of multiple genes in a cell free system. 

**WORK UNDER PROGRESS...**


You can download this repository as a zip file, or clone or pull it by using the command (from the command-line):


$ git pull https://github.com/


or you can clone it from :

$ git clone https://github.com/

User is required to install Julia. This can be done at: https://julialang.org/. Necessary packages will be installed when user runs the file Include.jl from the terminal the first time. 


The science behind this work can be found with a similar work: https://www.biorxiv.org/content/10.1101/2023.01.10.523462v1.abstract



|File|Description|	
| --- | --- |
|`RUN_ALL.jl`|Script to automate everything. Saves solutions in `simulated` directory.|
|`Parameter_Estimation_W_Splined_common_promoter.jl`|Solves the model equations for the ensemble of parameters sets for the test case, 0.00001mM gluconate. Saves solutions in the `simulated/poets_ensemble` directory.|
|`Updated_Driver_common_promoter.jl`|Gives a fit of the model's predicted solution versus experimental data obtained.|



