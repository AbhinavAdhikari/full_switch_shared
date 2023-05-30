using Shell

# DELETE GLUCONATE DYNAMICS SIMULATION FILES 
Shell.run("rm -f simulated/gluconate_dynamics/log-5.0mM/*");
Shell.run("rm -f simulated/gluconate_dynamics/log-3.0mM/*");
Shell.run("rm -f simulated/gluconate_dynamics/log-2.0mM/*");
Shell.run("rm -f simulated/gluconate_dynamics/log-1.0mM/*");
Shell.run("rm -f simulated/gluconate_dynamics/log0.0mM/*");
Shell.run("rm -f simulated/gluconate_dynamics/log1.0mM/*");

# DELETE POETS ENSEMBLE FILES
Shell.run("rm -f simulated/poets_ensemble/*");

# ESTIMATE PARAMETERS
include("Parameter_Estimation_W_Splined_common_promoter.jl");
include("Updated_Driver_common_promoter.jl");


# AFTER ESTIMATING PARAMETERS
include("Gluconate_dynamics.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T10.dat
include("dose_response_ensemble.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T10.dat
include("plot_dose_response_venus.jl")
include("plot_dose_response_bfp.jl")
include("plot_venus_protein.jl")
include("plot_venus_mrna.jl")
include("plot_bfp_protein.jl")
include("plot_bfp_mrna.jl")

# SENSITIVITY
include("sensitivity-gntr-gluconate-updated.jl")
include("plot_sensitivity_array.jl")

