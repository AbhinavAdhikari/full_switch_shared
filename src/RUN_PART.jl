using Shell
# Shell.run("rm -f simulated/gluconate_dynamics/log-5.0mM/*");
# Shell.run("rm -f simulated/gluconate_dynamics/log-3.0mM/*");
# Shell.run("rm -f simulated/gluconate_dynamics/log-2.0mM/*");
# Shell.run("rm -f simulated/gluconate_dynamics/log-1.0mM/*");
# Shell.run("rm -f simulated/gluconate_dynamics/log0.0mM/*");
# Shell.run("rm -f simulated/gluconate_dynamics/log1.0mM/*");

# Shell.run("rm -f simulated/poets_ensemble/*");

# Estimate parameters
# include("Parameter_Estimation_W_Splined_common_promoter.jl");
# include("Updated_Driver_common_promoter.jl");


# After estimating parameters
# include("Gluconate_dynamics.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T10.dat
include("dose_response_ensemble.jl") #MAKE SURE to confirm the simulation file being used, eg. PC_T10.dat
include("plot_dose_response_venus.jl")
include("plot_dose_response_bfp.jl")
include("plot_venus_protein.jl")
include("plot_venus_mrna.jl")
include("plot_bfp_protein.jl")
include("plot_bfp_mrna.jl")

# Sensitivity
# include("sensitivity-gntr-gluconate-updated.jl")
# include("visualize-sensitivity-array.jl")

