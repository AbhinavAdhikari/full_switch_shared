using LsqFit

# Define the Hill function
hill(x, p) = p[1] ./ (1 .+ (p[2] ./ x).^p[3])

# Define the x and y values
xdata = range(0.5, 2.5, length=100)
ydata = hill(xdata, [1.0, 1.0, 1.0]) .+ 0.1*randn(length(xdata))

# Define the fit model
model(x, p) = hill(x, p)

# Define the initial parameter guesses
p0 = [1.0, 1.0, 1.0]

# Fit the data
fit_result = curve_fit(model, xdata, ydata, p0)

# Print the fit result
println("Fit parameters: ", fit_result.param)

# Plot the data and the fit
using PyPlot
scatter(xdata, ydata, label="Data")
plot(xdata, model(xdata, fit_result.param), label="Fit")
savefig("hill_fit.png")
