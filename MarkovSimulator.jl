include("MarkovFunctions.jl")
using .MarkovFunctions
using Plots

order = 10000
steps = 100000

#Create a random stochastic matrix of the given order 
A = rand(order, order)
A = A ./ sum(A, dims=2)

#Take a row as the initial distribution
initial_A = A[1, :]

#Create the path
path = SimulateChain(initial_A, A, steps)

print("Finished calculating\n")

#Create the plot for the relative frequency of state 1
plot(1:steps, cumsum(==(1).(path[1:steps])) ./ (1:steps), ylims=(0, 1))