include("MarkovFunctions.jl")
include("QColorings.jl")
using .MarkovFunctions, .QColorings
using Graphs
using Plots

"""
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
"""

k = 60
a = 435

q = 31
g = Graph(k, a)

print("Starting estimation\n")
numereishon = estimate_num_colorings(g, q)
print(numereishon)
print("\nTotal num of configs\n")
denominereishon = big(q)^k
print(denominereishon)
print("\nRatio:\n")
print(numereishon/denominereishon)

print("\nDone")