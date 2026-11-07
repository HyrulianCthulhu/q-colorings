#MarkovFunctions.jl
module MarkovFunctions

export Initiate, Update, SimulateChain

function Initiate(x::Float64, initial::AbstractVector{Float64})
    b = initial[1]
    i = 1

    while b <= x < 1
        i += 1
        b += initial[i]
    end

    return i
end

function Update(x::Float64, current::Int64, P::AbstractMatrix{Float64})
    b = P[current, 1]
    i = 1

    while b <= x < 1
        i += 1
        b += P[current, i]
    end
    
    return i
end

function SimulateChain(initial::AbstractVector{Float64}, P::AbstractMatrix{Float64}, steps::Int64)
    # --- 1. Upfront Validations ---
    
    # Validate Initial Distribution
    if !(sum(initial) ≈ 1.0)
        throw(ArgumentError("The given initial distribution does not add to 1"))
    end

    # Validate Matrix Squareness
    if size(P)[1] != size(P)[2]
        throw(ArgumentError("The transition matrix must be square"))
    end

    # Validate All Matrix Rows
    for r in 1:rows
        if !(sum(@view P[r, :]) ≈ 1.0)
            throw(ArgumentError("Row $r in the transition matrix does not add to 1"))
        end
    end

    # --- 2. Begin the Simulation Loop ---
    
    path = Vector{Int64}(undef, steps)
    u = rand()
    current = Initiate(u, initial)
    path[1] = current

    @inbounds for i in 2:steps
        u = rand()
        new = Update(u, current, P)
        path[i] = new

        current = new
    end
        
    return path
end

end #End of module