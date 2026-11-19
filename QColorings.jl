# QColorings.jl
# Gibbs sampler para q-coloreos (Häggström, Ejemplo 7.3) y estimador del
# número de q-coloreos (Häggström, Cap. 9).
module QColorings

using Graphs

export is_proper, greedy_coloring, random_site_steps!, systematic_sweep!, estimate_num_colorings

# ¿Es un q-coloreo válido? (ninguna arista une dos vértices del mismo color)
function is_proper(g::AbstractGraph, coloring::Vector{Int})
    for e in edges(g)
        if coloring[src(e)] == coloring[dst(e)]
            return false
        end
    end
    return true
end

# Un q-coloreo válido inicial: recorre los vértices y a cada uno le da el
# primer color que no usan sus vecinos ya coloreados.
# Siempre funciona si q >= (máximo número de vecinos) + 1.
function greedy_coloring(g::AbstractGraph, q::Int)
    coloring = zeros(Int, nv(g))
    used = falses(q)
    for v in vertices(g)
        fill!(used, false)
        for w in neighbors(g, v)
            if coloring[w] != 0
                used[coloring[w]] = true
            end
        end
        c = findfirst(==(false), used)
        if c === nothing
            throw(ArgumentError("q = $q es muy pequeño para colorear el grafo de forma voraz (se necesita q >= grado máximo + 1)"))
        end
        coloring[v] = c
    end
    return coloring
end

#Takes the specified number of steps on the random site Gibbs sampler
function random_site_steps!(g::AbstractGraph, coloring::Vector{Int}, used::BitVector, steps::Int)
    for i in 1:steps
        #Choose a random vertex v
        v = rand(1:nv(g))
        fill!(used, false)

        #Check for used colors among neighbors, then color v with a random free one
        for w in neighbors(g, v)
            used[coloring[w]] = true
        end
        free = findall(==(false), used)
        coloring[v] = rand(free)
    end

    #Return mutated coloring, it only changed for vertex v
    return coloring
end

#One sweep (nv(g) steps) for the systematic Gibbs sampler
function systematic_sweep!(g::AbstractGraph, coloring::Vector{Int}, used::BitVector)
    for i in 1:nv(g)
        fill!(used, false)
        
        #Check for used colors among neighbors, then color vertex i with a random free one
        for w in neighbors(g, i)
            used[coloring[w]] = true
        end
        free = findall(==(false), used)
        coloring[i] = rand(free)
    end

    #Return mutated coloring
    return coloring
end

# Estima Z = número de q-coloreos de g, con el producto telescópico del libro:
#   Z = q^k * (Z_1/Z_0) * (Z_2/Z_1) * ... * (Z_m/Z_{m-1})
# donde G_j es g con solo las primeras j aristas. Cada razón Z_j/Z_{j-1} es la
# probabilidad de que, en un q-coloreo aleatorio de G_{j-1}, los dos extremos de
# la arista j tengan colores distintos; se estima con el Gibbs sampler.
#   burn_in  = barridos iniciales que se descartan (para "llegar al equilibrio")
#   nsamples = cuántos coloreos se muestrean por cada arista
#   gap      = barridos entre una muestra y la siguiente
function estimate_num_colorings(g::AbstractGraph, q::Int;
                                burn_in::Int = 30, nsamples::Int = 2000, gap::Int = 2)
    k = nv(g)
    start = greedy_coloring(g, q)     # coloreo válido de g (y de todo subgrafo)
    h = SimpleGraph(k)                # G_0: los mismos vértices, sin aristas
    used = falses(q)
    logZ = k * log(q)                 # log(Z_0) = log(q^k)

    for e in edges(g)
        x, y = src(e), dst(e)
        coloring = copy(start)
        for _ in 1:burn_in
            systematic_sweep!(h, coloring, used)
        end
        
        hits = 0
        for _ in 1:nsamples
            for _ in 1:gap
                systematic_sweep!(h, coloring, used)
            end
            if coloring[x] != coloring[y]
                hits += 1
            end
        end
        logZ += log(hits / nsamples)
        add_edge!(h, x, y)            # pasa de G_{j-1} a G_j
    end

    return exp(logZ)
end

end # module
