# QColorings.jl
# Gibbs sampler para q-coloreos (Häggström, Ejemplo 7.3) y estimador del
# número de q-coloreos (Häggström, Cap. 9).
module QColorings

using Graphs

export is_proper, greedy_coloring, random_site_sweep!, estimate_num_colorings

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

# Un "barrido": k = nv(g) pasos del Gibbs sampler del Ejemplo 7.3.
# Cada paso: (1) escoge un vértice v al azar, (2) le pone un color al azar
# entre los que NO usan sus vecinos, (3) deja los demás vértices igual.
function random_site_sweep!(g::AbstractGraph, coloring::Vector{Int}, q::Int, used::BitVector)
    for _ in 1:nv(g)
        v = rand(1:nv(g))
        fill!(used, false)
        for w in neighbors(g, v)
            used[coloring[w]] = true
        end
        free = findall(==(false), used)
        coloring[v] = rand(free)
    end
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
            random_site_sweep!(h, coloring, q, used)
        end
        hits = 0
        for _ in 1:nsamples
            for _ in 1:gap
                random_site_sweep!(h, coloring, q, used)
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
