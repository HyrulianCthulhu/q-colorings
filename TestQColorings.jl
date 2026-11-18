# TestQColorings.jl
# Validación: compara el estimador con conteos EXACTOS conocidos.
#   Camino/árbol con n vértices:  q (q-1)^(n-1)
#   Ciclo C_n:                    (q-1)^n + (-1)^n (q-1)
#   Grafo completo K_n:           q (q-1) (q-2) ... (q-n+1)
include("QColorings.jl")
using .QColorings
using Graphs
using Random

Random.seed!(1)

casos = [
    ("Ciclo C6, q=4",         cycle_graph(6),    4, 3^6 + 3),
    ("Camino de 7 nodos, q=3", path_graph(7),     3, 3 * 2^6),
    ("Completo K4, q=6",      complete_graph(4), 6, 6 * 5 * 4 * 3),
    ("Completo K5, q=6",      complete_graph(5), 6, 6 * 5 * 4 * 3 * 2),
]

for (nombre, g, q, exacto) in casos
    est = estimate_num_colorings(g, q)
    err = abs(est - exacto) / exacto
    println(nombre, ": exacto = ", exacto,
            ", estimado = ", round(est, digits = 1),
            ", error relativo = ", round(err, digits = 3))
end
