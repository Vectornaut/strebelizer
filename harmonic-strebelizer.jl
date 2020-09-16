module Strebelizer

using LinearAlgebra, Plots

export Mesh, edge, next, show, cot_to_next, distortion, example

## to do:
## - find distance functions to the two sides of a non-separating curve
## - calculate gradient flow carefully
## - figure out how to treat separating case?

struct Mesh
  pos_next::Vector{Int}
  neg_next::Vector{Int}
  pos_vertex::Vector{Int}
  neg_vertex::Vector{Int}
end

function Mesh(pos_next::Vector{Int}, neg_next::Vector{Int})
  # initialize the mesh's next-edge table
  mesh = Mesh(pos_next, neg_next, zeros(Int, n_edges), zeros(Int, n_edges))
  
  # `newfound[e]` is true if edge `e` points to a vertex we haven't seen yet
  n_edges = length(pos_next)
  newfound = trues(n_edges)
  
  v = 1 # the next vertex label to assign
  for edge in 1:n_edges
    neighbor = edge
    do
      neighbor = mesh.next(neighbor)
    while neighbor != edge
    if newfound[e]
      pos_vertex[e] = v
      pos_vertex[pos_next[e]] = v
      pos_vertex[pos_next[pos_next[e]]
      newfound[e] = false
    end
  end
end

Base.length(m::Mesh) = length(m.pos_edges)

# swing `e` clockwise around the vertex it points toward, and return the next
# edge you hit
next(m::Mesh, e::Int) = e > 0 ? m.pos_next[e] : m.neg_next[-e]

struct VertexFn{R}
  mesh::Mesh
  values::Vector{R}
end

(f::VertexFn)(v) = f.values[v]

struct EdgeFn{R}
  mesh::Mesh
  values::Vector{R}
end

(w::EdgeFn)(e) = w.values[e]

struct DirEdgeFn{R}
  mesh::Mesh
  pos_values::Vector{R}
  neg_values::Vector{R}
end

(w::DirEdgeFn)(e) = e > 0 ? w.pos_values[e] : w.neg_values[-e]

function show(m::Mesh)
  println("reverser: ", m.reverser)
  for e in 1:length(m)
    println(e, ": ", m.pos_edges[e])
  end
  for e in 1:length(m)
    println(e, " --> ", next(m, e), "\t\t", -e, " --> ", next(m, -e))
  end
end

##function refine{mesh::Mesh}
  ##
##end

function cot_to_next(flat_struct::DirEdgeFn, e::Int)
  ratio = -flat_struct(e) / edge(flat_struct, next(m, e))
  real(ratio) / imag(ratio)
end

# assume m is Delaunay (no checks for now)
function distortion(flat_struct::DirEdgeFn)
  m = flat_struct.mesh
  pos_distorts = [
    (cot_to_next(flat_struct, next(m, e)) + cot_to_next(m, next(m, -e)))/2
    for e in 1:length(m)
  ]
  Mesh(m.vertices, pos_distorts, m.pos_next, m.neg_next, m.pos_vertex, m.neg_vertex, 1)
end

# --- examples ---


function torus_mesh(a_edge::R, b_edge::R, a_res, b_res) where R
  N = a_res * b_res
  pos_edges = repeat([a_edge, b_edge, -a_edge - b_edge], N)
  
  pos_next = zeros(Int, 3*a_res*b_res)
  neg_next = zeros(Int, 3*a_res*b_res)
  for i in 0:(a_res-1), j in 0:(b_res-1)
    k = 3*(i + a_res*j)
    pos_next[k + 1] = k + 2
    pos_next[k + 2] = k + 3
    pos_next[k + 3] = k + 1
    
    e1 = k + 3
    e2 = 3a_res*mod(j+1, b_res) + 3i + 1
    e3 = 3a_res*j + 3mod(i-1, a_res) + 2
    neg_next[e1] = -e2
    neg_next[e2] = -e3
    neg_next[e3] = -e1
  end
  
  Mesh{R}(pos_edges, pos_next, neg_next, -1)
end

function dumb_example()
  m = Mesh{Complex{Float64}}(
    [1im, -1, 1-1im],
    [2, 3, 1],
    [-2, -3, -1],
    -1
  )
  show(m)
  println()
  show(distortion(m))
end

function example()
  m = torus_mesh(complex(1.0), -0.1 + 1.0im, 2, 3)
  show(m)
  println()
  show(distortion(m))
end

end
