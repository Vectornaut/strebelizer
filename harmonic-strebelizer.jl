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
  n_edges = length(pos_next)
  mesh = Mesh(pos_next, neg_next, zeros(Int, n_edges), zeros(Int, n_edges))
  
  # label the vertices
  v = 1 # the next vertex label to assign
  for edge in [1:n_edges; -n_edges:-1]
    neighbor = edge
    while (vertex(mesh, neighbor) == 0)
      if neighbor > 0
        mesh.pos_vertex[neighbor] = v
      else
        mesh.neg_vertex[-neighbor] = v
      end
      neighbor = next(mesh, neighbor)
      if neighbor == edge
        #  we're back where we started, and we're about to leave the loop
        v += 1
      end
    end
    if neighbor != edge
      error("Edge permutation is inconsistent")
    end
  end
  
  # return the finished mesh
  mesh
end

Base.length(m::Mesh) = length(m.pos_edges)

# the edge located clockwise of `e` around the vertex `e` points toward
next(m::Mesh, e::Int) = e > 0 ? m.pos_next[e] : m.neg_next[-e]

# the vertex `e` points toward
vertex(m::Mesh, e::Int) = e > 0 ? m.pos_vertex[e] : m.neg_vertex[-e]

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

function cot_to_next(flat_struct::DirEdgeFn, e)
  neighbor = next(flat_struct.mesh, e)
  ratio = -flat_struct(e) / flat_struct(neighbor)
  real(ratio) / imag(ratio)
end

# assume m is Delaunay (no checks for now)
function distortion(flat_struct::DirEdgeFn)
  m = flat_struct.mesh
  pos_distorts = [
    (
      cot_to_next(flat_struct, next(m, e))
      + cot_to_next(flat_struct, next(m, -e))
    ) / 2
    for e in 1:length(m)
  ]
  Mesh(m.vertices, pos_distorts, m.pos_next, m.neg_next, m.pos_vertex, m.neg_vertex, 1)
end

# --- tests ---

function test_mesh()
  tri_mesh = Mesh([-2, -3, -1], [3, 1, 2])
  tri_tests = [
    [
      length(unique(tri_mesh.pos_vertex)) == 3,
      length(unique(tri_mesh.neg_vertex)) == 3
    ];
    [vertex(tri_mesh, e) == vertex(tri_mesh, -mod1(e+1, 3)) for e in 1:3]
  ]
  
  key_mesh = Mesh([-2, -3, -1], [3, 2, 1])
  key_tests = [
    [
      length(unique(key_mesh.pos_vertex)) == 2,
      length(unique(key_mesh.neg_vertex)) == 2
    ];
    [
      vertex(key_mesh, -1) == vertex(key_mesh, 3),
      vertex(key_mesh, 3) == vertex(key_mesh, -1),
      vertex(key_mesh, 1) == vertex(key_mesh, -2),
      vertex(key_mesh, -2) == vertex(key_mesh, 2),
      vertex(key_mesh, 2) == vertex(key_mesh, -3)
    ]
  ]
  
  [tri_tests, key_tests]
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
