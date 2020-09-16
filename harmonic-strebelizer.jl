module Strebelizer

using LinearAlgebra, Formatting

export Mesh, edge, next, cot_to_next, distortion, example

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

Base.length(m::Mesh) = length(m.pos_next)

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

function Base.show(m::Mesh)
  e_form = string("{: 0", 1 + length(digits(length(m))), "d}")
  v_form = string("{:0", length(digits(maximum([m.pos_vertex; m.neg_vertex]))), "d}")
  to_endpts = FormatExpr(string(e_form, ": ", v_form, " --> ", v_form))
  next_form = string(e_form, " --> ", e_form)
  to_next = FormatExpr(string(next_form, "   ", next_form))
  for e in 1:length(m)
    printfmtln(to_endpts, e, vertex(m, -e), vertex(m, e))
  end
  println()
  for e in 1:length(m)
    printfmtln(to_next, e, next(m, e), -e, next(m, -e))
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

function cot_to_next(flat_struct::DirEdgeFn)
  m = flat_struct.mesh
  DirEdgeFn(
    m,
    [cot_to_next(flat_struct, e) for e in 1:length(m)],
    [cot_to_next(flat_struct, -e) for e in 1:length(m)]
  )
end

# assume m is Delaunay (no checks for now)
function energy_density(flat_struct::DirEdgeFn)
  m = flat_struct.mesh
  r = cot_to_next(flat_struct)
  EdgeFn(m, [abs2(flat_struct(e)) * (r(e) + r(-e)) / 2 for e in 1:length(m)])
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

function flat_torus(a_charge, b_charge, a_res, b_res)
  N = a_res * b_res
  pos_charge = repeat([a_charge, b_charge, -(a_charge + b_charge)], N)
  
  pos_next = zeros(Int, 3*a_res*b_res)
  neg_next = zeros(Int, 3*a_res*b_res)
  for i in 0:(a_res-1), j in 0:(b_res-1)
    k = 3*(i + a_res*j)
    pos_next[k + 1] = -(k + 2)
    pos_next[k + 2] = -(k + 3)
    pos_next[k + 3] = -(k + 1)
    
    e1 = k + 3
    e2 = 3a_res*mod(j+1, b_res) + 3i + 1
    e3 = 3a_res*j + 3mod(i-1, a_res) + 2
    neg_next[e1] = e2
    neg_next[e2] = e3
    neg_next[e3] = e1
  end
  
  DirEdgeFn(Mesh(pos_next, neg_next), pos_charge, -pos_charge)
end

function dumb_example()
  m = Mesh([-2, -3, -1], [3, 1, 2])
  pos_charge = [1, im, -(1+im)]
  flat_struct = DirEdgeFn(m, pos_charge, -pos_charge)
  show(m)
  println()
  show(energy_density(flat_struct).values)
end

function example()
  flat_struct = flat_torus(1, -0.1 + im, 3, 2)
  show(flat_struct.mesh)
  println()
  show(energy_density(flat_struct).values)
end

end
