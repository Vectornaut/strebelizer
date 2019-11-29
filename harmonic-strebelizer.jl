module Strebelizer

using LinearAlgebra, Plots

export Mesh, edge, next, show, cot_to_next, distortion, example

struct Mesh{S, R}
  vertices::Vector{S}
  pos_edges::Vector{R}
  pos_next::Vector{Int}
  neg_next::Vector{Int}
  reverser
end

Base.length(m::Mesh) = length(m.pos_edges)

function edge(m::Mesh, e::Int)
  if e > 0
    m.pos_edges[e]
  else
    m.reverser * m.pos_edges[-e]
  end
end

function next(m::Mesh, e::Int)
  if e > 0
    m.pos_next[e]
  else
    m.neg_next[-e]
  end
end

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

function cot_to_next(m::Mesh, e::Int)
  ratio = -edge(m, e) / edge(m, next(m, e))
  real(ratio) / imag(ratio)
end

# assume m is Delaunay (no checks for now)
function distortion(m::Mesh)
  pos_distorts = [
    (cot_to_next(m, next(m, e)) + cot_to_next(m, next(m, -e)))/2
    for e in 1:length(m)
  ]
  Mesh(pos_distorts, m.pos_next, m.neg_next, 1)
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
