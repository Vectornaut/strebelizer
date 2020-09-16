module Strebelizer

using LinearAlgebra, Formatting

export Mesh, edge, next, cot_to_next, distortion, example

## to do:
## - find distance functions to the two sides of a non-separating curve
## - calculate gradient flow carefully
## - figure out how to treat separating case?

struct RibbonGraph
  pos_next::Vector{Int}
  neg_next::Vector{Int}
  pos_vertex::Vector{Int}
  neg_vertex::Vector{Int}
end

function RibbonGraph(pos_next::Vector{Int}, neg_next::Vector{Int})
  # initialize the mesh's next-edge table
  n_edges = length(pos_next)
  graph = RibbonGraph(pos_next, neg_next, zeros(Int, n_edges), zeros(Int, n_edges))
  
  # label the vertices
  v = 1 # the next vertex label to assign
  for edge in [1:n_edges; -n_edges:-1]
    neighbor = edge
    while (vertex(graph, neighbor) == 0)
      if neighbor > 0
        graph.pos_vertex[neighbor] = v
      else
        graph.neg_vertex[-neighbor] = v
      end
      neighbor = next(graph, neighbor)
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
  graph
end

Base.length(g::RibbonGraph) = length(g.pos_next)

# the edge located clockwise of `e` around the vertex `e` points toward
next(g::RibbonGraph, e::Int) = e > 0 ? g.pos_next[e] : g.neg_next[-e]

# the vertex `e` points toward
vertex(g::RibbonGraph, e::Int) = e > 0 ? g.pos_vertex[e] : g.neg_vertex[-e]

struct VertexFn{R}
  graph::RibbonGraph
  values::Vector{R}
end

(f::VertexFn)(v) = f.values[v]

struct EdgeFn{R}
  graph::RibbonGraph
  values::Vector{R}
end

(w::EdgeFn)(e) = w.values[e]

struct DirEdgeFn{R}
  graph::RibbonGraph
  pos_values::Vector{R}
  neg_values::Vector{R}
end

(w::DirEdgeFn)(e) = e > 0 ? w.pos_values[e] : w.neg_values[-e]

function Base.show(g::RibbonGraph)
  e_form = string("{: 0", 1 + length(digits(length(g))), "d}")
  v_form = string("{:0", length(digits(maximum([g.pos_vertex; g.neg_vertex]))), "d}")
  to_endpts = FormatExpr(string(e_form, ": ", v_form, " --> ", v_form))
  next_form = string(e_form, " --> ", e_form)
  to_next = FormatExpr(string(next_form, "   ", next_form))
  for e in 1:length(g)
    printfmtln(to_endpts, e, vertex(g, -e), vertex(g, e))
  end
  println()
  for e in 1:length(g)
    printfmtln(to_next, e, next(g, e), -e, next(g, -e))
  end
end

##function refine{mesh::Mesh}
  ##
##end

function cot_to_next(flat_struct::DirEdgeFn, e)
  neighbor = next(flat_struct.graph, e)
  ratio = -flat_struct(e) / flat_struct(neighbor)
  real(ratio) / imag(ratio)
end

function cot_to_next(flat_struct::DirEdgeFn)
  graph = flat_struct.graph
  DirEdgeFn(
    graph,
    [cot_to_next(flat_struct, e) for e in 1:length(graph)],
    [cot_to_next(flat_struct, -e) for e in 1:length(graph)]
  )
end

# assume m is Delaunay (no checks for now)
function energy_density(flat_struct::DirEdgeFn)
  g = flat_struct.graph
  r = cot_to_next(flat_struct)
  EdgeFn(g, [abs2(flat_struct(e)) * (r(e) + r(-e)) / 2 for e in 1:length(g)])
end

# --- tests ---

function test_graphs()
  triangle = RibbonGraph([-2, -3, -1], [3, 1, 2])
  triangle_tests = [
    [
      length(unique(triangle.pos_vertex)) == 3,
      length(unique(triangle.neg_vertex)) == 3
    ];
    [vertex(triangle, e) == vertex(triangle, -mod1(e+1, 3)) for e in 1:3]
  ]
  
  heavy8 = RibbonGraph([-2, -3, -1], [3, 2, 1])
  heavy8_tests = [
    [
      length(unique(heavy8.pos_vertex)) == 2,
      length(unique(heavy8.neg_vertex)) == 2
    ];
    [
      vertex(heavy8, -1) == vertex(heavy8, 3),
      vertex(heavy8, 3) == vertex(heavy8, -1),
      vertex(heavy8, 1) == vertex(heavy8, -2),
      vertex(heavy8, -2) == vertex(heavy8, 2),
      vertex(heavy8, 2) == vertex(heavy8, -3)
    ]
  ]
  
  [triangle_tests, heavy8_tests]
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
  
  DirEdgeFn(RibbonGraph(pos_next, neg_next), pos_charge, -pos_charge)
end

function dumb_example()
  # not actually a flat structure, because the triangle 2, 4, -6 violates the
  # charge sum rule
  g = RibbonGraph([-2, -3, -1, 5, 6, 4], [-5, -6, -4, 2, 3, 1])
  pos_charge = [1+im, -1, -im, -1-im, 1, -im]
  flat_struct = DirEdgeFn(g, pos_charge, -pos_charge)
  
  show(g)
  println()
  show(energy_density(flat_struct).values)
end

function example()
  flat_struct = flat_torus(1, -0.1 + im, 3, 2)
  show(flat_struct.graph)
  println()
  show(energy_density(flat_struct).values)
end

end
