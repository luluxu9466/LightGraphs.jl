# Parts of this code were taken / derived from Graphs.jl. See LICENSE for
# licensing details.

immutable DijkstraHEntry
    vertex::Int
    dist
end

< (e1::DijkstraHEntry, e2::DijkstraHEntry) = e1.dist < e2.dist

abstract AbstractDijkstraVisitor

# invoked when a new vertex is first encountered
discover_vertex!(visitor::AbstractDijkstraVisitor, u, v, d) = nothing

# invoked when the distance of a vertex is determined
# (for each source vertex at the beginning, and when a vertex is popped from the heap)
# returns whether the algorithm should continue
include_vertex!(visitor::AbstractDijkstraVisitor, u, v, d) = true

# invoked when the distance to a vertex is updated (decreased)
update_vertex!(visitor::AbstractDijkstraVisitor, u, v, d) = nothing

# invoked when all neighbors of a vertex has been examined
close_vertex!(visitor::AbstractDijkstraVisitor, v) = nothing


# trivial visitor

type TrivialDijkstraVisitor <: AbstractDijkstraVisitor
end



abstract AbstractDijkstraState<:AbstractPathState
###################################################################
#
#   dijkstra_predecessor_and_distance functions
#
###################################################################

# standard dijkstra state - AbstractPathState is defined in core
type DijkstraState{T}<: AbstractDijkstraState
    parents::Vector{Int}
    dists::Vector{T}
    colormap::Vector{Int}
    heap::MutableBinaryHeap{DijkstraHEntry,DataStructures.LessThan}
    hmap::Vector{Int}
end

function create_dijkstra_state(g::AbstractGraph)
    n = nv(g)
    parents = zeros(Int, n)
    dists = fill(-1, n)
    colormap = zeros(Int, n)
    heap = mutable_binary_minheap(DijkstraHEntry)
    hmap = zeros(Int, n)
    DijkstraState(parents, dists, colormap, heap, hmap)
end

function set_source!{T}(state::DijkstraState{T}, g::AbstractGraph, s::Int)
    state.parents[s] = 0        # we are setting the parent of source to 0
    state.dists[s] = zero(T)
    state.colormap[s] = 2
end

function process_neighbors!{T}(
    state::DijkstraState,
    graph::AbstractGraph,
    edge_dists::AbstractArray{T,2},
    u::Int, du::T, visitor::AbstractDijkstraVisitor
)

    dists::Vector = state.dists
    parents::Vector{Int} = state.parents
    colormap::Vector{Int} = state.colormap
    heap::MutableBinaryHeap{DijkstraHEntry,DataStructures.LessThan} = state.heap
    hmap::Vector{Int} = state.hmap
    dv = zero(T)


    for e in out_edges(graph, u)
        v::Int = dst(e)
        v_color::Int = colormap[v]

        if v_color == 0

            edist = edge_dists[u, v]

            dists[v] = dv = du + edist
            parents[v] = u
            colormap[v] = 1
            discover_vertex!(visitor, u, v, dv)
            # push new vertex to the heap
            hmap[v] = push!(heap, DijkstraHEntry(v, dv))

        elseif v_color == 1
            dv = du + edge_dists[u, v]

            if dv < dists[v]
                dists[v] = dv
                parents[v] = u
                # update the value on the heap
                update_vertex!(visitor, u, v, dv)
                update!(heap, hmap[v], DijkstraHEntry(v, dv))
            end
        end
    end
end

function dijkstra_shortest_paths!{T}(
    graph::AbstractGraph,                # the graph
    sources::AbstractVector{Int},             # the sources
    edge_dists::AbstractArray{T,2},
    visitor::AbstractDijkstraVisitor,       # visitor object
    state::DijkstraState      # the state
    )

    # get state fields

    parents::Vector{Int} = state.parents
    dists::Vector = state.dists
    colormap::Vector{Int} = state.colormap
    heap::MutableBinaryHeap{DijkstraHEntry,DataStructures.LessThan} = state.heap
    hmap::Vector{Int} = state.hmap

    # initialize for sources

    d0 = zero(T)

    for s in sources
        set_source!(state, graph, s)
        if !include_vertex!(visitor, s, s, d0)
            return state
        end
    end

    # process direct neighbors of all sources

    for s in sources
        process_neighbors!(state, graph, edge_dists, s, d0, visitor)
        close_vertex!(visitor, s)
    end

    # main loop

    while !isempty(heap)
        # pick next vertex to include
        entry = pop!(heap)
        u::Int = entry.vertex
        du = entry.dist
        colormap[u] = 2
        if !include_vertex!(visitor, parents[u], u, du)
            return state
        end

        # process u's neighbors

        process_neighbors!(state, graph, edge_dists, u, du, visitor)
        close_vertex!(visitor, u)
    end

    state
end

function dijkstra_shortest_paths(
    graph::AbstractGraph,                # the graph
    sources::AbstractVector{Int},
    edge_dists::AbstractArray = DefaultDistance(),
    visitor::AbstractDijkstraVisitor=TrivialDijkstraVisitor()
)
    dt = eltype(edge_dists)
    state::DijkstraState = create_dijkstra_state(graph)
    dijkstra_shortest_paths!(graph, sources, edge_dists, visitor, state)
end

dijkstra_shortest_paths(graph::AbstractGraph, s::Int, edge_dists::AbstractArray) =
    dijkstra_shortest_paths(graph, [s;], edge_dists)

dijkstra_shortest_paths(graph::AbstractGraph, s::Int) =
    dijkstra_shortest_paths(graph, [s;])

















# This DijkstraState tracks predecessors and path counts
type DijkstraStateWithPred{T}<: AbstractDijkstraState
    parents::Vector{Int}
    dists::Vector{T}
    colormap::Vector{Int}
    pathcounts::Vector{Int}
    predecessors::Vector{Vector{Int}}
    heap::MutableBinaryHeap{DijkstraHEntry,DataStructures.LessThan}
    hmap::Vector{Int}
end


# Create Dijkstra state that tracks predecessors and path counts
function create_dijkstra_state_with_pred(g::AbstractGraph, dt::DataType)
    n = nv(g)
    parents = zeros(Int, n)
    dists = fill(typemax(dt), n)
    colormap = zeros(Int, n)
    pathcounts = zeros(Int, n)
    predecessors = Array(Vector{Int}, n)
    heap = mutable_binary_minheap(DijkstraHEntry)
    hmap = zeros(Int, n)

    for i = 1:n
        predecessors[i] = []
    end
    DijkstraStateWithPred(parents, dists, colormap, pathcounts, predecessors, heap, hmap)
end

function set_source_with_pred!(state::DijkstraStateWithPred, g::AbstractGraph, s::Int)
    state.parents[s] = 0        # we are setting the parent of source to 0
    state.dists[s] = 0.0
    state.colormap[s] = 2
    state.pathcounts[s] = 1
    state.predecessors[s] = []
end

function process_neighbors_with_pred!(
    state::DijkstraStateWithPred,
    graph::AbstractGraph,
    edge_dists,
    u::Int, du::Float64, visitor::AbstractDijkstraVisitor)

    dists::Vector{Float64} = state.dists
    parents::Vector{Int} = state.parents
    colormap::Vector{Int} = state.colormap
    pathcounts::Vector{Int} = state.pathcounts              # the # of paths from src to the vertex
    predecessors::Vector{Vector{Int}} = state.predecessors    # the vertex's predecessors
    heap::MutableBinaryHeap{DijkstraHEntry,DataStructures.LessThan} = state.heap
    hmap::Vector{Int} = state.hmap
    dv::Float64 = zero(Float64)

    # has_distances in distance.jl
    # use_dists = has_distances(edge_dists)

    for v in fadj(graph, u)
        v_color::Int = colormap[v]

        edist = edge_dists[u, v]
        if edist == 0.0
            edist = 1.0
        end
        if v_color == 0
            dists[v] = dv = du + edist
            parents[v] = u
            colormap[v] = 1
            discover_vertex!(visitor, u, v, dv)
            # increment pathcounts and add to predecessors
            # this ensures that changed pathcounts propagate
            pathcounts[v] += pathcounts[u]
            predecessors[v] = [u]


            # push new vertex to the heap
            hmap[v] = push!(heap, DijkstraHEntry(v, dv))

        elseif v_color == 1
            dv = du + edist
            if dv < dists[v]
                dists[v] = dv
                parents[v] = u
                # update the value on the heap
                update_vertex!(visitor, u, v, dv)
                update!(heap, hmap[v], DijkstraHEntry(v, dv))
            elseif (dv == dists[v])
                # increment pathcounts
                pathcounts[v] += pathcounts[u]
                push!(predecessors[v], u)
            end
        end
    end
end

function dijkstra_predecessor_and_distance!(
    graph::AbstractGraph,                # the graph
    edge_dists,
    sources::AbstractVector{Int},             # the sources
    visitor::AbstractDijkstraVisitor,       # visitor object
    state::DijkstraStateWithPred      # the state
)

    # get state fields

    parents::Vector{Int} = state.parents
    dists::Vector{Float64} = state.dists
    colormap::Vector{Int} = state.colormap
    heap::MutableBinaryHeap{DijkstraHEntry,DataStructures.LessThan} = state.heap
    hmap::Vector{Int} = state.hmap

    # initialize for sources

    d0 = zero(Float64)

    for s in sources
        set_source_with_pred!(state, graph, s)
        if !include_vertex!(visitor, s, s, d0)
            return state
        end
    end

    # process direct neighbors of all sources

    for s in sources
        process_neighbors_with_pred!(state, graph, edge_dists, s, d0, visitor)
        close_vertex!(visitor, s)
    end

    # main loop

    while !isempty(heap)
        # pick next vertex to include
        entry = pop!(heap)
        u::Int = entry.vertex
        du::Float64 = entry.dist
        colormap[u] = 2
        if !include_vertex!(visitor, parents[u], u, du)
            return state
        end

        # process u's neighbors

        process_neighbors_with_pred!(state, graph, edge_dists, u, du, visitor)
        close_vertex!(visitor, u)
    end

    state
end

function dijkstra_predecessor_and_distance(
    graph::AbstractGraph,
    sources::AbstractVector{Int},
    edge_dists = DefaultDistance(),
    visitor::AbstractDijkstraVisitor = TrivialDijkstraVisitor()
)
    state::DijkstraStateWithPred = create_dijkstra_state_with_pred(graph),
    dijkstra_predecessor_and_distance!(graph, edge_dists, sources, visitor, state)
end

dijkstra_predecessor_and_distance(graph::AbstractGraph, s::Int, edge_dists::AbstractArray) =
    dijkstra_predecessor_and_distance(graph, [s;], edge_dists)

dijkstra_predecessor_and_distance(graph::AbstractGraph, s::Int) =
    dijkstra_predecessor_and_distance(graph, [s;])

# function dijkstra_predecessor_and_distance(
#     graph::AbstractGraph,
#     s::Int;
#     visitor::AbstractDijkstraVisitor=TrivialDijkstraVisitor()
# )
#     state = create_dijkstra_state_with_pred(graph)
#     dijkstra_predecessor_and_distance!(graph, DefaultDistance(), [s], visitor, state)
# end
