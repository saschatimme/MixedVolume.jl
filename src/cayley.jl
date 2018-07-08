export CayleyConfiguration, CircuitTable, MixedCell, circuitdot
# package code goes here

struct CayleyConfiguration{I<:Integer} <: AbstractMatrix{I}
    A::Matrix{I}
    offsets::Vector{Int} # Maybe also SVector{N, Int}
end

function CayleyConfiguration(Aᵢs::Matrix{I}...) where {I<:Integer}
    n = size(Aᵢs[1], 1)
    @assert all(Aᵢ -> size(Aᵢ, 1) == n, Aᵢs) "Matrices do not have the same number of rows."

    A = fill(zero(I), 2n, sum(Aᵢ -> size(Aᵢ, 2), Aᵢs))

    offsets = Int[]
    offset = 0
    for (i, Aᵢ) in enumerate(Aᵢs)
        push!(offsets, offset)
        mᵢ = size(Aᵢ, 2)

        for j=1:mᵢ
            jj = offset + j
            for k=1:n
                A[k, jj] = Aᵢ[k, j]
            end
            A[n+i, jj] = one(I)
        end

        offset += mᵢ
    end
    CayleyConfiguration(A, offsets)
end

function linearindex(CC::CayleyConfiguration, configuration::Integer, col::Integer)
    CC.offsets[configuration] + col
end

Base.getindex(CC::CayleyConfiguration, i::Int) = getindex(CC.A, i)
Base.getindex(CC::CayleyConfiguration, I::Vararg{Int, 2}) = getindex(CC.A, I...)
Base.size(CC::CayleyConfiguration) = size(CC.A)
Base.size(CC::CayleyConfiguration, I) = size(CC.A, I)

const CellIndices = Vector{Tuple{Int,Int}} # Maybe also SVector{N, Tuple{T, T}}

struct CircuitTable{I<:Integer}
    # circuits is an 2n × m Matrix
    cellcircuits::Matrix{I}
    cellcircuit_indices::Vector{I} # length m
    # The entries of the columns *not* in the mixed cell, i.e.
    # for which we have an inequality.
    entries::BitVector
    γ::I # TODO: See circuitdot to store this in cellcircuits

end

function CircuitTable(A::CayleyConfiguration{I}, cell_indices::CellIndices) where {I<:Integer}
    # construct a circuit table freshly
    linear_cell_indices = Int[]
    for i=1:length(cell_indices)
        aᵢ, bᵢ = cell_indices[i]
        push!(linear_cell_indices, linearindex(A, i, aᵢ), linearindex(A, i, bᵢ))
    end

    n2, m = size(A)

    cellcircuits = fill(zero(I), n2, m) # 2n × m

    # construct 2n × 2n submatrix indexed by the cell columns
    D = A[:, linear_cell_indices]

    # For each column γ **not** in the cell we need to compute a circuit c
    # which represents a facet of the MixedCell cone. This is a generator
    # of the nullspace of [D A[:,γ]]. We compute this by
    # [I D^-1 A[:,γ]] = 0 and scale the result such that the first entry is det(D)
    LU = lufact(D)
    D_inv = inv(LU)
    γ = round(I, det(LU))

    k = 1
    entries = trues(m)
    for i=1:m
        if k ≤ n2 && i == linear_cell_indices[k]
            k += 1
            entries[i] = false
            continue
        end
        c = D_inv * A[:,i]
        scale!(c, abs(γ))
        for j=1:n2
            cellcircuits[j, i] = round(I, c[j])
        end
    end

    CircuitTable(cellcircuits, linear_cell_indices, entries, γ)
end



struct MixedCell{I<:Integer}
    indices::CellIndices
    table::CircuitTable{I}
end

function MixedCell(A::CayleyConfiguration, cell_indices::CellIndices)
    MixedCell(cell_indices, CircuitTable(A, cell_indices))
end

"""
    circuitdot(circuittable, j, ω)
    circuitdot(mixedcell, j, ω)

Compute the dot product of the `j`-th circuit with the vector `ω`.
"""
function circuitdot(table::CircuitTable{I}, j, ω::AbstractVector{I}) where {I}
    # TODO: This branch could be quite expensive, maybe store γ in an extra column??
    if !table.entries[j]
        return zero(I)
    end
    @inbounds sum = ω[j] * table.γ
    for k = 1:length(table.cellcircuit_indices)
        @inbounds i = table.cellcircuit_indices[k]
        @inbounds sum += table.cellcircuits[k,j] * ω[i]
    end
    sum
end
circuitdot(M::MixedCell, j, ω::AbstractVector) = circuitdot(M.table, j, ω)
