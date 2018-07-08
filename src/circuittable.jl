export CircuitTable

"""
    CircuitTable{I<:Integer}

Holding all circuits describing a mixed cell cone ``C_M`` of mixed cell ``M``.
"""
struct CircuitTable{I<:Integer}
    # circuits is an (2n + 1) × m Matrix
    # The first 2n rows are the entries corresponding to the columns of the mixed cell
    # the last row corresponds to the index of the choosen extra column.
    circuits::Matrix{I}
    # linear indices of the mixed cell which corresponds to this circuit table
    indices::Vector{I} # length 2n
    # The entries of the columns *not* in the mixed cell, i.e.
    # for which we have an inequality.
    entries::BitVector
end

function CircuitTable(A::CayleyConfiguration{I}, cellindices::CellIndices) where {I<:Integer}
    # construct a circuit table freshly
    indices = linearindices(A, cellindices)

    n2, m = size(A)

    circuits = fill(zero(I), n2+1, m) # (2n + 1) × m

    # construct 2n × 2n submatrix indexed by the cell columns
    D = A[:, indices]

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
        if k ≤ n2 && i == indices[k]
            k += 1
            entries[i] = false
            continue
        end
        c = D_inv * A[:,i]
        scale!(c, abs(γ))
        for j=1:n2
            circuits[j, i] = round(I, c[j])
        end
        circuits[n2+1, i] = γ
    end

    CircuitTable(circuits, indices, entries)
end

Base.size(T::CircuitTable) = size(T.circuits)
Base.size(T::CircuitTable, I) = size(T.circuits, I)

"""
    dot(circuittable, j, ω)

Compute the dot product of the `j`-th circuit with the vector `ω`.
"""
function Base.LinAlg.dot(table::CircuitTable{I}, j, ω::AbstractVector{I}) where {I}
    N = size(table, 1)
    sum = zero(I)
    for k = 1:N-1
        @inbounds sum = table.circuits[k,j] * ω[table.indices[k]]
    end
    sum + table.circuits[N, j] * ω[j]
end
