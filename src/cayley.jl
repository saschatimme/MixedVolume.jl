export CayleyConfiguration, CircuitTable, MixedCell, circuitdot

const CellIndices = Vector{Tuple{Int,Int}} # Maybe also SVector{N, Tuple{T, T}}

"""
    CayleyConfiguration{I<:Integer}

The Cayley matrix of a set of matrices ``A_1, …, A_n`` where ``A_i ∈ 𝐑^{n×mᵢ}``.
"""
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

"""
    linearindex(CC::CayleyConfiguration, i::Integer, col::Integer)

Get the linear index into the Caylex matrix corresponding to the configuration `i` and column `col`.
"""
function linearindex(CC::CayleyConfiguration, configuration::Integer, col::Integer)
    CC.offsets[configuration] + col
end

"""
    linearindices(CC::CayleyConfiguration, cellindices::CellIndices)

Get the linear indexes into the Caylex matrix corresponding to the cell indices.
"""
function linearindices(A::CayleyConfiguration, cellindices::CellIndices)
    linearindices = Int[]
    for i=1:length(cellindices)
        aᵢ, bᵢ = cellindices[i]
        push!(linearindices, linearindex(A, i, aᵢ), linearindex(A, i, bᵢ))
    end
    linearindices
end

Base.getindex(CC::CayleyConfiguration, i::Int) = getindex(CC.A, i)
Base.getindex(CC::CayleyConfiguration, I::Vararg{Int, 2}) = getindex(CC.A, I...)
Base.size(CC::CayleyConfiguration) = size(CC.A)
Base.size(CC::CayleyConfiguration, I) = size(CC.A, I)
