export MixedCell

struct MixedCell{I<:Integer}
    indices::CellIndices
    table::CircuitTable{I}
end

function MixedCell(A::CayleyConfiguration, cell_indices::CellIndices)
    MixedCell(cell_indices, CircuitTable(A, cell_indices))
end
