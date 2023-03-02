pathify(eds::Vector{E}) where E<:AbstractEdge = [src.(eds)..., dst(eds[end])]
