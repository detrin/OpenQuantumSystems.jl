

function take_el_part(A::Array, a, b, vibindices)
    a1 = vibindices[a][1]
    a2 = vibindices[a][end]
    b1 = vibindices[b][1]
    b2 = vibindices[b][end]

    return A[a1:a2, b1:b2]
end