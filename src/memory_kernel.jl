

function take_el_part(A::Array, a, b, vibindices)
    a1 = vibindices[a][1]
    a2 = vibindices[a][end]
    b1 = vibindices[b][1]
    b2 = vibindices[b][end]

    return A[a1:a2, b1:b2]
end

function MemoryKernel_1_traced(H_II_t::Array, H_II_s::Array, W_bath::Array, agg, FCProd, aggIndices, vibindices; groundState = false)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    if groundState
        elLen += 1
    end
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_t_s = H_II_t * H_II_s 

    if !groundState
        for a in 2:elLen+1
            for c in 2:elLen+1
                H_II_t_s_ac = take_el_part(H_II_t_s, a, c, vibindices)
                for d in 2:elLen+1
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    MK_big = H_II_t_s_ac * W_bath_cd
                    MemoryKernel[a-1, d-1, c-1, d-1] = trace_bath_part(MK_big, a, d, agg, FCProd, aggIndices, vibindices)
                end
            end
        end
    else
        for a in 1:elLen
            for c in 1:elLen
                H_II_t_s_ac = take_el_part(H_II_t_s, a, c, vibindices)
                for d in 1:elLen
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    MK_big = H_II_t_s_ac * W_bath_cd
                    MemoryKernel[a, d, c, d] = trace_bath_part(MK_big, a, d, agg, FCProd, aggIndices, vibindices)
                end
            end
        end
    end
    return MemoryKernel 
end

function MemoryKernel_2_traced(H_II_t::Array, H_II_s::Array, W_bath::Array, agg, FCProd, aggIndices, vibindices; groundState = false)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    if groundState
        elLen += 1
    end
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    if !groundState
        for a in 2:elLen+1
            for b in 2:elLen+1
                for c in 2:elLen+1
                    H_II_t_ac = take_el_part(H_II_t, a, c, vibindices)
                    for d in 2:elLen+1
                        H_II_s_db = take_el_part(H_II_s, d, b, vibindices)
                        W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                        MK_big = H_II_t_ac * W_bath_cd * H_II_s_db
                        MemoryKernel[a-1, b-1, c-1, d-1] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                    end
                end
            end
        end
    else
        for a in 1:elLen
            for b in 1:elLen
                for c in 1:elLen
                    H_II_t_ac = take_el_part(H_II_t, a, c, vibindices)
                    for d in 1:elLen
                        H_II_s_db = take_el_part(H_II_s, d, b, vibindices)
                        W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                        MK_big = H_II_t_ac * W_bath_cd * H_II_s_db
                        MemoryKernel[a, b, c, d] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                    end
                end
            end
        end
    end
    return MemoryKernel 
end

function MemoryKernel_3_traced(H_II_t::Array, H_II_s::Array, W_bath::Array, agg, FCProd, aggIndices, vibindices; groundState = false)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    if groundState
        elLen += 1
    end
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)

    if !groundState
        for a in 2:elLen+1
            for b in 2:elLen+1
                for c in 2:elLen+1
                    H_II_s_ac = take_el_part(H_II_s, a, c, vibindices)
                    for d in 2:elLen+1
                        H_II_t_db = take_el_part(H_II_t, d, b, vibindices)
                        W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                        MK_big = H_II_s_ac * W_bath_cd * H_II_t_db
                        MemoryKernel[a-1, b-1, c-1, d-1] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                    end
                end
            end
        end
    else
        for a in 1:elLen
            for b in 1:elLen
                for c in 1:elLen
                    H_II_s_ac = take_el_part(H_II_s, a, c, vibindices)
                    for d in 1:elLen
                        W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                        H_II_t_db = take_el_part(H_II_t, d, b, vibindices)
                        MK_big = H_II_s_ac * W_bath_cd * H_II_t_db
                        MemoryKernel[a, b, c, d] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                    end
                end
            end
        end
    end
    return MemoryKernel 
end

function MemoryKernel_4_traced(H_II_t::Array, H_II_s::Array, W_bath::Array, agg, FCProd, aggIndices, vibindices; groundState = false)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    if groundState
        elLen += 1
    end
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_s_t = H_II_s * H_II_t 

    if !groundState
        for a in 2:elLen+1
            for b in 2:elLen+1
                for d in 2:elLen+1
                    W_bath_ad = take_el_part(W_bath, a, d, vibindices)
                    H_II_s_t_db = take_el_part(H_II_s_t, d, b, vibindices)
                    MK_big = W_bath_ad * H_II_s_t_db
                    MemoryKernel[a-1, b-1, a-1, d-1] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                end
            end
        end
    else
        for a in 1:elLen
            for b in 1:elLen
                for d in 1:elLen
                    W_bath_ad = take_el_part(W_bath, a, d, vibindices)
                    H_II_s_t_db = take_el_part(H_II_s_t, d, b, vibindices)
                    MK_big = W_bath_ad * H_II_s_t_db
                    MemoryKernel[a, b, a, d] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                end
            end
        end
    end
    return MemoryKernel 
end

function MemoryKernel_traced(H_II_t::Array, H_II_s::Array, W_bath::Array, agg, FCProd, aggIndices, vibindices; groundState = false)
    aggIndLen = length(aggIndices)
    vibLen = length(vibindices[2])
    elLen = length(agg.molecules)
    if groundState
        elLen += 1
    end
    MemoryKernel = zeros(ComplexF64, elLen, elLen, elLen, elLen)
    H_II_s_t = H_II_s * H_II_t 
    H_II_t_s = H_II_t * H_II_s

    if !groundState
        for a in 2:elLen+1
            for c in 2:elLen+1
                H_II_t_s_ac = take_el_part(H_II_t_s, a, c, vibindices)
                H_II_s_ac = take_el_part(H_II_s, a, c, vibindices)
                H_II_t_ac = take_el_part(H_II_t, a, c, vibindices)
                for d in 2:elLen+1
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    for b in 2:elLen+1
                        H_II_s_db = take_el_part(H_II_s, d, b, vibindices)
                        
                        MK_big = - H_II_t_ac * W_bath_cd * H_II_s_db

                        H_II_t_db = take_el_part(H_II_t, d, b, vibindices)
                        MK_big[:, :] -= H_II_s_ac * W_bath_cd * H_II_t_db

                        if a == c
                            W_bath_ad = take_el_part(W_bath, a, d, vibindices)
                            H_II_s_t_db = take_el_part(H_II_s_t, d, b, vibindices)
                            MK_big[:, :] += W_bath_ad * H_II_s_t_db
                        end

                        if d == b
                            MK_big[:, :] += H_II_t_s_ac * W_bath_cd
                        end
                        MemoryKernel[a-1, b-1, c-1, d-1] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                    end
                end
            end
        end
    else
        for a in 1:elLen
            for c in 1:elLen
                H_II_t_s_ac = take_el_part(H_II_t_s, a, c, vibindices)
                H_II_s_ac = take_el_part(H_II_s, a, c, vibindices)
                H_II_t_ac = take_el_part(H_II_t, a, c, vibindices)
                for d in 1:elLen
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    W_bath_cd = take_el_part(W_bath, c, d, vibindices)
                    for b in 1:elLen
                        H_II_s_db = take_el_part(H_II_s, d, b, vibindices)
                        
                        MK_big = - H_II_t_ac * W_bath_cd * H_II_s_db

                        H_II_t_db = take_el_part(H_II_t, d, b, vibindices)
                        MK_big[:, :] -= H_II_s_ac * W_bath_cd * H_II_t_db

                        if a == c
                            W_bath_ad = take_el_part(W_bath, a, d, vibindices)
                            H_II_s_t_db = take_el_part(H_II_s_t, d, b, vibindices)
                            MK_big[:, :] += W_bath_ad * H_II_s_t_db
                        end

                        if d == b
                            MK_big[:, :] += H_II_t_s_ac * W_bath_cd
                        end
                        MemoryKernel[a, b, c, d] = trace_bath_part(MK_big, a, b, agg, FCProd, aggIndices, vibindices)
                    end
                end
            end
        end
    end
    return MemoryKernel 
end