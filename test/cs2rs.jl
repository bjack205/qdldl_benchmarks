using SparseArrays
using LinearAlgebra

function csc2csr(A)
    n = size(A,1)
    ia = zero(A.colptr)
    ja = zero(A.rowval)
    nzval = zero(A.nzval)
    cnt = 1
    ia[1] = 1
    for r = 1:n
        # search each column for the current row
        for c = 1:n
            istart = A.colptr[c]
            istop = A.colptr[c+1] - 1
            for i = istart:istop
                if A.rowval[i] == r
                    println("i = $i ($istart, $istop), r = $r, c = $c, cnt = $cnt, x = $(A.nzval[i])")
                    nzval[cnt] = A.nzval[i]
                    ja[cnt] = c
                    cnt += 1
                end
            end
        end
        ia[r+1] = cnt
    end
    return ia, ja, nzval
end

n = 5 
# A = triu(sprandn(n,n,0.5) + I)
A.rowval
for v in nzval
    print(round(v, digits=4), ", ")
end
println(A.nzval')
ia,ja,nzval = csc2csr(A)
L = SparseMatrixCSC(n, n, ia, ja, nzval)
A == L'