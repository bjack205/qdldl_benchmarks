using SparseArrays

struct LQRData{n,m,T}
    Q::Vector{Diagonal{T,Vector{T}}}
    R::Vector{Diagonal{T,Vector{T}}}
    q::Vector{Vector{T}}
    r::Vector{Vector{T}}
    c::Vector{T}
    A::Vector{Matrix{T}}
    B::Vector{Matrix{T}}
    C::Vector{Matrix{T}}
    D::Vector{Matrix{T}}
    d::Vector{Vector{T}}
    x0::Vector{T}
end

function Base.rand(::Type{<:LQRData{n,m}}, N::Integer; cond=1.0, implicit=false) where {n,m}
    Nx,Nu = n,m
    Q = [Diagonal(rand(Nx)) * 10^cond for k = 1:N]
    R = [Diagonal(rand(Nu)) for k = 1:N-1]
    q = [randn(Nx) for k = 1:N]
    r = [randn(Nu) for k = 1:N-1]
    A = [zeros(Nx,Nx) for k = 1:N-1]
    B = [zeros(Nx,Nu) for k = 1:N-1]
    C = [zeros(Nx,Nx) for k = 1:N-1]
    D = [zeros(Nx,Nu) for k = 1:N-1]
    for k = 1:N-1
        Ak, Bk = gencontrollable(n, m)
        A[k] .= Ak
        B[k] .= Bk
        if implicit
            Ck,Dk = gencontrollable(n ,m)
            C[k] .= Ck
            D[k] .= Dk
        else
            C[k] .= -I(n)
        end
    end
    d = [randn(Nx) for k = 1:N-1]
    x0 = randn(Nx)
    c = randn(N) * 10
    data = LQRData{n,m,Float64}(Q,R,q,r,c,A,B,C,D,d,x0)
end
Base.size(data::LQRData{n,m}) where {n,m} = (n,m, length(data.Q))

function Base.copy(data::LQRData)
    LQRData(
        deepcopy(data.Q),
        deepcopy(data.R),
        deepcopy(data.q),
        deepcopy(data.r),
        deepcopy(data.c),
        deepcopy(data.A),
        deepcopy(data.B),
        deepcopy(data.C),
        deepcopy(data.D),
        deepcopy(data.d),
        deepcopy(data.x0)
    )
end

function build_block_diagonal(blocks)
    n = 0
    m = 0
    for block in blocks
        n += size(block, 1)
        m += size(block, 2)
    end
    A = spzeros(n, m)
    off1 = 0
    off2 = 0
    for block in blocks
        inds1 = off1 .+ (1:size(block, 1))
        inds2 = off2 .+ (1:size(block, 2))
        A[inds1, inds2] .= block
        off1 += size(block, 1)
        off2 += size(block, 2)
    end
    return A
end

function stack_vectors(vectors)
    n = 0
    for vec in vectors 
        n += size(vec, 1)
    end
    b = zeros(eltype(vectors[1]), n)
    off = 0
    for vec in vectors 
        inds = off .+ (1:size(vec, 1))
        b[inds] .= vec 
        off += size(vec, 1)
    end
    return b
end

function build_Ab(data::LQRData, form; kwargs...)
    if form == "banded"
        build_Ab_banded(data; kwargs...)
    elseif form == "kkt"
        build_Ab_kkt(data; kwargs...)
    end
end

function build_Ab_banded(data::LQRData{n,m}; remove_x1::Bool=false, reg=0.0) where {n,m}
    N = length(data.Q)
    Q,R,q,r = data.Q, data.R, data.q, data.r
    A,B,d   = data.A, data.B, data.d
    C,D     = data.C, data.D
   
    Ds = [[
            Q[k] zeros(n,m) A[k]';
            zeros(m,n) R[k] B[k]';
            A[k] B[k] -I(n)*reg 
        ] for k = 1:N-1
    ]
    push!(Ds, Q[N])

    Is = [
        [
            zeros(n,n) C[k] D[k];
            C[k]' zeros(n,n) zeros(n,m);
            D[k]' zeros(m,n) zeros(m,m);
        ] for k = 1:N-1
    ]

    b = map(1:N) do k
        dk = k == 1 ? data.x0 : d[k-1]
        if k == N
            [dk; q[k]]
        else
            [dk; q[k]; r[k]]
        end
    end

    if remove_x1
        Is[1] = zeros(m,m)
        Ds[1] = Ds[1][n+1:end, n+1:end]
        b[1] = r[1]
        b[2] = [A[1]*data.x0 + d[1]; q[2]; r[2]]
    else
        pushfirst!(Ds, -I(n)*reg)
    end
    push!(Is, Is[end][1:2n,1:2n])

    Ds = build_block_diagonal(Ds)
    Is = build_block_diagonal(Is)
    b = Vector(stack_vectors(b))
    A = Ds + Is
    return A,-b
end

function build_Ab_kkt(data::LQRData{Nx,Nu}; reg=1e-8) where {Nx,Nu}
    N = length(data.Q)
    Q,R,q,r = data.Q, data.R, data.q, data.r
    A,B,d   = data.A, data.B, data.d
    C,D     = data.C, data.D
    Np = N*Nx + (N-1)*Nu
    Nd = N*Nx

    H = blockdiag(
        [blockdiag(sparse(Qk),sparse(Rk)) for (Qk,Rk) in zip(Q,R)]..., sparse(Q[end])
    )
    D1 = [
        sparse(-I, Nx, Np);
        blockdiag([sparse([A[k] B[k]]) for k = 1:N-1]...) spzeros(Nd-Nx, Nx)
    ]
    D2 = [
        spzeros(Nx, Np);
        spzeros(Nd-Nx, Nx + Nu) blockdiag([sparse(k == N-1 ? C[k] : [C[k] D[k]]) for k = 1:N-1]...)
    ]
    D = D1 + D2
    A = [H D'; D sparse(-I*reg, Nd, Nd)]

    bp = push!([Vector([qk; rk]) for (qk,rk) in zip(q,r)], q[end])
    bd = [data.x0; data.d]
    b = stack_vectors([bp; bd])
    return A,-b
end

function permvec(data::LQRData{Nx,Nu}) where {Nx,Nu}
    N = length(data.Q)
    Np = N*Nx + (N-1)*Nu
    Nd = N*Nx
    iz = [(1:Nx+Nu * (k < N)) .+ (k-1)*(Nx+Nu) for k = 1:N]
    iy = [Np .+ (1:Nx) .+ (k-1)*(Nx) for k = 1:N]
    stack_vectors([[y;z] for (z,y) in zip(iz,iy)])
end