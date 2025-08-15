"""
    Initialize a problem for a stochastic Schrödinger simulation.

    ω0s:            energies of states (angular units)
    ωs:             laser frequency components (angular units)
    sats:           saturation parameters for each frequency component
    pols:           polarizations of laser frequency components
    ψ0:             initial value for ψ
    d:              transition dipole array
    m:              mass of particle
    Γ:              linewidth of transition
    k:              k-value of transition
    params:         custom parameters to be passed to the solver
    add_terms_dψ:   function to add custom terms to dψ, with signature `add_terms_dψ(dψ, ψ, p, r, t)``
"""
mutable struct prob_params{T,N1,N2,N3,F1<:Function,F2<:Function}
    sim_params::T
    n_g::Int64
    n_e::Int64
    n_states::Int64
    u0::Vector{Float64}
    Γ::Float64
    ωs::Vector{Vector{Float64}}
    ω0s::MVector{N1, Float64}
    eiω0ts::StructVector{ComplexF64, @NamedTuple{re::MVector{N1, Float64}, im::MVector{N1, Float64}}, Int64}
    ϕs::Vector{Matrix{Float64}}
    as::Vector{StructArray{ComplexF64, 2, @NamedTuple{re::Matrix{Float64}, im::Matrix{Float64}}, Int64}}
    rs::Vector{Matrix{Float64}}
    kEs::Vector{StructArray{ComplexF64, 2, @NamedTuple{re::Matrix{Float64}, im::Matrix{Float64}}, Int64}}
    E::Vector{StructVector{ComplexF64, @NamedTuple{re::MVector{3, Float64}, im::MVector{3, Float64}}, Int64}}
    ϵs::Vector{StructArray{ComplexF64, 3, @NamedTuple{re::Array{Float64, 3}, im::Array{Float64, 3}}, Int64}}
    denom::Float64
    ψ::StructVector{ComplexF64, @NamedTuple{re::MVector{N1, Float64}, im::MVector{N1, Float64}}, Int64}
    dψ::StructVector{ComplexF64, @NamedTuple{re::MVector{N1, Float64}, im::MVector{N1, Float64}}, Int64}
    ψ_q::Vector{StructArray{ComplexF64, 2, @NamedTuple{re::MMatrix{N1, 3, Float64, N2}, im::MMatrix{N1, 3, Float64, N2}}, Int64}}
    d::MArray{Tuple{N1, N1, 3}, Float64, 3, N3}
    d_ge::Vector{Array{Float64, 3}}
    d_eg::Vector{Array{Float64, 3}}
    F::MVector{3, Float64}
    d_exp::StructVector{ComplexF64, @NamedTuple{re::MVector{3, Float64}, im::MVector{3, Float64}}, Int64}
    r::MVector{3, Float64}
    ks::Vector{Float64}
    r_idx::Int64
    v_idx::Int64
    F_idx::Int64
    m::Float64
    add_terms_dψ::F1
    update_params::F2
    decay_dist::Exponential{Float64}
    time_to_decay::Float64
    last_decay_time::Float64
    n_scatters::Float64
    diffusion_constant::MVector{3, Float64}
    add_spontaneous_decay_kick::Bool
    sats::Vector{Vector{Float64}}
    coupling_idxs::Vector{Vector{UnitRange{Int64}}}
end

function initialize_prob_multiple(
        sim_type,
        ω0s,
        ωs,
        sats,
        pols,
        ks,
        beam_radius,
        d,
        m,
        Γ,
        k,
        sim_params,
        update_params,
        add_terms_dψ,
        n_g,
        n_e,
        coupling_idxs
    )

    # get some initial constants
    n_states = length(ω0s)
    n_freqs = [length(ωs[i]) for i ∈ eachindex(coupling_idxs)]

    denom = sim_type((beam_radius*k)^2/2)

    eiω0ts = zeros(Complex{sim_type},n_states)

    k_dirs = 6

    as = [zeros(Complex{sim_type},k_dirs,n_freqs[i]) for i ∈ eachindex(coupling_idxs)]
    ϕs = [zeros(sim_type,3,3+n_freqs[i]) for i ∈ eachindex(coupling_idxs)]
    rs = [zeros(sim_type,2,3) for i ∈ eachindex(coupling_idxs)]
    kEs = [zeros(Complex{sim_type}, k_dirs, 3) for _ ∈ eachindex(coupling_idxs)]

    # define polarization array
    ϵs = [zeros(Complex{sim_type},k_dirs,n_freqs[i],3) for i ∈ eachindex(coupling_idxs)]
    for i ∈ eachindex(coupling_idxs)
        for j ∈ eachindex(pols[i])
            pol = pols[i][j]
            ϵs[i][1,j,:] .= rotate_pol(pol, x̂)
            ϵs[i][2,j,:] .= rotate_pol(pol, ŷ)
            ϵs[i][3,j,:] .= rotate_pol(flip(pol), ẑ)
            ϵs[i][4,j,:] .= rotate_pol(pol, -x̂)
            ϵs[i][5,j,:] .= rotate_pol(pol, -ŷ)
            ϵs[i][6,j,:] .= rotate_pol(flip(pol), -ẑ)
        end
    end

    # arrays related to energies
    # as = [StructArray(MMatrix{k_dirs,n_freqs[i]}(as[i])) for i ∈ eachindex(ωs)]
    # ωs = [MVector{size(ωs[i])...}(ωs[i]) for i ∈ eachindex(ωs)]
    # ϕs = [MMatrix{size(ϕs[i])...}(ϕs[i]) for i ∈ eachindex(ωs)]
    # rs = [MMatrix{size(rs[i])...}(rs[i]) for i ∈ eachindex(ωs)]
    # kEs = [StructArray(MMatrix{size(kEs[i])...}(kEs[i])) for i ∈ eachindex(ωs)]
    # ϵs = [StructArray(MArray{Tuple{k_dirs,n_freqs[i],3}}(ϵs[i])) for i ∈ eachindex(ωs)]

    as = [StructArray(as[i]) for i ∈ eachindex(ωs)]
    kEs = [StructArray(kEs[i]) for i ∈ eachindex(ωs)]
    ϵs = [StructArray(ϵs[i]) for i ∈ eachindex(ωs)]

    # arrays related to state energies
    ω0s = MVector{size(ω0s)...}(ω0s)
    eiω0ts = StructArray(MVector{size(eiω0ts)...}(eiω0ts))

    ψ = zeros(Complex{sim_type}, n_states)
    ψ = MVector{size(ψ)...}(ψ)
    ψ = StructArray(ψ)

    dψ = deepcopy(ψ)

    E = zeros(Complex{sim_type},3)
    E = MVector{size(E)...}(E)
    E = StructArray(E)
    E = [deepcopy(E) for _ ∈ eachindex(coupling_idxs)]

    # note that we take the negative to ensure that the Hamiltonian is -d⋅E
    d_ge = [.-real.(d[coupling_idxs[i][1],coupling_idxs[i][2],:]) for i ∈ eachindex(coupling_idxs)]
    # d_ge = [MArray{Tuple{size(d_ge[i])...}}(d_ge[i]) for i ∈ eachindex(d_ge)]
    d_eg = [permutedims(d_ge[i],(2,1,3)) for i ∈ eachindex(d_ge)]

    d = MArray{Tuple{n_states,n_states,3}}(sim_type.(real.(d)))

    # ψ_q = [zeros(Complex{sim_type}, size(d_ge[i],1)+size(d_ge[i],2)) for i ∈ eachindex(d)]
    # ψ_q = [MArray{Tuple{size(ψ_q[i])...,3}}(zeros(Complex{sim_type},size(ψ_q[i])...,3)) for i ∈ eachindex(d)]
    # ψ_q = [StructArray(ψ_q[i]) for i ∈ eachindex(d)]

    ψ_q = [MArray{Tuple{size(ψ)...,3}}(zeros(Complex{sim_type},size(ψ)...,3)) for i ∈ eachindex(coupling_idxs)]
    ψ_q = [StructArray(ψ_q[i]) for i ∈ eachindex(coupling_idxs)]

    F = MVector{3,sim_type}(zeros(3))

    d_exp = zeros(Complex{sim_type},3)
    d_exp = MVector{size(d_exp)...}(d_exp)
    d_exp = StructArray(d_exp)

    r = MVector{3}(zeros(sim_type,3))
    r_idx = 2n_states + n_e
    v_idx = r_idx + 3
    F_idx = v_idx + 3

    decay_dist = Exponential(one(sim_type))
    last_decay_time = zero(sim_type)

    diffusion_constant = MVector{3,sim_type}(zeros(3))
    add_spontaneous_decay_kick = false

    n_scatters = zero(sim_type)

    u0 = sim_type.([zeros(n_states)..., zeros(n_states)..., zeros(n_e)..., zeros(3)..., zeros(3)..., zeros(3)..., zeros(3)...])
    u0[1] = 1.0

    p = prob_params(
        sim_params,
        n_g,
        n_e,
        n_states,
        u0,
        Γ,
        ωs,
        ω0s,
        eiω0ts,
        ϕs,
        as,
        rs,
        kEs,
        E,
        ϵs,
        denom,
        ψ,
        dψ,
        ψ_q,
        d,
        d_ge,
        d_eg,
        F,
        d_exp,
        r,
        ks,
        r_idx,
        v_idx,
        F_idx,
        m,
        add_terms_dψ,
        update_params,
        decay_dist,
        rand(decay_dist),
        last_decay_time,
        n_scatters,
        diffusion_constant,
        add_spontaneous_decay_kick,
        sats,
        coupling_idxs
    )

    # p = MutableNamedTuple(
    #     sim_params=sim_params,
    #     n_g=n_g,
    #     n_e=n_e,
    #     n_states=n_states,
    #     u0=u0,
    #     Γ=Γ,
    #     ωs=ωs,
    #     ω0s=ω0s,
    #     eiω0ts=eiω0ts,
    #     ϕs=ϕs,
    #     as=as,
    #     rs=rs,
    #     kEs=kEs,
    #     E=E,
    #     ϵs=ϵs,
    #     denom=denom,
    #     ψ=ψ,
    #     dψ=dψ,
    #     ψ_q=ψ_q,
    #     d=d,
    #     d_ge=d_ge,
    #     d_eg=d_eg,
    #     F=F,
    #     d_exp=d_exp,
    #     r=r,
    #     ks=ks,
    #     r_idx=r_idx,
    #     v_idx=v_idx,
    #     F_idx=F_idx,
    #     m=m,
    #     add_terms_dψ=add_terms_dψ,
    #     update_params=update_params,
    #     decay_dist=decay_dist,
    #     time_to_decay=rand(decay_dist),
    #     last_decay_time=last_decay_time,
    #     n_scatters=n_scatters,
    #     diffusion_constant=diffusion_constant,
    #     add_spontaneous_decay_kick=add_spontaneous_decay_kick,
    #     sats=sats,
    #     coupling_idxs=coupling_idxs
    # )

    return p
end
export initialize_prob_multiple

function ψ_fast_multiple!(du, u, p, t)
    
    normalize_u!(u, p.n_states)

    update_r!(u, p.r, p.r_idx)

    p.update_params(p, p.r, t)

    update_ψ!(p.ψ, u, p.n_states)

    update_eiωt_new!(p.eiω0ts, p.ω0s, t)

    Heisenberg_turbo_state!(p.ψ, p.eiω0ts, -1)

    zero_array!(p.dψ)
    @inline for i ∈ eachindex(p.coupling_idxs)

        update_fields_fast_multiple!(p.denom, p.rs[i], p.ϕs[i], p.ωs[i], p.as[i], p.sats[i], p.kEs[i], p.ϵs[i], p.E[i], p.ks[i], p.r, t)

        update_ψq_multiple!(p.ψ_q[i], p.d_ge[i], p.d_eg[i], p.ψ, p.coupling_idxs[i][1][1]-1, p.coupling_idxs[i][2][1]-1)
        
        update_dψ_multiple!(p.dψ, p.ψ_q[i], p.E[i], p.d_ge[i], p.coupling_idxs[i][1][1]-1, p.coupling_idxs[i][2][1]-1)

    end

    update_d_exp_multiple!(p.d_exp, p.ψ, p.ψ_q[1])

    update_force!(p.F, p.d_exp, p.kEs[1])

    p.add_terms_dψ(p.dψ, p.ψ, p, p.r, t) # custom terms to add to dψ

    Heisenberg_turbo_state!(p.dψ, p.eiω0ts, +1)

    update_du!(du, u, p.dψ, p.ψ, p.n_states, p.n_g, p.n_e, p.r_idx, p.F, p.v_idx, p.F_idx, p.m)

    return nothing
end
export ψ_fast_multiple!

function ψ_fast_multiple_ballistic!(du, u, p, t)
    
    normalize_u!(u, p.n_states)

    update_r!(u, p.r, p.r_idx)

    p.update_params(p, p.r, t)

    update_ψ!(p.ψ, u, p.n_states)

    update_eiωt_new!(p.eiω0ts, p.ω0s, t)

    Heisenberg_turbo_state!(p.ψ, p.eiω0ts, -1)

    zero_array!(p.dψ)
    @inline for i ∈ eachindex(p.coupling_idxs)

        update_fields_fast_multiple!(p.denom, p.rs[i], p.ϕs[i], p.ωs[i], p.as[i], p.sats[i], p.kEs[i], p.ϵs[i], p.E[i], p.ks[i], p.r, t)

        update_ψq_multiple!(p.ψ_q[i], p.d_ge[i], p.d_eg[i], p.ψ, p.coupling_idxs[i][1][1]-1, p.coupling_idxs[i][2][1]-1)
        
        update_dψ_multiple!(p.dψ, p.ψ_q[i], p.E[i], p.d_ge[i], p.coupling_idxs[i][1][1]-1, p.coupling_idxs[i][2][1]-1)
        
    end

    update_d_exp_multiple!(p.d_exp, p.ψ, p.ψ_q[1])

    update_force!(p.F, p.d_exp, p.kEs[1])

    p.add_terms_dψ(p.dψ, p.ψ, p, p.r, t) # custom terms to add to dψ

    Heisenberg_turbo_state!(p.dψ, p.eiω0ts, +1)

    update_du_ballistic!(du, u, p.dψ, p.ψ, p.n_states, p.n_g, p.n_e, p.r_idx, p.F, p.v_idx, p.F_idx, p.m)

    return nothing
end
export ψ_fast_multiple_ballistic!

@inline function update_d_exp_multiple!(d_exp, ψ, ψ_q)
    g_idx = 4
    @turbo for q ∈ 1:3
        re = zero(eltype(ψ.re))
        im = zero(eltype(ψ.im))
        for i ∈ 1:12
            ψ_re = ψ.re[i+g_idx] # take conjugate
            ψ_im = -ψ.im[i+g_idx]
            ψq_re = ψ_q.re[i+g_idx,q]
            ψq_im = ψ_q.im[i+g_idx,q]
            re += ψ_re * ψq_re - ψ_im * ψq_im
            im += ψ_re * ψq_im + ψq_re * ψ_im
        end
        d_exp.re[q] = re
        d_exp.im[q] = im
    end
    return nothing
end

function zero_array!(A)
    @turbo for i ∈ eachindex(A)
        A.re[i] = zero(eltype(A.re))
        A.im[i] = zero(eltype(A.im))
    end
    return nothing
end

function update_ψq_multiple!(ψ_q, d_ge, d_eg, ψ, g_idx, e_idx)
    # ground states
    @turbo for q ∈ 1:3
        for i ∈ axes(d_ge,1)
            ψq_re_i = zero(eltype(ψ_q.re))
            ψq_im_i = zero(eltype(ψ_q.im))
            for j ∈ axes(d_ge,2)
                d_q_ij = d_ge[i,j,q]
                ψ_re_j = ψ.re[j+e_idx]
                ψ_im_j = ψ.im[j+e_idx]
                ψq_re_i += d_q_ij * ψ_re_j
                ψq_im_i += d_q_ij * ψ_im_j
            end
            ψ_q.re[i+g_idx,q] = ψq_re_i
            ψ_q.im[i+g_idx,q] = ψq_im_i
        end
    end
    # excited states
    @turbo for q ∈ 1:3
        for i ∈ axes(d_eg,1)
            ψq_re_i = zero(eltype(ψ_q.re))
            ψq_im_i = zero(eltype(ψ_q.im))
            for j ∈ axes(d_eg,2)
                d_q_ij = d_eg[i,j,q]
                ψ_re_j = ψ.re[j+g_idx]
                ψ_im_j = ψ.im[j+g_idx]
                ψq_re_i += d_q_ij * ψ_re_j
                ψq_im_i += d_q_ij * ψ_im_j
            end
            ψ_q.re[i+e_idx,q] = ψq_re_i
            ψ_q.im[i+e_idx,q] = ψq_im_i
        end
    end
    return nothing
end

function update_dψ_multiple!(dψ, ψ_q, E, d_ge, g_idx, e_idx)
    @turbo for i ∈ axes(d_ge,1)
        dψ_i_re = zero(eltype(dψ.re))
        dψ_i_im = zero(eltype(dψ.im))
        for q ∈ 1:3
            E_q_re = E.re[q]
            E_q_im = E.im[q]
            ψ_q_re = ψ_q.re[i+g_idx,q]
            ψ_q_im = ψ_q.im[i+g_idx,q]
            
            dψ_i_re += ψ_q_re * E_q_re - ψ_q_im * E_q_im
            dψ_i_im += ψ_q_re * E_q_im + ψ_q_im * E_q_re
        end
        dψ.re[i+g_idx] += dψ_i_im # multiply by -im
        dψ.im[i+g_idx] += -dψ_i_re
    end
    @turbo for i ∈ axes(d_ge,2)
        dψ_i_re = zero(eltype(dψ.re))
        dψ_i_im = zero(eltype(dψ.im))
        for q ∈ 1:3
            E_q_re = E.re[q]
            E_q_im = -E.im[q] # conjugate for the excited states
            ψ_q_re = ψ_q.re[i+e_idx,q]
            ψ_q_im = ψ_q.im[i+e_idx,q]
            
            dψ_i_re += ψ_q_re * E_q_re - ψ_q_im * E_q_im
            dψ_i_im += ψ_q_re * E_q_im + ψ_q_im * E_q_re
        end
        dψ.re[i+e_idx] += dψ_i_im # multiply by -im
        dψ.im[i+e_idx] += -dψ_i_re
    end
    return nothing
end
export update_dψ!