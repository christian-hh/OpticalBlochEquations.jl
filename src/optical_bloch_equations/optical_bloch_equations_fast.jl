function initialize_obes_prob(
    sim_type,
    ω0s,
    ωs,
    sats,
    pols,
    beam_radius,
    d,
    m,
    Γ,
    k,
    sim_params,
    update_params,
    add_terms_dρ,
    should_round_freqs=true,
    freq_res=1e-2
)

    # get some initial constants
    n_states = length(ω0s)
    n_g = find_n_g(d)
    n_e = n_states - n_g
    n_freqs = length(ωs)

    # set the integer type for the simulation
    intT = get_int_type(sim_type)

    denom = sim_type((beam_radius*k)^2/2)

    eiω0ts = zeros(Complex{sim_type},n_states)

    k_dirs = 6

    as = zeros(Complex{sim_type},k_dirs,n_freqs)
    ϕs = zeros(sim_type,3,3+n_freqs)
    rs = zeros(sim_type,2,3)
    kEs = zeros(Complex{sim_type}, k_dirs, 3)

    idxs = intT.(reshape(collect(1:(3n_freqs)),3,n_freqs))

    # define polarization array
    ϵs = zeros(Complex{sim_type},k_dirs,n_freqs,3)
    for i ∈ eachindex(pols)
        pol = pols[i]
        ϵs[1,i,:] .= rotate_pol(pol, x̂)
        ϵs[2,i,:] .= rotate_pol(pol, ŷ)
        ϵs[3,i,:] .= rotate_pol(flip(pol), ẑ)
        ϵs[4,i,:] .= rotate_pol(pol, -x̂)
        ϵs[5,i,:] .= rotate_pol(pol, -ŷ)
        ϵs[6,i,:] .= rotate_pol(flip(pol), -ẑ)
    end

    # arrays related to energies
    as = StructArray(MMatrix{k_dirs,n_freqs}(as))
    ωs = MVector{size(ωs)...}(ωs)
    ϕs = MMatrix{size(ϕs)...}(ϕs)
    rs = MMatrix{size(rs)...}(rs)
    kEs = StructArray(MMatrix{size(kEs)...}(kEs))
    ϵs = StructArray(MArray{Tuple{k_dirs,n_freqs,3}}(ϵs))
    idxs = MMatrix{size(idxs)...}(idxs)

    # arrays related to state energies
    ω0s = MVector{size(ω0s)...}(ω0s)
    eiω0ts = StructArray(MVector{size(eiω0ts)...}(eiω0ts))

    ρ = zeros(Complex{sim_type}, n_states, n_states)
    ρ = MMatrix{size(ρ)...}(ρ)
    ρ = StructArray(ρ)

    dρ = deepcopy(ρ)
    dρ_adj = deepcopy(ρ)

    decay_terms = deepcopy(ρ)

    # ρ_q = deepcopy(ρ)
    ρ_q = MArray{Tuple{size(ρ)...,3}}(zeros(Complex{sim_type},size(ρ)...,3))
    ρ_q = StructArray(ρ_q)

    H = StructArray(MArray{Tuple{size(ρ)...}}(zeros(Complex{sim_type},size(ρ)...)))

    # dge_ρ = StructArray(ρ_q)
    # deg_ρ = StructArray(ρ_q)

    dge_ρ = zeros(Complex{sim_type},12,16,3)
    deg_ρ = zeros(Complex{sim_type},4,16,3)
    dge_ρ = MArray{Tuple{size(dge_ρ)...}}(dge_ρ)
    deg_ρ = MArray{Tuple{size(deg_ρ)...}}(deg_ρ)
    dge_ρ = StructArray(dge_ρ)
    deg_ρ = StructArray(deg_ρ)

    E_total = zeros(Complex{sim_type},3)
    E_total = MVector{size(E_total)...}(E_total)
    E_total = StructArray(E_total)

    # note that we take the negative to ensure that the Hamiltonian is -d⋅E
    d_ge = sim_type.(real.(-d[1:n_g,(n_g+1):n_states,:]))
    d_ge = MArray{Tuple{size(d_ge)...}}(d_ge)
    d_eg = permutedims(d_ge,(2,1,3))

    F = MVector{3,sim_type}(zeros(3))

    d = MArray{Tuple{n_states,n_states,3}}(sim_type.(real.(d)))

    d_exp = zeros(Complex{sim_type},3)
    d_exp = MVector{size(d_exp)...}(d_exp)
    d_exp = StructArray(d_exp)

    d_exp_split = zeros(Complex{sim_type},3,2)
    d_exp_split = MMatrix{size(d_exp_split)...}(d_exp_split)
    d_exp_split = StructArray(d_exp_split)

    r = MVector{3}(zeros(sim_type,3))
    r_idx = 2n_states^2 + n_e
    v_idx = r_idx + 3
    F_idx = v_idx + 3

    n_scatters = zero(sim_type)

    u0 = sim_type.([zeros(n_states^2)..., zeros(n_states^2)..., zeros(4)..., zeros(3)..., zeros(3)..., zeros(3)..., zeros(3)...])
    u0[1] = 1.0
    
    # Create jumps corresponding to spontaneous decay
    Js = Array{Jump}(undef, 0)
    ds = [Float64[], Float64[], Float64[]]
    ds_state1 = [Int64[], Int64[], Int64[]]
    ds_state2 = [Int64[], Int64[], Int64[]]
    for s′ in 1:n_states, s in s′:n_states, q in (-1,0,+1)
        dme = d[s′, s, q+2]
        if abs(dme) > 1e-10 && (ω0s[s′] < ω0s[s]) # only energy-allowed jumps are generated
            push!(ds_state1[q+2], s)
            push!(ds_state2[q+2], s′)
            push!(ds[q+2], dme)
            J = Jump(s, s′, q, dme)
            push!(Js, J)
        end
    end
    ds = [ds[1], ds[2], ds[3]]

    n_transitions = length(ds[1])
    ds = [SVector{n_transitions}(ds[1]), SVector{n_transitions}(ds[2]), SVector{n_transitions}(ds[3])]
    ds_state1 = [SVector{n_transitions}(ds_state1[1]), SVector{n_transitions}(ds_state1[2]), SVector{n_transitions}(ds_state1[3])]
    ds_state2 = [SVector{n_transitions}(ds_state2[1]), SVector{n_transitions}(ds_state2[2]), SVector{n_transitions}(ds_state2[3])]

    return MutableNamedTuple(
        u0=u0,
        Γ=Γ,
        ωs=ωs,
        ω0s=ω0s,
        eiω0ts=eiω0ts,
        ϕs=ϕs,
        as=as,
        rs=rs,
        kEs=kEs,
        E_total=E_total,
        ϵs=ϵs,
        idxs=idxs,
        denom=denom,
        ρ=ρ,
        dρ=dρ,
        dρ_adj=dρ_adj,
        ρ_q=ρ_q,
        dge_ρ=dge_ρ,
        deg_ρ=deg_ρ,
        sim_params=sim_params,
        d_ge=d_ge,
        d_eg=d_eg,
        F=F,
        d=d,
        d_exp=d_exp,
        d_exp_split=d_exp_split,
        r=r,
        r_idx=r_idx,
        v_idx=v_idx,
        F_idx=F_idx,
        n_g=n_g,
        n_e=n_e,
        n_states=n_states,
        m=m,
        add_terms_dρ=add_terms_dρ,
        update_params=update_params,
        n_scatters=n_scatters,
        sats=sats,
        Js=Js,
        ds=ds, ds_state1=ds_state1, ds_state2=ds_state2, decay_terms=decay_terms,
        H=H
    )
end
export initialize_obes_prob

"""
    Time step update function for optical Bloch equations (OBEs) simulation.
"""
function update_obes_fast!(du, u, p, t)

    update_r!(u, p.r, p.r_idx)

    p.update_params(p, p.r, t)

    update_ρ!(p.ρ, u, p.n_states)

    update_fields_fast!(p, p.r, t)
    
    update_eiωt_new!(p.eiω0ts, p.ω0s, t)

    Heisenberg_turbo!(p.ρ, p.eiω0ts, -1)

    update_H_fast!(p.H, p.d_ge, p.E_total, p.n_g, p.n_e)

    im_commutator_fast!(p.dρ, p.H, p.ρ, p.dρ_adj)

    # update_dρ_exp!(p.d_exp, p.ρ_q)

    # update_force!(p.F, p.d_exp, p.kEs)

    obe_decay_terms!(p.dρ, p.ρ, p.ds, p.ds_state1, p.ds_state2)

    # p.add_terms_dρ(p.dρ, p.ρ, p, p.r, t) # custom terms to add to dρ

    Heisenberg_turbo!(p.dρ, p.eiω0ts, +1)

    update_du_obes!(du, u, p.dρ, p.ρ, p.n_states, p.n_g, p.n_e, p.r_idx, p.F, p.v_idx, p.F_idx, p.m)

    return nothing
end
export update_obes_fast!

function im_commutator_fast!(C, A, B, tmp)
    # Multiply C = A * B -- (H * ρ)
    # Subtract adjoint C = A * B - B† * A†
    @turbo for j ∈ axes(B,2), i ∈ axes(A,1)
        Cre = zero(eltype(C.re))
        Cim = zero(eltype(C.im))
        for k ∈ axes(A,2)
            Aik_re = A.re[i,k]
            Aik_im = A.im[i,k]
            Bkj_re = B.re[k,j]
            Bkj_im = B.im[k,j]
            Cre += Aik_re * Bkj_re - Aik_im * Bkj_im
            Cim += Aik_re * Bkj_im + Aik_im * Bkj_re
        end
        # mul by -im
        C.re[i,j] = +Cim
        C.im[i,j] = -Cre
    end

    @turbo for j ∈ axes(C,2), i ∈ axes(C,1)
        tmp.re[j,i] = +C.re[i,j]
        tmp.im[j,i] = -C.im[i,j]
    end

    # add adjoint
    @turbo for j ∈ axes(C,2), i ∈ axes(C,1)
        C.re[i,j] += +tmp.re[i,j]
        C.im[i,j] += +tmp.im[i,j]
    end

    return nothing
end
export im_commutator_fast!

@inline function update_ρ!(ρ, u, n_states)
    c_idx = n_states^2
    @turbo for i ∈ eachindex(ρ)
        ρ.re[i] = u[i]
        ρ.im[i] = u[i+c_idx]
    end
    return nothing
end
export update_ρ!

function update_H_fast!(H, d_ge, E, n_g, n_e)

    @turbo for j ∈ axes(d_ge,2), i ∈ axes(d_ge,1)
        H_re = zero(eltype(H.re))
        H_im = zero(eltype(H.re))
        for q ∈ 1:3
            E_q_re = E.re[q]
            E_q_im = E.im[q]
            H_re += E_q_re * d_ge[i,j,q]
            H_im += E_q_im * d_ge[i,j,q]
        end
        H.re[i,j+n_g] = H_re
        H.im[i,j+n_g] = H_im
    end

    # add hermitian conjugate
    @turbo for j ∈ axes(d_ge,2), i ∈ axes(d_ge,1)
        H.re[j+n_g,i] = H.re[i,j+n_g]
        H.im[j+n_g,i] = -H.im[i,j+n_g]
    end

    # add decays -i∑ᵢJᵢ†Jᵢ/2
    @inbounds @fastmath for j ∈ 1:n_e
        H.im[j+n_g,j+n_g] = -1/2
    end

    return nothing
end
export update_H_fast!

"""
    Take the trace of Hρ to get the expectation value of H.
"""
@inline function update_H_exp!(d_exp, ρ_q)
    for q ∈ 1:3
        re = zero(eltype(ρ_q.re))
        for i ∈ axes(ρ_q,1)
            re += ρ_q.re[i,i,q]
        end
        d_exp[q] = re
    end
    return nothing
end
export update_dρ_exp!

@inline function update_du_obes!(du, u, dρ, ρ, n_states, n_g, n_e, r_idx, F, v_idx, F_idx, m)

    c_idx = n_states^2

    @turbo for i ∈ eachindex(dρ)
        du[i] = dρ.re[i]
        du[i+c_idx] = dρ.im[i]
    end

    @inbounds @fastmath for k ∈ 1:3
        du[r_idx+k] = u[v_idx+k]
        du[v_idx+k] = F[k] / m
        u[F_idx+k] = F[k]
        du[F_idx+k+3] = F[k] # integrated force
    end

    # integrated excited state population
    # @turbo for i ∈ 1:4 # combine this with loop below?
    #     ρ_i_pop = u[n_g+i]
    #     du[i+2c_idx] = ρ_i_pop
    # end

    # @turbo for i ∈ 1:n_e
    #     # non-hermitian part of Hamiltonian, -im/2, but multiplied by -im also
    #     du[n_g+i] -= u[n_g+i]/2
    #     du[n_states+n_g+i] -= u[n_states+n_g+i]/2
    # end

    return nothing
end

"""
    Decay terms like ∑ᵢJᵢρJᵢ† terms, where J = |g⟩⟨e|.
    Could probably use the dρ matrix as a starting point here.
"""
function obe_decay_terms!(dρ, ρ, ds, ds_state1, ds_state2)

    @inbounds @fastmath for q ∈ 1:3
        ds_q = ds[q]
        ds_state1_q = ds_state1[q]
        ds_state2_q = ds_state2[q]
        for i ∈ eachindex(ds_q)
            s1_i = ds_state1_q[i]
            s2_i = ds_state2_q[i]
            ds_q_i = ds_q[i]

            # ground state accumulation
            dρ.re[s2_i, s2_i] += ds_q_i^2 * ρ.re[s1_i, s1_i]

            for j ∈ (i+1):length(ds_q)
                ds_q_j = ds_q[j]
                s1_j = ds_state1_q[j]
                s2_j = ds_state2_q[j]
                ds_q_ij = ds_q_i * ds_q_j
                val_re = ds_q_ij * ρ.re[s1_i, s1_j]
                val_im = ds_q_ij * ρ.im[s1_i, s1_j]
                dρ.re[s2_i, s2_j] += val_re
                dρ.im[s2_i, s2_j] += val_im
                dρ.re[s2_j, s2_i] += val_re
                dρ.im[s2_j, s2_i] -= val_im
            end
        end

    end

    return nothing
end
export obe_decay_terms!

#################################################

@inline function update_ρq_full!(ρ_q, d_ge, d_eg, E, ρ, n_g)
    C = ρ_q
    A1 = d_ge
    A2 = d_eg
    B = ρ
    @turbo for q ∈ 1:3
        E_q_re = E.re[q]
        E_q_im = E.im[q]
        for m ∈ axes(A1,1), n ∈ axes(B,2)
            Cmn_re = zero(eltype(C))
            Cmn_im = zero(eltype(C))
            for k ∈ axes(A1,2) # iterate excited states
                A1_mk_re = E_q_re * A1[m,k,q]
                A1_mk_im = E_q_im * A1[m,k,q]
                B_kn_re = B.re[k+n_g,n]
                B_kn_im = B.im[k+n_g,n]
                Cmn_re += A1_mk_re * B_kn_re - A1_mk_im * B_kn_im
                Cmn_im += A1_mk_re * B_kn_im + A1_mk_im * B_kn_re
            end
            C.re[m,n,q] = Cmn_re
            C.im[m,n,q] = Cmn_im
        end
    end
    @turbo for q ∈ 1:3
        E_q_re = E.re[q]
        E_q_im = E.im[q]
        for m ∈ axes(A2,1), n ∈ axes(B,2)
            Cmn_re = zero(eltype(D))
            Cmn_im = zero(eltype(D))
            for k ∈ axes(A2,2)
                A2_mk_re = E_q_re * A2[m,k,q]
                A2_mk_im = E_q_im * A2[m,k,q]
                B_kn_re = B.re[k,n]
                B_kn_im = B.im[k,n]
                Cmn_re += A2_mk_re * B_kn_re - A2_mk_im * B_kn_im
                Cmn_im += A2_mk_re * B_kn_im + A2_mk_im * B_kn_re
            end
            C.re[m+n_g,n,q] = Cmn_re
            C.im[m+n_g,n,q] = Cmn_im
        end
    end
    return nothing
end
export update_ρq_full!

@inline function update_dρ!(dρ, dρ_adj, ρ_q, E)
    @turbo for j ∈ axes(dρ,2), i ∈ axes(dρ,1)
        dρ_i_re = zero(eltype(dρ.re))
        dρ_i_im = zero(eltype(dρ.im))
        for q ∈ 1:3
            E_q_re = E.re[q]
            E_q_im = E.im[q]
            ρ_q_re = ρ_q.re[i,j,q]
            ρ_q_im = ρ_q.im[i,j,q]
            
            dρ_i_re += ρ_q_re * E_q_re - ρ_q_im * E_q_im
            dρ_i_im += ρ_q_re * E_q_im + ρ_q_im * E_q_re
        end
        dρ.re[i,j] = dρ_i_re
        dρ.im[i,j] = dρ_i_im
    end

    # add negative adjoint to make commutator ρ*d
    @turbo for j ∈ axes(dρ,2), i ∈ axes(dρ,1)
        dρ_adj.re[i,j] = +dρ.re[j,i]
        dρ_adj.im[i,j] = -dρ.im[j,i]
    end

    # subtract it
    @turbo for i ∈ eachindex(dρ)
        dρ.re[i] -= dρ_adj.re[i]
        dρ.im[i] -= dρ_adj.im[i]
    end

    return nothing

end
export update_dρ!