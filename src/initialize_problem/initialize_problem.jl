function initialize_p(
    sim_type,
    energies,
    freqs,
    sats,
    pols,
    beam_radius,
    d,
    m,
    Γ,
    k,
    sim_params,
    update_params,
    add_terms_dψ
)

    # get some initial constants
    n_states = length(energies)
    n_ground = find_n_g(d)
    n_excited = n_states - n_ground
    n_freqs = length(freqs)
    n_ks = 6

    # set the integer type for the simulation
    sim_type_int = get_int_type(sim_type)

    denom = sim_type((beam_radius*k)^2/2)

    eiω0ts = zeros(Complex{sim_type}, n_states)

    # amplitudes of electric fields
    as = zeros(Complex{sim_type}, n_ks, n_freqs)
    as = StructArray(MMatrix{n_ks, n_freqs}(as))

    # define polarization array
    ϵs = zeros(Complex{sim_type}, n_ks, n_freqs, 3)
    for i ∈ eachindex(pols)
        pol = pols[i]
        ϵs[1,i,:] .= rotate_pol(pol, x̂)
        ϵs[2,i,:] .= rotate_pol(pol, ŷ)
        ϵs[3,i,:] .= rotate_pol(flip(pol), ẑ)
        ϵs[4,i,:] .= rotate_pol(pol, -x̂)
        ϵs[5,i,:] .= rotate_pol(pol, -ŷ)
        ϵs[6,i,:] .= rotate_pol(flip(pol), -ẑ)
    end
    ϵs = StructArray(MArray{Tuple{n_ks, n_freqs, 3}}(ϵs))

    # arrays related to energies
    energies = MVector{size(energies)...}(energies)

    # phases of laser frequencies
    ϕs = zeros(sim_type, 3, 3 + n_freqs)
    ϕs = MMatrix{size(ϕs)...}(ϕs)

    # position
    rs = zeros(sim_type,2,3)
    rs = MMatrix{size(rs)...}(rs)

    #
    kEs = zeros(Complex{sim_type}, k_dirs, 3)
    kEs = StructArray(MMatrix{size(kEs)...}(kEs))

    # total electric field
    E_total = zeros(Complex{sim_type}, 3)
    E_total = MVector{size(E_total)...}(E_total)
    E_total = StructArray(E_total)

    # note that we take the negative to ensure that the Hamiltonian is -d⋅E
    d_ge = sim_type.(real.(-d[1:n_g,(n_g+1):n_states,:]))
    d_ge = MArray{Tuple{size(d_ge)...}}(d_ge)
    d_eg = permutedims(d_ge,(2,1,3))

    # expectation value of d_exp (only the positive component)
    d_exp = zeros(Complex{sim_type}, 3)
    d_exp = MVector{size(d_exp)...}(d_exp)
    d_exp = StructArray(d_exp)

    # force
    F = MVector{3, sim_type}(zeros(3))

    n_scatters = zero(sim_type_int)

    return MutableNamedTuple(
        Γ=Γ,
        energies=energies,
        freqs=freqs,
        sats=sats,
        eiω0ts=eiω0ts,
        ϕs=ϕs,
        as=as,
        rs=rs,
        kEs=kEs,
        E_total=E_total,
        ϵs=ϵs,
        denom=denom,
        sim_params=sim_params,
        d_ge=d_ge,
        d_eg=d_eg,
        F=F,
        d_exp=d_exp,
        r=r,
        r_idx=r_idx,
        v_idx=v_idx,
        F_idx=F_idx,
        n_ground=n_ground,
        n_excited=n_excited,
        n_states=n_states,
        m=m,
        update_p=update_p,
        n_scatters=n_scatters
    )

function initialize_schrodinger_prob(
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
    add_terms_dψ
)

    # get some initial constants
    n_states = length(ω0s)
    n_g = find_n_g(d)
    n_excited = n_states - n_g
    n_freqs = length(ωs)
    n_ks = 6

    # set the integer type for the simulation
    intT = get_int_type(sim_type)

    denom = sim_type((beam_radius*k)^2/2)

    eiω0ts = zeros(Complex{sim_type}, n_states)

    # amplitudes of electric fields
    as = zeros(Complex{sim_type}, n_ks, n_freqs)
    as = StructArray(MMatrix{n_ks, n_freqs}(as))

    # define polarization array
    ϵs = zeros(Complex{sim_type}, n_ks, n_freqs, 3)
    for i ∈ eachindex(pols)
        pol = pols[i]
        ϵs[1,i,:] .= rotate_pol(pol, x̂)
        ϵs[2,i,:] .= rotate_pol(pol, ŷ)
        ϵs[3,i,:] .= rotate_pol(flip(pol), ẑ)
        ϵs[4,i,:] .= rotate_pol(pol, -x̂)
        ϵs[5,i,:] .= rotate_pol(pol, -ŷ)
        ϵs[6,i,:] .= rotate_pol(flip(pol), -ẑ)
    end
    ϵs = StructArray(MArray{Tuple{n_ks, n_freqs, 3}}(ϵs))

    # arrays related to energies
    ωs = MVector{size(ωs)...}(ωs)
    idxs = MMatrix{size(idxs)...}(idxs)

    # phases of laser frequencies
    ϕs = zeros(sim_type, 3, 3 + n_freqs)
    ϕs = MMatrix{size(ϕs)...}(ϕs)

    # position
    rs = zeros(sim_type,2,3)
    rs = MMatrix{size(rs)...}(rs)

    #
    kEs = zeros(Complex{sim_type}, k_dirs, 3)
    kEs = StructArray(MMatrix{size(kEs)...}(kEs))

    # total electric field
    E_total = zeros(Complex{sim_type}, 3)
    E_total = MVector{size(E_total)...}(E_total)
    E_total = StructArray(E_total)

    # note that we take the negative to ensure that the Hamiltonian is -d⋅E
    d_ge = sim_type.(real.(-d[1:n_g,(n_g+1):n_states,:]))
    d_ge = MArray{Tuple{size(d_ge)...}}(d_ge)
    d_eg = permutedims(d_ge,(2,1,3))

    d = MArray{Tuple{n_states,n_states,3}}(sim_type.(real.(d)))

    d_exp = zeros(Complex{sim_type},3)
    d_exp = MVector{size(d_exp)...}(d_exp)
    d_exp = StructArray(d_exp)

    d_exp_split = zeros(Complex{sim_type},3,2)
    d_exp_split = MMatrix{size(d_exp_split)...}(d_exp_split)
    d_exp_split = StructArray(d_exp_split)

    # force
    F = MVector{3, sim_type}(zeros(3))

    return MutableNamedTuple(

    )

function initialize_prob(
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
    update_params
)

    # get some initial constants
    n_states = length(ω0s)
    n_g = find_n_g(d)
    n_excited = n_states - n_g
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

    ψ = zeros(Complex{sim_type}, n_states)
    ψ = MVector{size(ψ)...}(ψ)
    ψ = StructArray(ψ)

    dψ = deepcopy(ψ)

    ψ_q = deepcopy(ψ)
    ψ_q = MArray{Tuple{size(ψ)...,3}}(zeros(Complex{sim_type},size(ψ)...,3))
    ψ_q = StructArray(ψ_q)

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
    r_idx = 2n_states + n_excited
    v_idx = r_idx + 3
    F_idx = v_idx + 3

    decay_dist = Exponential(one(sim_type))
    last_decay_time = zero(sim_type)

    diffusion_constant = MVector{3,sim_type}(zeros(3))
    add_spontaneous_decay_kick = false

    n_scatters = zero(sim_type)

    u0 = sim_type.([zeros(n_states)..., zeros(n_states)..., zeros(4)..., zeros(3)..., zeros(3)..., zeros(3)..., zeros(3)...])
    u0[1] = 1.0

    p = MutableNamedTuple(
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
        ψ=ψ,
        dψ=dψ,
        ψ_q=ψ_q,
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
        n_excited=n_excited,
        n_states=n_states,
        m=m,
        add_terms_dψ=add_terms_dψ,
        update_params=update_params,
        decay_dist=decay_dist,
        time_to_decay=rand(decay_dist),
        last_decay_time=last_decay_time,
        n_scatters=n_scatters,
        diffusion_constant=diffusion_constant,
        add_spontaneous_decay_kick=add_spontaneous_decay_kick,
        sats=sats
    )

    return p
end