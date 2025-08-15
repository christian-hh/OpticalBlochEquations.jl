import MutableNamedTuples: MutableNamedTuple
import StructArrays: StructArray, StructVector
import StaticArrays: @SVector
import LinearAlgebra: norm, ⋅, adjoint!, diag
import LoopVectorization: @turbo
import LinearAlgebra: det
export schrödinger, obe

"""
Structure for quantum jumps from state `s` to state `s′` with rate `r`.
"""
struct Jump
    s ::Int64
    s′::Int64
    q::Int64
    r ::ComplexF64
end

# Round `val` to the nearest multiple of `prec`
round_to_mult(val, prec) = round.(val ./ prec) .* prec

function round_freq(ω, freq_res)
    ω_min = freq_res
    return round_to_mult(ω, ω_min)
end
export round_freq

function round_vel(v, freq_res)
    v_min = freq_res
    return round_to_mult(v, v_min)
end
export round_vel

function round_v_by_norm(v, digits)
    v_norm = norm(v)
    return v ./ (v_norm / round(v_norm, digits=digits))
end

function round_freqs!(states, fields, freq_res)
    """
    Rounds frequencies of state energies and fields by a common denominator.
    
    freq_res::Float: all frequencies are rounded by this value (in units of Γ)
    """
    for i in eachindex(fields)
        fields.ω[i] = round_freq(fields.ω[i], freq_res)
    end
    for i in eachindex(states)
        states.E[i] = round_freq(states.E[i], freq_res)
    end
    return nothing
end

function round_params(p)
    round_freqs!(p.states, p.fields, p.freq_res)
    p.v .= round_vel(p.v, p.freq_res)
end
export round_params

function obe(ρ0, particle, states, fields, d, should_round_freqs, include_jumps; 
    sim_params=nothing, extra_data=nothing, λ=1.0, Γ=2π, freq_res=1e-2, update_H_and_∇H=update_H_and_∇H)

    period = 2π / freq_res

    n_states = length(states)
    n_fields = length(fields)

    states = StructArray(states)
    if n_fields > 0
        fields = StructArray(fields)
    end
    
    k = 2π / λ
    # particle.r0 *= 2π #(1 / k)  # `r` is in units of 1/k
    # particle.v /= (Γ / k) # velocity is in units of Γ/k
    # Convert to angular frequencies
    for i ∈ eachindex(fields)
        fields.ω[i] /= Γ
    end
    for i ∈ eachindex(states)
        states.E[i] *= 2π
        states.E[i] /= Γ
    end

    if should_round_freqs
        round_freqs!(states, fields, freq_res)
        particle.v = round_vel(particle.v, freq_res)
    end

    r0 = particle.r0
    r = particle.r
    v = particle.v

    type_complex = ComplexF64

    H = StructArray( zeros(type_complex, n_states, n_states) )

    ω = [state.E for state ∈ states]
    eiωt = StructArray(zeros(type_complex, n_states))

    ρ_soa = StructArray(zeros(ComplexF64, n_states, n_states))
    dρ_soa = deepcopy(ρ_soa)

    # Allocate some temporary arrays
    H₀ = deepcopy(ρ_soa)
    tmp = deepcopy(ρ_soa)

    # Compute cartesian indices to indicate nonzero transition dipole moments in `d`
    # Indices below the diagonal of the Hamiltonian are removed, since those are defined via taking the conjugate
    d_nnz_m = [cart_idx for cart_idx ∈ findall(d[:,:,1] .!= 0) if cart_idx[2] >= cart_idx[1]]
    d_nnz_0 = [cart_idx for cart_idx ∈ findall(d[:,:,2] .!= 0) if cart_idx[2] >= cart_idx[1]]
    d_nnz_p = [cart_idx for cart_idx ∈ findall(d[:,:,3] .!= 0) if cart_idx[2] >= cart_idx[1]]
    d_nnz = [d_nnz_m, d_nnz_0, d_nnz_p]

    # Create jumps corresponding to spontaneous decay
    Js = Array{Jump}(undef, 0)
    ds = [Complex{Float64}[], Complex{Float64}[], Complex{Float64}[]]
    ds_state1 = [Int64[], Int64[], Int64[]]
    ds_state2 = [Int64[], Int64[], Int64[]]
    for s′ in eachindex(states), s in s′:n_states, q in qs
        dme = d[s′, s, q+2]
        # println(dme, " ", s′, " ", s)
        # println(states[s′].E, " ", states[s].E)
        if abs(dme) > 1e-10 && (states[s′].E < states[s].E) # only energy-allowed jumps are generated
        # if (states[s′].E < states[s].E) # only energy-allowed jumps are generated
            push!(ds_state1[q+2], s)
            push!(ds_state2[q+2], s′)
            push!(ds[q+2], dme)
            J = Jump(s, s′, q, dme)
            # J = Jump(s, s′, q, norm(dme)) # is this needed? probably not!
            push!(Js, J)
        end
    end
    ds = [StructArray(ds[1]), StructArray(ds[2]), StructArray(ds[3])]

    # The last 3 indices are for the integrated force
    populations = diag(ρ0)
    ρ0_vec = [[ρ0[i] for i ∈ eachindex(ρ0)]; populations; zeros(3)]

    force_last_period = SVector(0.0, 0.0, 0.0)

    # Some additional arrays to hold information about fields
    # fields_ϵ = [SVector{3, ComplexF64}(0.,0.,0.) for _ ∈ eachindex(fields)]
    # fields_kr = zeros(length(fields))
    # fields_re = zeros(length(fields))
    # fields_im = zeros(length(fields))

    E = @SVector Complex{Float64}[0,0,0]
    E_k = [@SVector Complex{Float64}[0,0,0] for _ ∈ 1:3]

    p = MutableNamedTuple(
        H=H, particle=particle, ρ0=ρ0, ρ0_vec=ρ0_vec, ρ_soa=ρ_soa, dρ_soa=dρ_soa, Js=Js, eiωt=eiωt, ω=ω,
        states=states, fields=fields, r0=r0, r=r, v=v, Γ=Γ, tmp=tmp, λ=λ,
        period=period, k=k, freq_res=freq_res, H₀=H₀,
        force_last_period=force_last_period, populations=populations,
        d=d, d_nnz=d_nnz,
        E=E, E_k=E_k,
        ds=ds, ds_state1=ds_state1, ds_state2=ds_state2,
        sim_params=sim_params,
        extra_data=extra_data,
        update_H_and_∇H=update_H_and_∇H)

    return p
end
export obe


function obe_old(ρ0, particle, states, fields, d, should_round_freqs, include_jumps; 
    λ=1.0, Γ=2π, freq_res=1e-2, extra_p=extra_p)

    period = 2π / freq_res

    n_states = length(states)
    # n_fields = length(fields)

    states = StructArray(states)
    # if n_fields > 0
        fields = StructArray(fields)
    # end
    
    k = 2π / λ
    # particle.r0 *= 2π #(1 / k)  # `r` is in units of 1/k
    # particle.v /= (Γ / k) # velocity is in units of Γ/k
    # Convert to angular frequencies
    for i ∈ eachindex(fields)
        fields.ω[i] /= Γ
    end
    for i ∈ eachindex(states)
        states.E[i] *= 2π
        states.E[i] /= Γ
    end

    if should_round_freqs
        round_freqs!(states, fields, freq_res)
        particle.v = round_vel(particle.v, freq_res)
    end

    r0 = particle.r0
    r = particle.r
    v = particle.v

    type_complex = ComplexF64

    H = StructArray( zeros(type_complex, n_states, n_states) )

    ω = [state.E for state ∈ states]
    eiωt = StructArray(zeros(type_complex, n_states))

    ρ_soa = StructArray(zeros(ComplexF64, n_states, n_states))
    dρ_soa = deepcopy(ρ_soa)

    # Allocate some temporary arrays
    H₀ = deepcopy(ρ_soa)
    tmp = deepcopy(ρ_soa)

    # Compute cartesian indices to indicate nonzero transition dipole moments in `d`
    # Indices below the diagonal of the Hamiltonian are removed, since those are defined via taking the conjugate
    d_nnz_m = [cart_idx for cart_idx ∈ findall(d[:,:,1] .!= 0) if cart_idx[2] >= cart_idx[1]]
    d_nnz_0 = [cart_idx for cart_idx ∈ findall(d[:,:,2] .!= 0) if cart_idx[2] >= cart_idx[1]]
    d_nnz_p = [cart_idx for cart_idx ∈ findall(d[:,:,3] .!= 0) if cart_idx[2] >= cart_idx[1]]
    d_nnz = [d_nnz_m, d_nnz_0, d_nnz_p]

    # Create jumps corresponding to spontaneous decay
    Js = Array{Jump}(undef, 0)
    ds = [Complex{Float64}[], Complex{Float64}[], Complex{Float64}[]]
    ds_state1 = [Int64[], Int64[], Int64[]]
    ds_state2 = [Int64[], Int64[], Int64[]]
    for s′ in eachindex(states), s in s′:n_states, q in qs
        dme = d[s′, s, q+2]
        # println(dme, " ", s′, " ", s)
        # println(states[s′].E, " ", states[s].E)
        if abs(dme) > 1e-10 && (states[s′].E < states[s].E) # only energy-allowed jumps are generated
        # if (states[s′].E < states[s].E) # only energy-allowed jumps are generated
            push!(ds_state1[q+2], s)
            push!(ds_state2[q+2], s′)
            push!(ds[q+2], dme)
            J = Jump(s, s′, q, dme)
            push!(Js, J)
        end
    end
    ds = [StructArray(ds[1]), StructArray(ds[2]), StructArray(ds[3])]

    # The last 3 indices are for the integrated force
    populations = diag(ρ0)
    ρ0_vec = [[ρ0[i] for i ∈ eachindex(ρ0)]; populations; zeros(3)]

    force_last_period = SVector(0.0, 0.0, 0.0)

    # Some additional arrays to hold information about fields
    # fields_ϵ = [SVector{3, ComplexF64}(0.,0.,0.) for _ ∈ eachindex(fields)]
    # fields_kr = zeros(length(fields))
    # fields_re = zeros(length(fields))
    # fields_im = zeros(length(fields))

    E = @SVector Complex{Float64}[0,0,0]
    E_k = [@SVector Complex{Float64}[0,0,0] for _ ∈ 1:3]

    p = MutableNamedTuple(
        H=H, particle=particle, ρ0=ρ0, ρ0_vec=ρ0_vec, ρ_soa=ρ_soa, dρ_soa=dρ_soa, Js=Js, eiωt=eiωt, ω=ω,
        states=states, fields=fields, r0=r0, r=r, v=v, Γ=Γ, tmp=tmp, λ=λ,
        period=period, k=k, freq_res=freq_res, H₀=H₀,
        force_last_period=force_last_period, populations=populations,
        d=d, d_nnz=d_nnz,
        E=E, E_k=E_k,
        ds=ds, ds_state1=ds_state1, ds_state2=ds_state2,
        extra_p=extra_p
        )

    return p
end
export obe_old

function soa_to_base!(ρ::Array{<:Complex}, ρ_soa::StructArray{<:Complex})
    @inbounds for i in eachindex(ρ_soa)
        ρ[i] = ρ_soa.re[i] + im * ρ_soa.im[i]
    end
    return nothing
end

function base_to_soa!(ρ::Array{<:Complex}, ρ_soa::StructArray{<:Complex})
    @inbounds for i in eachindex(ρ_soa)
        ρ_soa.re[i] = real(ρ[i])
        ρ_soa.im[i] = imag(ρ[i])
    end
    return nothing
end
export base_to_soa!

function base_to_soa_vec!(ρ::Array{<:Complex}, ρ_soa::StructArray{<:Complex})
    @inbounds for i in eachindex(ρ, ρ_soa)
        ρ_soa.re[i] = real(ρ[i])
        ρ_soa.im[i] = imag(ρ[i])
    end
    return nothing
end

function update_eiωt!(eiωt::StructArray{<:Complex}, ω::Array{<:Real}, τ::Real)
    @turbo for i ∈ eachindex(ω)
        eiωt.im[i], eiωt.re[i] = sincos( ω[i] * τ )
    end
    return nothing
end

function update_eiωt_nocomplex!(eiωt_re::Array{<:Real}, eiωt_im::Array{<:Real}, ω::Array{<:Real}, τ::Real)
    @turbo for i ∈ eachindex(ω)
        eiωt_im[i], eiωt_re[i] = sincos( ω[i] * τ )
    end
    return nothing
end

"""
    Apply the transformation (A)_(ij) -> (A)_(ij) * exp(-iω_it) * exp(+iω_jt)
"""
function Heisenberg_obes!(A::StructArray{<:Complex}, eiωt::StructArray{<:Complex}, im_prefactor)
    @inbounds for j ∈ 1:size(A, 2)
        jre = eiωt.re[j]
        jim = eiωt.im[j]
        for i ∈ 1:size(A, 1)
            ire = eiωt.re[i]
            iim = -eiωt.im[i]
            cisim = im_prefactor * (iim * jre + ire * jim)
            cisre = ire * jre - iim * jim
            Are_i = A.re[i,j]
            Aim_i = A.im[i,j]
            A.re[i,j] = Are_i * cisre - Aim_i * cisim
            A.im[i,j] = Are_i * cisim + Aim_i * cisre
        end
    end
    return nothing
end
export Heisenberg!

"""
    Apply the transformation (A)_(ij) -> (A)_(ij) * exp(-iω_it) * exp(+iω_jt)
"""
function Heisenberg!(A::StructArray{<:Complex}, eiωt::StructArray{<:Complex})
    @inbounds for j ∈ axes(A, 2)
        jre = eiωt.re[j]
        jim = -eiωt.im[j] # negative sign to get exp(-iω_it) rather than exp(iω_it)
        for i ∈ axes(A, 1)
            ire = eiωt.re[i]
            # iim = -eiωt.im[i] # added negative sign on 10/20/23 --> shouldn't the negative sign be here?
            iim = eiωt.im[i]
            cisim = iim * jre + ire * jim
            cisre = ire * jre - iim * jim
            Are_i = A.re[i,j]
            Aim_i = A.im[i,j]
            A.re[i,j] = Are_i * cisre - Aim_i * cisim
            A.im[i,j] = Are_i * cisim + Aim_i * cisre
        end
    end
    return nothing
end
export Heisenberg!

function Heisenberg_turbo!(A::StructArray{<:Complex}, eiωt::StructArray{<:Complex}, sign=+1)
    @turbo for j ∈ axes(A, 2)
        jre = eiωt.re[j]
        jim = sign * (-eiωt.im[j]) # negative sign to get exp(-iω_it) rather than exp(iω_it)
        for i ∈ axes(A, 1)
            ire = eiωt.re[i]
            iim = sign * eiωt.im[i]
            cisim = iim * jre + ire * jim
            cisre = ire * jre - iim * jim
            Are_i = A.re[i,j]
            Aim_i = A.im[i,j]
            A.re[i,j] = Are_i * cisre - Aim_i * cisim
            A.im[i,j] = Are_i * cisim + Aim_i * cisre
        end
    end
    return nothing
end
export Heisenberg_turbo!

function Heisenberg_turbo_state!(v::StructArray{<:Complex}, eiωt::StructArray{<:Complex}, sign=+1)
    @turbo for i ∈ eachindex(eiωt)
        eiωt_re = eiωt.re[i]
        eiωt_im = sign * eiωt.im[i]
        vre_i = v.re[i]
        vim_i = v.im[i]
        v.re[i] = vre_i * eiωt_re - vim_i * eiωt_im
        v.im[i] = vre_i * eiωt_im + vim_i * eiωt_re
    end
    return nothing
end
export Heisenberg_turbo_state!

function Heisenberg_nocomplex!(A_re::Array{<:Real}, A_im::Array{<:Real}, eiωt_re::Array{<:Real}, eiωt_im::Array{<:Real})
    @inbounds for j ∈ 1:size(A_re, 2)
        jre = eiωt_re[j]
        jim = -eiωt_im[j] # negative sign to get exp(-iω_it) rather than exp(iω_it)
        for i ∈ 1:size(A_re, 1)
            ire = eiωt_re[i]
            iim = eiωt_im[i]
            cisim = iim * jre + ire * jim
            cisre = ire * jre - iim * jim
            Are_i = A_re[i,j]
            Aim_i = A_im[i,j]
            A_re[i,j] = Are_i * cisre - Aim_i * cisim
            A_im[i,j] = Are_i * cisim + Aim_i * cisre
        end
    end
    return nothing
end

# function Heisenberg!(A::StructArray{<:Complex}, eiωt::StructArray{<:Complex}, im_factor=1)
#     @inbounds for j ∈ 1:size(A, 2)
#         jre = eiωt.re[j]
#         jim = eiωt.im[j]
#         for i ∈ 1:size(A, 1)
#             ire = eiωt.re[i]
#             iim = eiωt.im[i]
#             cisim = im_factor * (iim * jre - ire * jim)
#             cisre = ire * jre + iim * jim
#             Are_i = A.re[i,j]
#             Aim_i = A.im[i,j]
#             A.re[i,j] = Are_i * cisre - Aim_i * cisim
#             A.im[i,j] = Are_i * cisim + Aim_i * cisre
#         end
#     end
#     return nothing
# end
# export Heisenberg!

function Heisenberg_ψ!(ψ::StructArray{<:Complex}, eiωt::StructArray{<:Complex}, im_factor=1)
    @turbo for i ∈ 1:size(ψ, 1)
        cisre = eiωt.re[i]
        cisim = im_factor * eiωt.im[i]
        ψ_i_re = ψ.re[i]
        ψ_i_im = ψ.im[i]
        ψ.re[i] = ψ_i_re * cisre - ψ_i_im * cisim
        ψ.im[i] = ψ_i_re * cisim + ψ_i_im * cisre
    end
    return nothing
end
export Heisenberg_ψ!

function im_commutator!(C, A, B, tmp)
    # Multiply C = A * B -- (H * ρ)
    # Add adjoint C = A * B + B† * A†
    # Multiply by "i"
    @turbo for i ∈ 1:size(A,1), j ∈ 1:size(B,2)
        Cre = zero(eltype(C.re))
        Cim = zero(eltype(C.im))
        for k ∈ 1:size(A,2)
            Aik_re = A.re[i,k]
            Aik_im = A.im[i,k]
            Bkj_re = B.re[k,j]
            Bkj_im = B.im[k,j]
            # Multiply by -1
            Cre -= Aik_re * Bkj_re - Aik_im * Bkj_im
            Cim -= Aik_re * Bkj_im + Aik_im * Bkj_re
        end
        C.re[i,j] = Cre
        C.im[i,j] = Cim # not sure about this minus sign...
    end

    # Add adjoint, which is equivalent to ρ * H, then compute [H, ρ]
    adjoint!(tmp, C)
    C_add_A!(C, tmp, -1)
    mul_by_im!(C)

end
export im_commutator!

"""
    ρ!(dρ, ρ, p, τ)

Evaluates the change in the density matrix `dρ` given the current density matrix `ρ` for parameters `p` and time `τ`.
"""
function ρ!(dρ, ρ, p, τ)

    @unpack H, H₀, E, E_k, dρ_soa, ρ_soa, tmp, Js, eiωt, ω, fields, ds, ds_state1, ds_state2, Γ, r, r0, v, Js = p

    r .= r0 .+ v * τ

    base_to_soa!(ρ, ρ_soa)

    # Update the Hamiltonian according to the new time τ
    p.update_H_and_∇H(H₀, p, r, τ)
    
    update_H_obes!(p, τ, r, H₀, fields, H, E_k, ds, ds_state1, ds_state2, Js)

    # Apply a transformation to go to the Heisenberg picture
    update_eiωt!(eiωt, ω, τ)
    Heisenberg_obes!(ρ_soa, eiωt, +1.0)

    # Compute coherent evolution terms
    im_commutator!(dρ_soa, H, ρ_soa, tmp)

    # Add the terms ∑ᵢJᵢρJᵢ†
    @inbounds for i ∈ eachindex(Js)
        J = Js[i]
        dρ_soa[J.s′, J.s′] += norm(J.r)^2 * ρ_soa[J.s, J.s]
        @inbounds for j ∈ (i+1):length(Js)
            J′ = Js[j]
            if J.q == J′.q
                val = conj(J.r) * J′.r * ρ_soa[J.s, J′.s]
                dρ_soa[J.s′, J′.s′] += val
                dρ_soa[J′.s′, J.s′] += conj(val)
            end
        end
    end

    # The left-hand side also needs to be transformed into the Heisenberg picture
    # To do this, we require the transpose of the `ω` matrix
    Heisenberg_obes!(dρ_soa, eiωt, -1.0)
    soa_to_base!(dρ, dρ_soa)

    n = length(ρ_soa)
    for i ∈ axes(ρ_soa, 1)
        dρ[n+i] = ρ_soa[i,i]
    end
    f = force_noupdate(E_k, ds, ds_state1, ds_state2, ρ_soa)
    dρ[end-2:end] = f

    return nothing
end
export ρ!

function ρ_updated!(dρ, ρ, p, t)

    @unpack H, H₀, E, E_k, dρ_soa, ρ_soa, tmp, Js, eiωt, ω, fields, ds, ds_state1, ds_state2, Γ, r, r0, v, Js = p

    r .= r0 .+ v * t

    base_to_soa!(ρ, ρ_soa)

    # Update the Hamiltonian according to the new time τ
    update_H_obes!(p, t, r, H₀, fields, H, E_k, ds, ds_state1, ds_state2, Js)

    # add additional part of Hamiltonian that is not light-molecule part
    ∇H = p.update_H_and_∇H(H₀, p, r, t) # in schrodinger picutre
    # Heisenberg!(H₀, eiωt) # convert to interaction picture
    
    @turbo for i ∈ eachindex(H)
        H.re[i] += H₀.re[i]
        H.im[i] += H₀.im[i]
    end

    # Apply a transformation to go to the interaction picture
    update_eiωt!(eiωt, ω, t)
    Heisenberg_obes!(ρ_soa, eiωt, +1.0)

    # Compute coherent evolution terms
    im_commutator!(dρ_soa, H, ρ_soa, tmp)

    # Add the terms ∑ᵢJᵢρJᵢ†
    # We assume jumps take the form Jᵢ = sqrt(Γ)|g⟩⟨e| such that JᵢρJᵢ† = Γ^2|g⟩⟨g|ρₑₑ
    @inbounds for i ∈ eachindex(Js)
        J = Js[i]
        dρ_soa[J.s′, J.s′] += norm(J.r)^2 * ρ_soa[J.s, J.s]
        @inbounds for j ∈ (i+1):length(Js)
            J′ = Js[j]
            if J.q == J′.q
                val = conj(J.r) * J′.r * ρ_soa[J.s, J′.s]
                dρ_soa[J.s′, J′.s′] += val
                dρ_soa[J′.s′, J.s′] += conj(val)
            end
        end
    end

    # The left-hand side also needs to be transformed into the Heisenberg picture
    # To do this, we require the transpose of the `ω` matrix
    Heisenberg_obes!(dρ_soa, eiωt, -1.0)
    soa_to_base!(dρ, dρ_soa)

    n = length(ρ_soa)
    for i ∈ axes(ρ_soa, 1)
        dρ[n+i] = ρ_soa[i,i]
    end

    f = force_noupdate(E_k, ds, ds_state1, ds_state2, ρ_soa)
    H₀_expectation = H_exp(H₀, ρ_soa)
    f += ∇H .* (-H₀_expectation) # force due to conservative potential

    dρ[end-2:end] = f

    return nothing
end
export ρ_updated!

function mat_to_vec!(ρ, ρ_vec)
    @turbo for i in eachindex(ρ)
        ρ_vec[i] = ρ[i]
    end
    return nothing
end
export mat_to_vec!

function mat_to_vec_minus1!(ρ, ρ_vec)
    @turbo for i in 1:(length(ρ)-1)
        ρ_vec[i] = ρ[i]
    end
    return nothing
end
export mat_to_vec_minus1!

function C_copy_AplusB!(C::Array{<:Real}, A::Array{<:Real}, B::Array{<:Real}, α=1, β=1)
    @turbo for i in eachindex(A, B, C)
        C[i] = α * A[i] + β * B[i]
    end
end

function C_add_AplusB!(C::Array{<:Real}, A::Array{<:Real}, B::Array{<:Real}, α=1, β=1)
    @turbo for i in eachindex(A, B, C)
        C[i] += α * A[i] + β * B[i]
    end
end

function C_add_AplusB!(C::StructArray{<:Complex}, A::StructArray{<:Complex}, B::StructArray{<:Complex},
    α=1, β=1)
    @turbo for i in eachindex(A, B, C)
        C.re[i] += α * A.re[i] + β * B.re[i]
        C.im[i] += α * A.im[i] + β * B.im[i]
    end
end

function C_add_A!(C::StructArray{<:Complex}, A::StructArray{<:Complex}, factor=1)
    @turbo for i in eachindex(C, A)
        C.re[i] += factor * A.re[i]
        C.im[i] += factor * A.im[i]
    end
end

function C_copy_T₁T₂!(C::StructArray{<:Complex}, T1::Array{<:Real}, T2::Array{<:Real})
    @turbo for i ∈ eachindex(C)
        C.re[i] = T2[i] - T1[i]
        C.im[i] = T2[i] + T1[i]
    end
end

function mul_by_im!(C::StructArray{<:Complex})
    @turbo for i ∈ eachindex(C)
        a = C.re[i]
        C.re[i] = -C.im[i]
        C.im[i] = a
    end
end

function mul_by_im_minus!(C::StructArray{<:Complex})
    @turbo for i ∈ eachindex(C)
        a = C.re[i]
        C.re[i] = C.im[i]
        C.im[i] = -a
    end
end

function mul_diagonal!(C::StructArray{<:Complex}, A::StructArray{<:Complex}, B::Array{<:Real})
    """
    Computes A × B, where B is diagonal and real.
    """
    @turbo for j in axes(A,2)
        d = B[j]
        for i in axes(A,1)
            C.re[i,j] = d * A.re[i,j]
            C.im[i,j] = d * A.im[i,j]
        end
    end
end
export mul_diagonal!

function update_T₁T₂!(T1::Array{<:Real}, T2::Array{<:Real}, A::StructArray{<:Complex}, B::StructArray{<:Complex})
    @turbo for i ∈ 1:size(A,1), j ∈ 1:size(B,2)
        C1 = zero(eltype(A))
        C2 = zero(eltype(A))
        for k ∈ 1:size(A,2)
            C1 += A.re[i,k] * B.re[k,j]
            C2 += A.im[i,k] * B.im[k,j]
        end
        T1[i,j] = C1
        T2[i,j] = C2
    end
end

function D(cosβ, sinβ, α, γ) 
    return [
        (1/2)*(1 + cosβ)*exp(im*α)*exp(im*γ) -(1/√2)*sinβ*exp(im*α) (1/2)*(1 - cosβ)*exp(im*α)*exp(-im*γ);
        (1/√2)*sinβ*exp(im*γ) cosβ -(1/√2)*sinβ*exp(-im*γ);
        (1/2)*(1 - cosβ)*exp(-im*α)*exp(im*γ) (1/√2)*sinβ*exp(-im*α) (1/2)*(1 + cosβ)*exp(-im*α)*exp(-im*γ)
    ]
end
export D

using LinearAlgebra: cross

function rotate_pol(pol, k)::SVector{3, Complex{Float64}}
    # Rotates polarization `pol` onto the quantization axis `k`
    k /= norm(k)

    # Find the axis-angle rotation corresponding to the rotation
    x,y,z = cross(k, ẑ)
    θ = acos(k ⋅ ẑ)

    # k = SVector(k[1], k[2], k[3]) # phase convention for the rotation

    # x = -x
    y = -y
    # z = -z
    # println(k)
    # println(x)
    # println(y)
    # println(z)
    # println(θ)

    if θ ≈ 0
        α = 0.
        β = 0.
        γ = 0.
    elseif θ ≈ π
        α = 0.
        β = π
        γ = 0.
    else

        A33 = (1 - cos(θ)) * z^2 + cos(θ)
        A31 = (1 - cos(θ)) * z * x - y * sin(θ)
        A32 = (1 - cos(θ)) * z * y + x * sin(θ)
        A13 = (1 - cos(θ)) * x * z + y * sin(θ)
        A23 = (1 - cos(θ)) * y * z - x * sin(θ)

        α = atan(A23, A13)
        β = atan(sqrt(1 - A33^2), A33)
        γ = atan(A32, -A31)

        A12 = x * y * (1 - cos(θ)) - z * sin(θ)
        A21 = y * x * (1 - cos(θ)) + z * sin(θ)
        A22 = cos(θ) + y^2 * (1 - cos(θ))

        # α = atan(A12, A32)
        # β = acos(A22)
        # γ = atan(A21, -A23)

        A11 = cos(θ) + x^2 * (1 - cos(θ))

        # α = π + atan(-A23, A33)
        # β = -asin(A13)
        # γ = atan(-A12, A11)

        # α = atan(-A31, A11)
        # β = asin(A21)
        # γ = atan(-A23, A22)

        # α = atan(A21, A11)
        # β = asin(-A31)
        # γ = atan(A32, A33)

        # α = atan(A13, -A23)
        # β = acos(A33)
        # γ = atan(A31, A32)

    end

    # k = k / norm(k)
    # β = acos(k[3])
    # cosβ = cos(β)
    # sinβ = sqrt(1 - cosβ^2)
    # α = 4atan(k[1])
    # if abs(cosβ) < 1
    #     γ = atan(k[2], k[1])
    # else
    #     γ = 0.0
    # end
    # # γ = π
    # if γ > 0
    #     pol *= im # phase convention
    # end

    # print(α, " ", β, " ", γ)

    # print(γ)
    # return inv(D(cosβ, sinβ, α, γ)) * pol

    return inv(D(cos(β), sin(β), α, γ)) * pol
end
export rotate_pol

# Multiplication using `@turbo` from LoopVectorization
"""
    C <- AB
"""
function mul_turbo!(C, A, B)
    @turbo for m ∈ 1:size(A,1), n ∈ 1:size(B,2)
        Cmn_re = zero(eltype(C))
        Cmn_im = zero(eltype(C))
        for k ∈ 1:size(A,2)
            A_mk_re = A.re[m,k]
            A_mk_im = A.im[m,k]
            B_kn_re = B.re[k,n]
            B_kn_im = B.im[k,n]
            Cmn_re += A_mk_re * B_kn_re - A_mk_im * B_kn_im
            Cmn_im += A_mk_re * B_kn_im + A_mk_im * B_kn_re
        end
        C.re[m,n] = Cmn_re
        C.im[m,n] = Cmn_im
    end
end
export mul_turbo!

# Same as above but with an extra -i factor
function mul_turbo_im_minus!(C, A, B)
    @turbo for m ∈ 1:size(A,1), n ∈ 1:size(B,2)
        Cmn_re = zero(eltype(C))
        Cmn_im = zero(eltype(C))
        for k ∈ 1:size(A,2)
            A_mk_re = A.re[m,k]
            A_mk_im = A.im[m,k]
            B_kn_re = B.im[k,n] # B -> B * (-i)
            B_kn_im = -B.re[k,n] # B -> B * (-i)
            Cmn_re += A_mk_re * B_kn_re - A_mk_im * B_kn_im
            Cmn_im += A_mk_re * B_kn_im + A_mk_im * B_kn_re
        end
        C.re[m,n] = Cmn_re
        C.im[m,n] = Cmn_im
    end
end
export mul_turbo!

function mul_turbo_transposeA!(C, A, B)
    @turbo for m ∈ 1:size(A,1), n ∈ 1:size(B,2)
        Cmn_re = zero(eltype(C))
        Cmn_im = zero(eltype(C))
        for k ∈ 1:size(A,2)
            A_mk_re = A.re[k,m]
            A_mk_im = -A.im[k,m]
            B_kn_re = B.re[k,n]
            B_kn_im = B.im[k,n]
            Cmn_re += A_mk_re * B_kn_re - A_mk_im * B_kn_im
            Cmn_im += A_mk_re * B_kn_im + A_mk_im * B_kn_re
        end
        C.re[m,n] = Cmn_re
        C.im[m,n] = Cmn_im
    end
end
export mul_turbo!

function mul_turbo_conjA!(C, A, B)
    @turbo for m ∈ 1:size(A,1), n ∈ 1:size(B,2)
        Cmn_re = zero(eltype(C))
        Cmn_im = zero(eltype(C))
        for k ∈ 1:size(A,2)
            A_mk_re = A.re[m,k]
            A_mk_im = -A.im[m,k]
            B_kn_re = B.re[k,n]
            B_kn_im = B.im[k,n]
            Cmn_re += A_mk_re * B_kn_re - A_mk_im * B_kn_im
            Cmn_im += A_mk_re * B_kn_im + A_mk_im * B_kn_re
        end
        C.re[m,n] = Cmn_re
        C.im[m,n] = Cmn_im
    end
end
export mul_turbo_conjA!

function jgemturbo!(y, A, x)
    @turbo for i ∈ eachindex(y)
        yi_re = zero(eltype(y))
        yi_im = zero(eltype(y))
        for j ∈ eachindex(x)
            A_ij_re = A.re[i,j]
            A_ij_im = A.im[i,j]
            x_j_re = x.re[j]
            x_j_im = x.im[j]
            yi_re += A_ij_re * x_j_re - A_ij_im * x_j_im
            yi_im += A_ij_re * x_j_im + A_ij_im * x_j_re
        end
        y.re[i] = yi_re
        y.im[i] = yi_im
    end
    return nothing
end
export jgemturbo!

function mul_turbo_real!(C, A, B)
    @turbo for m ∈ axes(A,1), n ∈ axes(B,2)
        Cmn = zero(eltype(C))
        for k ∈ axes(A,2)
            Cmn += A[m,k] * B[k,n]
        end
        C[m,n] = Cmn
    end
end
export mul_turbo_real!

function mul_turbo_nocomplex!(n_states, dψ, H_re, H_im, ψ)
    @inbounds for i ∈ 1:n_states
        dψ_re_i = zero(eltype(dψ))
        dψ_im_i = zero(eltype(dψ))
        for j ∈ 1:n_states
            ψ_re = ψ[j]
            ψ_im = ψ[j+n_states]
            H_ij_re = H_re[i,j]
            H_ij_im = H_im[i,j]
            dψ_re_i += H_ij_re * ψ_re - H_ij_im * ψ_im
            dψ_im_i += H_ij_re * ψ_im + H_ij_im * ψ_re
        end
        dψ[i] = dψ_re_i
        dψ[i+n_states] = dψ_im_i
    end
    return nothing
end

function mul_turbo_nocomplex_v0!(n_states, dψ, H, ψ)
    H_re = H.re
    H_im = H.im
    @inbounds for i ∈ 1:n_states
        dψ_re_i = zero(eltype(dψ))
        dψ_im_i = zero(eltype(dψ))
        for j ∈ 1:n_states
            ψ_re = ψ[j]
            ψ_im = ψ[j+n_states]
            H_ij_re = H_re[i,j]
            H_ij_im = H_im[i,j]
            dψ_re_i += H_ij_re * ψ_re - H_ij_im * ψ_im
            dψ_im_i += H_ij_re * ψ_im + H_ij_im * ψ_re
        end
        dψ[i] = dψ_re_i
        dψ[i+n_states] = dψ_im_i
    end
    return nothing
end

function mul_turbo_nocomplex_v00!(n_states, C, A, B)
    @inbounds for m ∈ 1:size(A,1), n ∈ 1:size(B,2)
        Cmn_re = 0.0
        Cmn_im = 0.0
        for k ∈ 1:size(A,2)
            A_mk_re = A.re[m,k]
            A_mk_im = A.im[m,k]
            B_kn_re = B[k,n]
            B_kn_im = B[k+n_states,n]
            Cmn_re += A_mk_re * B_kn_re - A_mk_im * B_kn_im
            Cmn_im += A_mk_re * B_kn_im + A_mk_im * B_kn_re
        end
        C[m,n] = Cmn_re
        C[m+n_states,n] = Cmn_im
    end
end

#### GPU FUNCTIONS ####

# function base_to_soa_gpu!(ρ::Array{ComplexF64}, ρ_soa::Array{ComplexF64})
#     for i in eachindex(ρ_soa)
#         ρ_soa.re[i] = ρ[i]
#         ρ_soa.im[i] = ρ[i]
#     end
#     return nothing
# end