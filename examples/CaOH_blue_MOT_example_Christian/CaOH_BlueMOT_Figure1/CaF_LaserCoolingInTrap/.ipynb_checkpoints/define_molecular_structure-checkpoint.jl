using Serialization
states = deserialize("states_cooling_sim.jl")
ground_states = states[1:12]
excited_states = states[13:16]

d = zeros(ComplexF64, 16, 16, 3)
d[1:12, 13:16, :] .= tdms_between_states(ground_states, excited_states)
d[13:16, 1:12, :] .= tdms_between_states(excited_states, ground_states)

# Define constants for the laser cooling transition
@everywhere begin
    @consts begin
        λ = 626e-9             # wavelength of transition
        Γ = 2π * 6.4e6         # transition linewidth
        m = @with_unit 57 "u"  # mass of CaF
        k = 2π / λ             # wavenumber of transition
    end
end

@everywhere import LoopVectorization: @turbo
@everywhere function add_terms_dψ!(dψ, ψ, p, r, t)
    @turbo for i ∈ 1:16
        dψ_i_re = zero(eltype(dψ.re))
        dψ_i_im = zero(eltype(dψ.im))
        for j ∈ 1:16
            ψ_i_re = ψ.re[j]
            ψ_i_im = ψ.im[j]
            
            H_re = -p.sim_params.trap_scalar * p.sim_params.H_ODT_matrix[i,j]
            H_im = 0.0
            
            dψ_i_re += ψ_i_re * H_re - ψ_i_im * H_im
            dψ_i_im += ψ_i_re * H_im + ψ_i_im * H_re
            
        end
        dψ.re[i] += dψ_i_im
        dψ.im[i] -= dψ_i_re
    end
    return nothing
end