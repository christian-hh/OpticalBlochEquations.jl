# DEFINE STATES #
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +0.44 * 4.8
δ1 = -0.264
δ2 = +0.000

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)

f1 = energy(states[end]) - energy(states[1]) + Δ1
f2 = energy(states[end]) - energy(states[10]) + Δ2

freqs = [f1, f2] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
beam_radius = 5e-3
Isat = π*h*c*Γ/(3λ^3)
P = @with_unit 16 "mW"
I = 2P / (π * beam_radius^2)

total_sat = I / Isat
s1 = 2 #total_sat
s2 = 2 #total_sat

sats = [s1, s2]

# DEFINE POLARIZATIONS #
pols = [σ⁺, σ⁻]

# DEFINE FUNCTION TO UPDATE PARAMETERS DURING SIMULATION #
@everywhere function update_p!(p, r, t)

    # set ramped scale factor for saturation parameters
    s_scalar = min(t / p.sim_params.s_ramp_time, 1.0)
    s_factor = p.sim_params.total_sat * (p.sim_params.s_factor_start + (p.sim_params.s_factor_end - p.sim_params.s_factor_start) * s_scalar)
    p.sats[1] = p.sim_params.s1_ratio * s_factor
    p.sats[2] = p.sim_params.s2_ratio * s_factor

    # set ramped scale factor for B field
    B_scalar = min(t / p.sim_params.B_ramp_time, 1.0)
    B_grad = p.sim_params.B_grad_start + (p.sim_params.B_grad_end - p.sim_params.B_grad_start) * B_scalar
    p.sim_params.Bx = +r[1] * B_grad * 1e2 / k / 2
    p.sim_params.By = +r[2] * B_grad * 1e2 / k / 2
    p.sim_params.Bz = -r[3] * B_grad * 1e2 / k
    
    return nothing
end

@everywhere function update_p_diffusion!(p, r, t)
    s_factor = p.sim_params.total_sat * p.sim_params.s_factor_end
    p.sats[1] = p.sim_params.s1_ratio * s_factor
    p.sats[2] = p.sim_params.s2_ratio * s_factor
end