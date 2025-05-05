### SIMULATION PARAMETERS: 3-FREQUENCY 1+2 BLUE MOT  ###

# DEFINE STATES #
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +9.0
δ1 = +0.00
δ2 = +3.00
δ3 = -0.60
δ4 = +0.60

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)
Δ3 = 1e6 * (detuning + δ3)
Δ4 = 1e6 * (detuning + δ4)

f1 = energy(states[end]) - energy(states[1]) + Δ1
f2 = energy(states[end]) - energy(states[5]) + Δ2
f3 = energy(states[end]) - energy(states[10]) + Δ3
f4 = energy(states[end]) - energy(states[10]) + Δ4

freqs = [f1, f2, f3, f4] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
beam_radius = 5e-3
Isat = π*h*c*Γ/(3λ^3)
P = @with_unit 5 "mW"
I = 2P / (π * beam_radius^2)

total_sat = I / Isat

# size is ~33-35 um
s1 = 0.22total_sat
s2 = 0.58total_sat
s3 = 0.10total_sat
s4 = 0.10total_sat

sats = [s1, s2, s3, s4]

# DEFINE POLARIZATIONS #
pols = [σ⁺, σ⁻, σ⁺, σ⁻]

# DEFINE FUNCTION TO UPDATE PARAMETERS DURING SIMULATION #
@everywhere function update_p!(p, r, t)

    # set ramped scale factor for saturation parameters
    s_scalar = min(t / p.sim_params.s_ramp_time, 1.0)
    s_factor = p.sim_params.total_sat * (p.sim_params.s_factor_start + (p.sim_params.s_factor_end - p.sim_params.s_factor_start) * s_scalar)
    p.sats[1] = p.sim_params.s1_ratio * s_factor
    p.sats[2] = p.sim_params.s2_ratio * s_factor
    p.sats[3] = p.sim_params.s3_ratio * s_factor
    p.sats[4] = p.sim_params.s4_ratio * s_factor

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
    p.sats[3] = p.sim_params.s3_ratio * s_factor
    p.sats[4] = p.sim_params.s4_ratio * s_factor
    return nothing
end