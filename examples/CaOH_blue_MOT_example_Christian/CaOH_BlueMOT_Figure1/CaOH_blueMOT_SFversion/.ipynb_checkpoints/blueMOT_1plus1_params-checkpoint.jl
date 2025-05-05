### SIMULATION PARAMETERS: 3-FREQUENCY 1+1 BLUE MOT  ###

# DEFINE STATES #
energy_offset = (2π / Γ) * energy(states[13])
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +7.6
δ1 = +0.00
δ2 = -0.75

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)

f1 = energy(states[end]) - energy(states[1]) + Δ1
f2 = energy(states[end]) - energy(states[10]) + Δ2

freqs = [f1, f2] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
beam_radius = 5e-3
Isat = π*h*c*Γ/(3λ^3)
P = 0.60 * 13.1e-3 # 13.1 mW/1 V at 1.0 V, factor of 0.55 to match scattering rates
I = 2P / (π * beam_radius^2)

total_sat = I / Isat
s1 = 0.5625total_sat
s2 = 0.4375total_sat

sats = [s1, s2]

# DEFINE POLARIZATIONS #
pols = [σ⁻, σ⁺]

# DEFINE FUNCTION TO UPDATE PARAMETERS DURING SIMULATION #
@everywhere function update_p_1plus1!(p, r, t)
    
    s = p.sim_params.total_sat
    s_scalar = min(t / p.sim_params.s_ramp_time, 1.0)
    s_factor = p.sim_params.s_factor_start + (p.sim_params.s_factor_end - p.sim_params.s_factor_start) * s_scalar
    p.sats[1] = 0.5625s * s_factor
    p.sats[2] = 0.4375s * s_factor
    
    # set ramped scale factor for B field
    B_scalar = min(t / p.sim_params.B_ramp_time, 1.0)
    B_grad = p.sim_params.B_grad_start + (p.sim_params.B_grad_end - p.sim_params.B_grad_start) * B_scalar
    p.sim_params.Bx = +r[1] * B_grad * 1e2 / k / 2
    p.sim_params.By = +r[2] * B_grad * 1e2 / k / 2
    p.sim_params.Bz = -r[3] * B_grad * 1e2 / k
    
    return nothing
end

@everywhere function update_p_1plus1_diffusion!(p, r, t)
    s = p.sim_params.total_sat
    s_factor = p.sim_params.s_factor_end
    p.sats[1] = 0.5625s * s_factor
    p.sats[2] = 0.4375s * s_factor
    return nothing
end