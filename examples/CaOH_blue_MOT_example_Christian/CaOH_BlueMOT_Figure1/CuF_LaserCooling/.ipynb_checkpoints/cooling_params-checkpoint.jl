# # DEFINE STATES #
# energy_offset = (2π / Γ) * energy(states[13])
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +0.0
δ1 = +10.00
δ2 = +0.00

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)

f1 = energy(states[end]) - energy(states[1]) + Δ1
f2 = energy(states[end]) - energy(states[10]) + Δ2

freqs = [f1, f2] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
beam_radius = 5e-3
Isat = π*h*c*Γ/(3λ^3)
P = 1000e-3 # mW
I = 2P / (π * beam_radius^2)

total_sat = I / Isat
s1 = 0.5total_sat
s2 = 0.5total_sat

sats = [s1, s2]

# DEFINE POLARIZATIONS #
pols = [σ⁻, σ⁺]

# DEFINE FUNCTION TO UPDATE PARAMETERS DURING SIMULATION #
@everywhere function update_p!(p, r, t)
    return nothing
end

@everywhere function update_p_diffusion!(p, r, t)
    return nothing
end
