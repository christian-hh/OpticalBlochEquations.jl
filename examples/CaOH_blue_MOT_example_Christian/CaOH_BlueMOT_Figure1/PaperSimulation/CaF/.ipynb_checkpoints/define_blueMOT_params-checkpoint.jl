# DEFINE STATES #
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +23.8
δ1 = -0.75
δ2 = +0.00

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)

f1 = energy(states[end]) - energy(states[1]) + Δ1
f2 = energy(states[end]) - energy(states[10]) + Δ2

freqs = [f1, f2] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
beam_radius = 6e-3
Isat = π*h*c*Γ/(3λ^3)
P = @with_unit 5.0 "mW"
I = 2P / (π * beam_radius^2)

total_sat = I / Isat
s1 = 1.0total_sat/2.4
s2 = 1.4total_sat/2.4

sats = [s1, s2]

# DEFINE POLARIZATIONS #
pols = [σ⁺, σ⁻]

# DEFINE FUNCTION TO UPDATE PARAMETERS DURING SIMULATION #
@everywhere function update_p!(p, r, t)
    B_grad = p.sim_params.B_grad
    p.sim_params.Bx = +r[1] * B_grad * 1e2 / k / 2
    p.sim_params.By = +r[2] * B_grad * 1e2 / k / 2
    p.sim_params.Bz = -r[3] * B_grad * 1e2 / k
    return nothing
end

@everywhere function update_p_diffusion!(p, r, t)
    return nothing
end

import MutableNamedTuples: MutableNamedTuple
sim_params = MutableNamedTuple(
    Zeeman_Hx = MMatrix{size(Zeeman_x_mat)...}(Zeeman_x_mat),
    Zeeman_Hy = MMatrix{size(Zeeman_y_mat)...}(Zeeman_y_mat),
    Zeeman_Hz = MMatrix{size(Zeeman_z_mat)...}(Zeeman_z_mat),

    photon_budget = rand(Geometric(1/1000000)),
    
    x_dist = Normal(0, 200e-6),
    y_dist = Normal(0, 200e-6),
    z_dist = Normal(0, 200e-6),
    
    vx_dist = Normal(0, sqrt(kB*35e-6/2m)),
    vy_dist = Normal(0, sqrt(kB*35e-6/2m)),
    vz_dist = Normal(0, sqrt(kB*35e-6/2m)),

    f_z = StructArray(zeros(ComplexF64, 16, 16)),

    B_grad = 14.6,
    
    Bx = 0.,
    By = 0.,
    Bz = 0.,

    dt_diffusion = 1e-7 / (1/Γ)
)