# DEFINE STATES #
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +0.44 * 4.8
δ1 = -0.264
δ2 = +0.000

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)

f1 = energy(states[end]) - energy(states[5]) + Δ1
f2 = energy(states[end]) - energy(states[14]) + Δ2
f_mw = energy(states[5]) - energy(states[1])

# freqs = [f1, f2] .* (2π / Γ)
# freqs = [[f1, f2]] .* (2π / Γ)
freqs = [[f1, f2], [f_mw]] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
beam_radius = 5e-3
Isat = π*h*c*Γ/(3λ^3)

s1 = 1.93
s2 = 1.93

s_mw = 10.

# sats = [s1, s2]
# sats = [[s1, s2]]
sats = [[s1, s2], [s_mw]]

# DEFINE POLARIZATIONS #
# pols = [σ⁺, σ⁻]
# pols = [[σ⁺, σ⁻]]
pols = [[σ⁺, σ⁻], [σ⁰]]
# pols = [[σ⁺, σ⁻], [σ⁺]]

ks = [1, f_mw/f1]

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
    
    x_dist = Normal(0, 800e-6),
    y_dist = Normal(0, 800e-6),
    z_dist = Normal(0, 800e-6),
    
    vx_dist = Normal(0, sqrt(kB*50e-6/2m)),
    vy_dist = Normal(0, sqrt(kB*50e-6/2m)),
    vz_dist = Normal(0, sqrt(kB*50e-6/2m)),

    f_z = StructArray(zeros(ComplexF64, 20, 20)),

    B_grad = 4,
    
    Bx = 0.,
    By = 0.,
    Bz = 0.,

    dt_diffusion = 1e-7 / (1/Γ)
)