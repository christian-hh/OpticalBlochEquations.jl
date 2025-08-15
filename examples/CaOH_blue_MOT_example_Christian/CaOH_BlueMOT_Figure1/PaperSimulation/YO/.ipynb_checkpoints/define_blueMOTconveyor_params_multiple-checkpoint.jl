### SIMULATION PARAMETERS: 3-FREQUENCY 1+2 BLUE MOT  ###

# DEFINE STATES #
energies = energy.(states) .* (2π / Γ)

# DEFINE FREQUENCIES #
detuning = +5.0
δ1 = +0.00
δ2 = +0.00
δ3 = -0.60
δ4 = +0.60

Δ1 = 1e6 * (detuning + δ1)
Δ2 = 1e6 * (detuning + δ2)
Δ3 = 1e6 * (detuning + δ3)
Δ4 = 1e6 * (detuning + δ4)

f1 = energy(states[end]) - energy(states[14]) + Δ1
f2 = energy(states[end]) - energy(states[13]) + Δ2
f3 = energy(states[end]) - energy(states[5]) + Δ3
f4 = energy(states[end]) - energy(states[5]) + Δ4

f_mw = energy(states[5]) - energy(states[1])

ks = [1, 1/48000]

freqs = [[f1, f2, f3, f4], [f_mw]] .* (2π / Γ)

# DEFINE SATURATION INTENSITIES #
beam_radius = 5e-3
Isat = π*h*c*Γ/(3λ^3)
P = @with_unit 2 "mW"
Ip = 2P / (π * beam_radius^2)

total_sat = Ip / Isat
# s1 = 0.15total_sat
# s2 = 0.00total_sat
# s3 = 0.25total_sat
# s4 = 0.25total_sat

s1 = 0.28total_sat
s2 = 0.16total_sat
s3 = 0.28total_sat
s4 = 0.28total_sat

s_mw = 0

sats = [[s1, s2, s3, s4], [s_mw]]

# DEFINE POLARIZATIONS #
pols = [[σ⁺, σ⁻, σ⁺, σ⁻], [σ⁰]]

# DEFINE FUNCTION TO UPDATE PARAMETERS DURING SIMULATION #
@everywhere function update_p!(p, r, t)

    # set ramped scale factor for saturation parameters
    s_scalar = min(t / p.sim_params.s_ramp_time, 1.0)
    s_factor = p.sim_params.total_sat * (p.sim_params.s_factor_start + (p.sim_params.s_factor_end - p.sim_params.s_factor_start) * s_scalar)
    p.sats[1][1] = p.sim_params.s1_ratio * s_factor
    p.sats[1][2] = p.sim_params.s2_ratio * s_factor
    p.sats[1][3] = p.sim_params.s3_ratio * s_factor
    p.sats[1][4] = p.sim_params.s4_ratio * s_factor

    # set ramped scale factor for B field
    B_scalar = min(t / p.sim_params.B_ramp_time, 1.0)
    B_grad = p.sim_params.B_grad * B_scalar
    p.sim_params.Bx = +r[1] * B_grad * 1e2 / k / 2
    p.sim_params.By = +r[2] * B_grad * 1e2 / k / 2
    p.sim_params.Bz = -r[3] * B_grad * 1e2 / k
    
    return nothing
end

@everywhere function update_p_diffusion!(p, r, t)
    s_factor = p.sim_params.total_sat * p.sim_params.s_factor_end
    p.sats[1][1] = p.sim_params.s1_ratio * s_factor
    p.sats[1][2] = p.sim_params.s2_ratio * s_factor
    p.sats[1][3] = p.sim_params.s3_ratio * s_factor
    p.sats[1][4] = p.sim_params.s4_ratio * s_factor
    return nothing
end

import MutableNamedTuples: MutableNamedTuple
sim_params = MutableNamedTuple(
    Zeeman_Hx = MMatrix{size(Zeeman_x_mat)...}(Zeeman_x_mat),
    Zeeman_Hy = MMatrix{size(Zeeman_y_mat)...}(Zeeman_y_mat),
    Zeeman_Hz = MMatrix{size(Zeeman_z_mat)...}(Zeeman_z_mat),
    
    B_ramp_time = 4e-3 / (1/Γ),
    B_grad = 74.,

    s_ramp_time = 1e-9 / (1/Γ),
    s_factor_start = 1.0,
    s_factor_end = 1.0,

    photon_budget = rand(Geometric(1/1000000)),
    
    x_dist = Normal(0, 500e-6),
    y_dist = Normal(0, 500e-6),
    z_dist = Normal(0, 500e-6),
    
    vx_dist = Normal(0, sqrt(kB*35e-6/2m)),
    vy_dist = Normal(0, sqrt(kB*35e-6/2m)),
    vz_dist = Normal(0, sqrt(kB*35e-6/2m)),

    # f_z = StructArray(zeros(ComplexF64, 16, 16)),
    f_z = StructArray(zeros(ComplexF64, 20, 20)),
    
    total_sat = total_sat,
    s1_ratio = s1,
    s2_ratio = s2,
    s3_ratio = s3,
    s4_ratio = s4,
    
    Bx = 0.,
    By = 0.,
    Bz = 0.,

    dt_diffusion = 1e-7 / (1/Γ)
)