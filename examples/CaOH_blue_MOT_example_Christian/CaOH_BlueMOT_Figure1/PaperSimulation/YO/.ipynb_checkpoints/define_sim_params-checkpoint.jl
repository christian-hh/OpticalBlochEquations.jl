### DEFINE OTHER PARAMETERS FOR THE SIMULATION ###

sim_type = Float64

σx_initial = 25e-6
σy_initial = 25e-6
σz_initial = 25e-6
Tx_initial = 50e-6
Ty_initial = 50e-6
Tz_initial = 50e-6

import MutableNamedTuples: MutableNamedTuple
sim_params = MutableNamedTuple(
    Zeeman_Hx = MMatrix{size(Zeeman_x_mat)...}(sim_type.(Zeeman_x_mat)),
    Zeeman_Hy = MMatrix{size(Zeeman_y_mat)...}(sim_type.(Zeeman_y_mat)),
    Zeeman_Hz = MMatrix{size(Zeeman_z_mat)...}(sim_type.(Zeeman_z_mat)),
    
    B_ramp_time = 2e-3 / (1/Γ),
    B_grad_start = 0.,
    B_grad_end = 150.,

    s_ramp_time = 2e-3 / (1/Γ),
    s_factor_start = 2.0,
    s_factor_end = 0.4,

    photon_budget = rand(Geometric(1/1000000)),
    
    x_dist = Normal(0, σx_initial),
    y_dist = Normal(0, σy_initial),
    z_dist = Normal(0, σz_initial),
    
    vx_dist = Normal(0, sqrt(kB*Tx_initial/2m)),
    vy_dist = Normal(0, sqrt(kB*Ty_initial/2m)),
    vz_dist = Normal(0, sqrt(kB*Tz_initial/2m)),

    f_z = StructArray(zeros(Complex{sim_type}, 16, 16)),
    
    total_sat = 0.,
    s1_ratio = s1,
    s2_ratio = s2,
    s3_ratio = s3,
    s4_ratio = s4,
    
    Bx = 0.,
    By = 0.,
    Bz = 0.,

    dt_diffusion = 1e-7 / (1/Γ)
)