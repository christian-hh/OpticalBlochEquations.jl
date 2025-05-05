### DEFINE OTHER PARAMETERS FOR THE SIMULATION ###

sim_type = Float64

σx_initial = 50e-6
σy_initial = 50e-6
σz_initial = 50e-6
Tx_initial = 100e-6
Ty_initial = 100e-6
Tz_initial = 100e-6

w = 20e-6
P = 10
I0_trap = 2P / (π * w^2)

using Serialization
H_ODT_matrix = deserialize("H_ODT_matrix.jl")

import MutableNamedTuples: MutableNamedTuple
sim_params = MutableNamedTuple(
    trap_scalar = 2π * 0.03 * 2I0_trap / (ε0 * c),
    H_ODT_matrix = MMatrix{size(H_ODT_matrix)...}(sim_type.(H_ODT_matrix)),
    
    x_dist = Normal(0, σx_initial),
    y_dist = Normal(0, σy_initial),
    z_dist = Normal(0, σz_initial),
    
    vx_dist = Normal(0, sqrt(kB*Tx_initial/2m)),
    vy_dist = Normal(0, sqrt(kB*Ty_initial/2m)),
    vz_dist = Normal(0, sqrt(kB*Tz_initial/2m)),

    f_z = StructArray(zeros(Complex{sim_type}, 16, 16)),

    dt_diffusion = 1e-7 / (1/Γ)
)