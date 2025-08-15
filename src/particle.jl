export Particle, sample_direction,update_initial_position!, update_initial_velocity!

@with_kw mutable struct Particle
    r0::MVector{3, Float64} = MVector(0.0, 0.0, 0.0)
    r::MVector{3, Float64}  = MVector(0.0, 0.0, 0.0)
    v::MVector{3, Float64}  = MVector(0.0, 0.0, 0.0)
end

function sample_direction(r=1.0)
    θ = 2π * rand()
    z = rand() * 2 - 1
    return (r * sqrt(1 - z^2) * cos(θ), r * sqrt(1 - z^2) * sin(θ), r * z)
end

function cdf_maxwell_boltzmann(v)
    return nothing
end

function update_initial_velocity!(prob, v)
    prob.u0[2prob.p.n_states + prob.p.n_e + 4] = v[1]
    prob.u0[2prob.p.n_states + prob.p.n_e + 5] = v[2]
    prob.u0[2prob.p.n_states + prob.p.n_e + 6] = v[3]
    return nothing
end

function update_initial_position!(prob, r)
    prob.u0[2prob.p.n_states + prob.p.n_e + 1] = r[1]
    prob.u0[2prob.p.n_states + prob.p.n_e + 2] = r[2]
    prob.u0[2prob.p.n_states + prob.p.n_e + 3] = r[3]
    return nothing
end

import OrdinaryDiffEqCore: ODEIntegrator
function update_initial_velocity!(integrator::ODEIntegrator, v)
    integrator.u[2integrator.p.n_states + integrator.p.n_e + 4] = v[1]
    integrator.u[2integrator.p.n_states + integrator.p.n_e + 5] = v[2]
    integrator.u[2integrator.p.n_states + integrator.p.n_e + 6] = v[3]
    return nothing
end

function update_initial_position!(integrator::ODEIntegrator, r)
    integrator.u[2integrator.p.n_states + integrator.p.n_e + 1] = r[1]
    integrator.u[2integrator.p.n_states + integrator.p.n_e + 2] = r[2]
    integrator.u[2integrator.p.n_states + integrator.p.n_e + 3] = r[3]
    return nothing
end