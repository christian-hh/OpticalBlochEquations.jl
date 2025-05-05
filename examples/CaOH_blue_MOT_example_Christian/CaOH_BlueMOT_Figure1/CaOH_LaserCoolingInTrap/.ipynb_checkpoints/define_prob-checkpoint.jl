# PROBLEM TO CALCULATE TRAJECTORIES #
include("cooling_params.jl")

t_start = 0.0
t_end   = 10e-3
t_span  = (t_start, t_end) ./ (1/Γ)

p = initialize_prob(sim_type, energies, freqs, sats, pols, beam_radius, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p!, add_terms_dψ!)

cb1 = DiscreteCallback(condition_discrete, stochastic_collapse_new!, save_positions=(false,false))
cb2 = DiscreteCallback(terminate_condition, terminate!)
cbs = CallbackSet(cb1, cb2)

# kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=100000000, callback=cbs)
kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=100000000, callback=cb1)
prob = ODEProblem(ψ_fast!, p.u0, sim_type.(t_span), p; kwargs...)

# DIFFUSION PROBLEM #
p_diffusion = initialize_prob(sim_type, energies, freqs, sats, pols, Inf, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p_diffusion!, add_terms_dψ!)
cbs_diffusion = CallbackSet(cb1)

kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=100000000, callback=cbs_diffusion)
prob_diffusion = ODEProblem(ψ_fast_ballistic!, p_diffusion.u0, sim_type.(t_span), p_diffusion; kwargs...)

# PROBLEM WITH PERIODIC CALLBACK FOR DIFFUSION #
@everywhere function diffusion_kick(integrator)
    p = integrator.p
    n_states = p.n_states
    n_excited = p.n_excited
    @inbounds @fastmath for i ∈ 1:3
        kick = sqrt( 2p.diffusion_constant[i] * p.sim_params.dt_diffusion ) / p.m
        integrator.u[2n_states + n_excited + 3 + i] += rand((-1,1)) * kick
    end
    return nothing
end
cb_periodic = PeriodicCallback(diffusion_kick, p.sim_params.dt_diffusion, save_positions=(false,false))
cbs_periodic = CallbackSet(cb1, cb_periodic)

kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=100000000, callback=cbs_periodic)
prob_periodic = ODEProblem(ψ_fast!, p.u0, sim_type.(t_span), p; kwargs...)