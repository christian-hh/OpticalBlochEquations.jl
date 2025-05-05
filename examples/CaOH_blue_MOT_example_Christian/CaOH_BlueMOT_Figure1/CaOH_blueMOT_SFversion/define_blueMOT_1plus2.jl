### SIMULATION PROBLEMS DEFINITION: 1+2 BLUE MOT ###

# PROBLEM TO CALCULATE TRAJECTORIES #
include("blueMOT_1plus2_params.jl")

t_start = 0.0
t_end   = 10e-3
t_span  = (t_start, t_end) ./ (1/Γ)

p_1plus2 = initialize_prob(sim_type, energies, freqs, sats, pols, beam_radius, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p_1plus2!, add_terms_dψ!)

@everywhere function diffusion_kick(integrator)
    p = integrator.p
    n_states = p.n_states
    n_excited = p.n_excited
    @inbounds @fastmath for i ∈ 1:3
        kick = sqrt(2p.diffusion_constant[i] * p.sim_params.dt_diffusion) / p.m
        integrator.u[2n_states + n_excited + 3 + i] += rand((-1,1)) * kick
    end   
    return nothing
end

diffusion_times = (0:1e-7:t_end) ./ (1/Γ)
cb_diffusion = PresetTimeCallback(diffusion_times, diffusion_kick, save_positions=(false,false))
cb1 = DiscreteCallback(condition_discrete, stochastic_collapse_nokick!, save_positions=(false,false))
cb2 = DiscreteCallback(terminate_condition, terminate!)
cbs = CallbackSet(cb1, cb2, cb_diffusion)
kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=100000000, callback=cbs)
prob_1plus2 = ODEProblem(ψ_fast!, p_1plus2.u0, sim_type.(t_span), p_1plus2; kwargs...)

# PROBLEM TO COMPUTE DIFFUSION #
p_1plus2_diffusion = initialize_prob(sim_type, energies, freqs, sats, pols, Inf, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p_1plus2_diffusion!, add_terms_dψ!)
cbs_diffusion = CallbackSet(cb1, cb_diffusion)
kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=100000000, callback=cbs_diffusion)
prob_1plus2_diffusion = ODEProblem(ψ_fast_ballistic!, p_1plus2_diffusion.u0, sim_type.(t_span), p_1plus2_diffusion; kwargs...)

# set the total saturation
prob_1plus2.p.sim_params.total_sat = sum(sats)
prob_1plus2_diffusion.p.sim_params.total_sat = sum(sats)

# set individual saturation ratios
prob_1plus2.p.sim_params.s1_ratio = 0.37
prob_1plus2.p.sim_params.s2_ratio = 0.28
prob_1plus2.p.sim_params.s3_ratio = 0.35
prob_1plus2_diffusion.p.sim_params.s1_ratio = 0.37
prob_1plus2_diffusion.p.sim_params.s2_ratio = 0.28
prob_1plus2_diffusion.p.sim_params.s3_ratio = 0.35