### SIMULATION PROBLEMS DEFINITION: 1+2 BLUE MOT ###

# PROBLEM TO CALCULATE TRAJECTORIES #
t_start = 0.0
t_end   = 10e-3
t_span  = (t_start, t_end) ./ (1/Γ)

p = initialize_prob(sim_type, energies, freqs, sats, pols, beam_radius, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p!, add_terms_dψ!)

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
kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=1000000000, callback=cbs)
prob = ODEProblem(ψ_fast!, p.u0, sim_type.(t_span), p; kwargs...)

# PROBLEM TO COMPUTE DIFFUSION #
p_diffusion = initialize_prob(sim_type, energies, freqs, sats, pols, Inf, d, m/(ħ*k^2/Γ), Γ, k, sim_params, update_p_diffusion!, add_terms_dψ!)
cbs_diffusion = CallbackSet(cb1, cb_diffusion)
kwargs = (alg=DP5(), reltol=1e-4, abstol=1e-5, saveat=1000, maxiters=1000000000, callback=cbs_diffusion)
prob_diffusion = ODEProblem(ψ_fast_ballistic!, p_diffusion.u0, sim_type.(t_span), p_diffusion; kwargs...)

# set the total saturation
prob.p.sim_params.total_sat = sum(sats)
prob_diffusion.p.sim_params.total_sat = sum(sats)

# set individual saturation ratios
prob.p.sim_params.s1_ratio = s1
prob.p.sim_params.s2_ratio = s2
prob.p.sim_params.s3_ratio = s3
prob.p.sim_params.s4_ratio = s4
prob_diffusion.p.sim_params.s1_ratio = s1
prob_diffusion.p.sim_params.s2_ratio = s2
prob_diffusion.p.sim_params.s3_ratio = s3
prob_diffusion.p.sim_params.s4_ratio = s4

function prob_func_YO!(prob)
    prob.p.n_scatters = 0.
    prob.p.last_decay_time = 0.
    prob.p.time_to_decay = rand(Exponential(1))
    prob.p.sim_params.photon_budget = 1000000
    update_initial_position!(prob, sample_position(prob.p.sim_params) ./ (1/k))
    update_initial_velocity!(prob, sample_velocity(prob.p.sim_params) ./ (Γ/k))
    update_phases!(prob)
    return nothing
end

function prob_func_diffusion_YO!(prob)
    prob_func_YO!(prob)
    prob.p.sim_params.Bx = +x(prob.u0) * prob.p.sim_params.B_grad_end * 1e2 / 2
    prob.p.sim_params.By = +y(prob.u0) * prob.p.sim_params.B_grad_end * 1e2 / 2
    prob.p.sim_params.Bz = -z(prob.u0) * prob.p.sim_params.B_grad_end * 1e2
    return nothing
end