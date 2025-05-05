@everywhere function scan_nothing(prob, scan_value)
    return nothing
end

function compute_trajectories_with_diffusion(
        prob, prob_func, prob_diffusion, prob_func_diffusion, 
        n_trajectories1, n_trajectories2, n_trajectories_diffusion, 
        n_times, diffusion_t_end, diffusion_τ_total, scan_func=scan_nothing, scan_values=[0]
    )
    """
    1) Compute 
    2) 
    3)
    4)
    5)
    """

    prob_diffusion.p.add_spontaneous_decay_kick = false

    # 1) Trajectory
    prob.p.diffusion_constant[1] = 0.
    prob.p.diffusion_constant[2] = 0.
    prob.p.diffusion_constant[3] = 0.
    prob.p.add_spontaneous_decay_kick = true
    all_sols_no_diffusion = distributed_solve(n_trajectories1, prob, prob_func, scan_func, scan_values)
    
    σxs = σx_fit.(all_sols_no_diffusion)
    σys = σy_fit.(all_sols_no_diffusion)
    σzs = σz_fit.(all_sols_no_diffusion)
    Txs = Tx_fit.(all_sols_no_diffusion)
    Tys = Ty_fit.(all_sols_no_diffusion)
    Tzs = Tz_fit.(all_sols_no_diffusion)
    
    # 2) Diffusion
    scan_values_with_σ_and_T = zip(scan_values, zip(σxs, σys, σzs, Txs, Tys, Tzs))
    # prob_diffusion.p.add_spontaneous_decay_kick = true

    (diffusions, diffusion_errors, diffusions_over_time1) = distributed_compute_diffusion(
        prob_diffusion, prob_func_diffusion, n_trajectories_diffusion, diffusion_t_end, diffusion_τ_total, n_times, scan_func_with_initial_conditions!(scan_func), scan_values_with_σ_and_T
        )

    # 3) Trajectory
    prob.p.add_spontaneous_decay_kick = false
    prob_diffusion.p.add_spontaneous_decay_kick = false

    scan_values_with_diffusion = zip(scan_values, diffusions)
    all_sols_with_diffusion1 = distributed_solve(n_trajectories2, prob, prob_func, scan_func_with_diffusion!(scan_func), scan_values_with_diffusion)
    
    σxs = σx_fit.(all_sols_with_diffusion1)
    σys = σy_fit.(all_sols_with_diffusion1)
    σzs = σz_fit.(all_sols_with_diffusion1)
    Txs = Tx_fit.(all_sols_with_diffusion1)
    Tys = Ty_fit.(all_sols_with_diffusion1)
    Tzs = Tz_fit.(all_sols_with_diffusion1)
    
    # 4) Diffusion
    scan_values_with_σ_and_T = zip(scan_values, zip(σxs, σys, σzs, Txs, Tys, Tzs))
    (diffusions, diffusion_errors, diffusions_over_time2) = distributed_compute_diffusion(
        prob_diffusion, prob_func_diffusion, n_trajectories_diffusion, diffusion_t_end, diffusion_τ_total, n_times, scan_func_with_initial_conditions!(scan_func), scan_values_with_σ_and_T
        )

    # scan_values_with_σ_T_and_diffusion = zip(scan_values, zip(σxs, σys, σzs, Txs, Tys, Tzs), diffusions)
    # (diffusions, diffusion_errors, diffusions_over_time2) = distributed_compute_diffusion(
    #     prob_diffusion, prob_func_diffusion, n_trajectories_diffusion, diffusion_t_end, diffusion_τ_total, n_times, scan_func_with_diffusion_and_initial_conditions!(scan_func), scan_values_with_σ_T_and_diffusion
    #     )

    # 5) Trajectory
    scan_values_with_diffusion = zip(scan_values, diffusions)
    all_sols_with_diffusion2 = distributed_solve(n_trajectories2, prob, prob_func, scan_func_with_diffusion!(scan_func), scan_values_with_diffusion)
    
    if length(scan_values) == 1
        return (all_sols_no_diffusion[1], all_sols_with_diffusion1[1], all_sols_with_diffusion2[1], diffusions[1], diffusion_errors[1], diffusions_over_time1[1], diffusions_over_time2[1])
    else
        return (all_sols_no_diffusion, all_sols_with_diffusion1, all_sols_with_diffusion2, diffusions, diffusion_errors, diffusions_over_time1, diffusions_over_time2)
    end
end

function compute_trajectories_with_diffusion_iterations(
        prob, prob_func, prob_diffusion, prob_func_diffusion, 
        n_trajectories1, n_trajectories2, n_trajectories_diffusion, 
        n_times, diffusion_t_end, diffusion_τ_total, scan_func=scan_nothing, scan_values=[0]
    )
    """
    1) Compute 
    2) 
    3)
    4)
    5)
    """

    # 1) Trajectory
    prob.p.diffusion_constant[1] = 0.
    prob.p.diffusion_constant[2] = 0.
    prob.p.diffusion_constant[3] = 0.
    prob.p.add_spontaneous_decay_kick = true
    all_sols_no_diffusion = distributed_solve(n_trajectories1, prob, prob_func, scan_func, scan_values)
    
    σxs = σx_fit.(all_sols_no_diffusion)
    σys = σy_fit.(all_sols_no_diffusion)
    σzs = σz_fit.(all_sols_no_diffusion)
    Txs = Tx_fit.(all_sols_no_diffusion)
    Tys = Ty_fit.(all_sols_no_diffusion)
    Tzs = Tz_fit.(all_sols_no_diffusion)
    
    # 2) Diffusion
    scan_values_with_σ_and_T = zip(scan_values, zip(σxs, σys, σzs, Txs, Tys, Tzs))
    prob_diffusion.p.add_spontaneous_decay_kick = true

    (diffusions, diffusion_errors, diffusions_over_time1) = distributed_compute_diffusion(
        prob_diffusion, prob_func_diffusion, n_trajectories_diffusion, diffusion_t_end, diffusion_τ_total, n_times, scan_func_with_initial_conditions!(scan_func), scan_values_with_σ_and_T
        )

    # 3) Trajectory
    prob.p.add_spontaneous_decay_kick = false
    prob_diffusion.p.add_spontaneous_decay_kick = false

    scan_values_with_diffusion = zip(scan_values, diffusions)
    all_sols_with_diffusion1 = distributed_solve(n_trajectories2, prob, prob_func, scan_func_with_diffusion!(scan_func), scan_values_with_diffusion)
    
    σxs = σx_fit.(all_sols_with_diffusion1)
    σys = σy_fit.(all_sols_with_diffusion1)
    σzs = σz_fit.(all_sols_with_diffusion1)
    Txs = Tx_fit.(all_sols_with_diffusion1)
    Tys = Ty_fit.(all_sols_with_diffusion1)
    Tzs = Tz_fit.(all_sols_with_diffusion1)
    
    # 4) Diffusion
    # scan_values_with_σ_and_T = zip(scan_values, zip(σxs, σys, σzs, Txs, Tys, Tzs))
    # (diffusions, diffusion_errors, diffusions_over_time) = distributed_compute_diffusion(
    #     prob_diffusion, prob_func_diffusion, n_trajectories_diffusion, diffusion_t_end, diffusion_τ_total, n_times, scan_func_with_initial_conditions!(scan_func), scan_values_with_σ_and_T
    #     )

    scan_values_with_σ_T_and_diffusion = zip(scan_values, zip(σxs, σys, σzs, Txs, Tys, Tzs), diffusions)
    (diffusions, diffusion_errors, diffusions_over_time2) = distributed_compute_diffusion(
        prob_diffusion, prob_func_diffusion, n_trajectories_diffusion, diffusion_t_end, diffusion_τ_total, n_times, scan_func_with_diffusion_and_initial_conditions!(scan_func), scan_values_with_σ_T_and_diffusion
        )

    # 5) Trajectory
    scan_values_with_diffusion = zip(scan_values, diffusions)
    all_sols_with_diffusion2 = distributed_solve(n_trajectories2, prob, prob_func, scan_func_with_diffusion!(scan_func), scan_values_with_diffusion)
    
    if length(scan_values) == 1
        return (all_sols_no_diffusion[1], all_sols_with_diffusion1[1], all_sols_with_diffusion2[1], diffusions[1], diffusion_errors[1], diffusions_over_time1[1], diffusions_over_time2[1])
    else
        return (all_sols_no_diffusion, all_sols_with_diffusion1, all_sols_with_diffusion2, diffusions, diffusion_errors, diffusions_over_time1, diffusions_over_time2)
    end
end

@everywhere function scan_func_with_initial_conditions!(scan_func)
    (prob, scan_value) -> begin
        scan_func(prob, scan_value[1])
        σx, σy, σz, Tx, Ty, Tz = scan_value[2]
        prob.p.sim_params.x_dist = Normal(0, σx)
        prob.p.sim_params.y_dist = Normal(0, σy)
        prob.p.sim_params.z_dist = Normal(0, σz)
        prob.p.sim_params.vx_dist = Normal(0, Tx)
        prob.p.sim_params.vy_dist = Normal(0, Ty)
        prob.p.sim_params.vz_dist = Normal(0, Tz)
        return nothing
    end
end

@everywhere function scan_func_with_diffusion!(scan_func)
    (prob, scan_value) -> begin
        scan_func(prob, scan_value[1])
        diffusion = scan_value[2]
        prob.p.diffusion_constant[1] = diffusion
        prob.p.diffusion_constant[2] = diffusion
        prob.p.diffusion_constant[3] = diffusion
        return nothing
    end
end

@everywhere function scan_func_with_diffusion_and_initial_conditions!(scan_func)
    (prob, scan_value) -> begin
        scan_func(prob, scan_value[1])
        σx, σy, σz, Tx, Ty, Tz = scan_value[2]
        diffusion = scan_value[3]
        prob.p.sim_params.x_dist = Normal(0, σx)
        prob.p.sim_params.y_dist = Normal(0, σy)
        prob.p.sim_params.z_dist = Normal(0, σz)
        prob.p.sim_params.vx_dist = Normal(0, Tx)
        prob.p.sim_params.vy_dist = Normal(0, Ty)
        prob.p.sim_params.vz_dist = Normal(0, Tz)
        prob.p.diffusion_constant[1] = diffusion
        prob.p.diffusion_constant[2] = diffusion
        prob.p.diffusion_constant[3] = diffusion
        return nothing
    end
end

    
    