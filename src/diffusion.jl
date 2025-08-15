export compute_diffusion

"""
    Compute the z component of the force and store the result in `p.f_z`.
"""
function compute_fz(prob)

    p = prob.p
    
    f_z = p.sim_params.f_z
    kEs = p.kEs[1]
    d_ge = p.d
    n_g = p.n_g

    k = 3
    @turbo for j ∈ axes(d_ge,2)
        for i ∈ axes(d_ge,1)
            f_z_ij_re = zero(eltype(f_z.re))
            f_z_ij_im = zero(eltype(f_z.im))
            for q ∈ 1:3
                d = d_ge[i,j,q]
                E_kq_re = -(kEs.im[k,q] - kEs.im[k+3,q]) # multiply by -i
                E_kq_im = +(kEs.re[k,q] - kEs.re[k+3,q])
                val_x_re = d * E_kq_re
                val_x_im = d * E_kq_im
                f_z_ij_re += val_x_re
                f_z_ij_im += val_x_im
            end
            f_z.re[i,j+n_g] = f_z_ij_re
            f_z.im[i,j+n_g] = f_z_ij_im
            f_z.re[j+n_g,i] = f_z_ij_re
            f_z.im[j+n_g,i] = -f_z_ij_im
        end
    end

    return nothing
end
export compute_fz

"""
    Compute the two-time correlation for `prob` using the Monte Carlo wavefunction method.

    `prob_func`: function that updates `prob` for each sample.
    `n_avgs`: number of samples.
    `t_end`: time 
    `τ_total`: total time to compute diffusion out to
    `n_times`: total number of times to save

    This implementation follows reference [ref].
"""
function compute_diffusion(prob, prob_func, n_avgs, t_end, τ_total, n_times, channel)
    
    n_states = prob.p.n_states
    n_excited = prob.p.n_e
    F_idx = prob.p.F_idx
    coord_idx = 3

    # data arrays
    Cs = zeros(Float64, n_times)
    fτ_fts = zeros(Float64, n_times)

    Cs_integrated = zeros(Float64, n_avgs)
    fτ_fts_integrated = zeros(Float64, n_avgs)
    
    # initialize with equal population between all states
    ψ₀ = StructArray(zeros(ComplexF64, n_states))
    ψ₀[1:n_states] .= 1.0
    ψ₀ ./= norm(ψ₀)

    ϕ = deepcopy(ψ₀)

    for i ∈ 1:n_avgs

        put!(channel, true)

        # set times to simulate
        t_start = 0.0
        t_end′  = t_end + 1e-6 * rand()
        t_span  = (t_start, t_end′)

        τ_start = t_end′
        τ_end   = τ_start + τ_total
        τ_span  = (τ_start, τ_end)
        τ_times = range(τ_start, τ_end, n_times)

        # round all times
        Γ = prob.p.Γ
        t_span  = round.(t_span ./ (1/Γ), digits=9)
        τ_span  = round.(τ_span ./ (1/Γ), digits=9)
        τ_times = round.(τ_times ./ (1/Γ), digits=9)

        dτ = τ_times[2] - τ_times[1]
        τ_span = (τ_span[1], τ_span[2] + dτ)

        # create the new problem
        prob.u0 .= 0.
        update_u!(prob.u0, ψ₀, n_states)
        prob_func(prob)
        prob.p.last_decay_time = 0.
        prob.p.time_to_decay = rand(prob.p.decay_dist)

        # solve the problem
        sol_ϕ = solve(prob, tspan=t_span)
        ut = sol_ϕ.u[end]

        last_decay_time = sol_ϕ.prob.p.last_decay_time
        time_to_decay = sol_ϕ.prob.p.time_to_decay
        
        ϕ.re .= ut[1:n_states]
        ϕ.im .= ut[(n_states+1):2n_states]

        compute_fz(sol_ϕ.prob)
        Heisenberg!(sol_ϕ.prob.p.sim_params.f_z, sol_ϕ.prob.p.eiω0ts)
        f = sol_ϕ.prob.p.sim_params.f_z

        χm = ϕ .- f*ϕ
        χp = ϕ .+ f*ϕ
        μm = norm(χm)
        μp = norm(χp)
        χm ./= μm
        χp ./= μp

        ft = real(ϕ' * f * ϕ)
        fχm = real(χm' * f * χm)
        fχp = real(χp' * f * χp)

        # solve the problem up to time `t`
        prob.u0 .= ut
        prob.p.last_decay_time = last_decay_time
        prob.p.time_to_decay = time_to_decay
        prob.u0[F_idx + coord_idx + 3] = 0.
        sol_ϕτ = solve(prob, tspan=τ_span, saveat=τ_times)

        # solve χm
        prob.u0 .= ut
        prob.p.time_to_decay = rand(prob.p.decay_dist)
        for i ∈ 1:n_excited
            prob.u0[2n_states+i] = 0.
        end
        prob.p.last_decay_time = last_decay_time
        update_u!(prob.u0, χm, n_states)
        prob.u0[F_idx + coord_idx] = fχm
        prob.u0[F_idx + coord_idx + 3] = 0.
        sol_χm = solve(prob, tspan=τ_span, saveat=τ_times)

        # solve χp
        prob.u0 .= ut
        prob.p.time_to_decay = rand(prob.p.decay_dist)
        for i ∈ 1:n_excited
            prob.u0[2n_states+i] = 0.
        end
        prob.p.last_decay_time = last_decay_time
        update_u!(prob.u0, χp, n_states)
        prob.u0[F_idx + coord_idx] = fχp
        prob.u0[F_idx + coord_idx + 3] = 0.
        sol_χp = solve(prob, tspan=τ_span, saveat=τ_times)
        
        # using the integrated version of the force
        cm = sol_χm.u[end][F_idx + coord_idx + 3]
        cp = sol_χp.u[end][F_idx + coord_idx + 3]
        C = (1/4) * (μp^2 * cp - μm^2 * cm)

        fτ = sol_ϕτ.u[end][F_idx + coord_idx + 3]

        Cs_integrated[i] = C
        fτ_fts_integrated[i] = ft * fτ

        for j ∈ 1:n_times

            cm = sol_χm.u[j][F_idx + coord_idx + 3]
            cp = sol_χp.u[j][F_idx + coord_idx + 3]
    
            C = (1/4) * (μp^2 * cp - μm^2 * cm)
            Cs[j] += C
            
            fτ = sol_ϕτ.u[j][F_idx + coord_idx + 3]
            fτ_ft = fτ * ft
            fτ_fts[j] += fτ_ft

        end
        
    end
    Cs ./= n_avgs
    fτ_fts ./= n_avgs
    
    return Cs, fτ_fts, Cs_integrated, fτ_fts_integrated
end
export compute_diffusion

# import DifferentialEquations: init, set_t!, step!
# function compute_diffusion(prob, prob_func, n_avgs, t_total, τ_total, n_times, channel)
    
#     Γ = prob.p.Γ

#     integrator = init(prob)

#     n_states = prob.p.n_states
#     n_excited = prob.p.n_e
#     F_idx = prob.p.F_idx
#     coord_idx = 3

#     # data arrays
#     Cs = zeros(Float64, n_times)
#     fτ_fts = zeros(Float64, n_times)

#     Cs_integrated = zeros(Float64, n_avgs)
#     fτ_fts_integrated = zeros(Float64, n_avgs)

#     # initialize with equal population between all states
#     ψ₀ = StructArray(zeros(ComplexF64, n_states))
#     ψ₀[1:n_states] .= 1.0
#     ψ₀ ./= norm(ψ₀)

#     ϕ = deepcopy(ψ₀)
#     ut = deepcopy(integrator.u)

#     fτs = zeros(n_times)
#     cps = zeros(n_times)
#     cms = zeros(n_times)

#     for i ∈ 1:n_avgs
        
#         put!(channel, true)
        
#         # set times to simulate
#         t_end′ = (t_total + 1e-6 * rand()) / (1/Γ)
#         dτ = τ_total / n_times / (1/Γ)

#         # create the new problem, solve it up to `t`
#         set_t!(integrator, 0)
#         update_u!(integrator.u, ψ₀, n_states)
#         prob_func(integrator)
#         integrator.p.last_decay_time = 0.
#         integrator.p.time_to_decay = rand(integrator.p.decay_dist)
#         step!(integrator, t_end′, true)
#         ut .= integrator.u

#         last_decay_time = integrator.p.last_decay_time
#         time_to_decay = integrator.p.time_to_decay
        
#         ϕ.re .= ut[1:n_states]
#         ϕ.im .= ut[(n_states+1):2n_states]

#         compute_fz(integrator)
#         Heisenberg!(integrator.p.sim_params.f_z, integrator.p.eiω0ts)
#         f = integrator.p.sim_params.f_z

#         χm = ϕ .- f*ϕ
#         χp = ϕ .+ f*ϕ
#         μm = norm(χm)
#         μp = norm(χp)
#         χm ./= μm
#         χp ./= μp

#         ft = real(ϕ' * f * ϕ)
#         fχm = real(χm' * f * χm)
#         fχp = real(χp' * f * χp)

#         # solve the problem up to time `τ`
#         integrator.p.last_decay_time = last_decay_time
#         integrator.p.time_to_decay = time_to_decay
#         integrator.u[F_idx + coord_idx + 3] = 0.

#         for j ∈ 1:n_times
#             step!(integrator, dτ, true)
#             fτs[j] = integrator.u[F_idx + coord_idx + 3]
#         end
        
#         # solve χm
#         set_t!(integrator, t_end′)
#         integrator.u .= ut
#         update_u!(integrator.u, χm, n_states)
#         for i ∈ 1:n_excited
#             integrator.u[2n_states+i] = 0.
#         end
#         integrator.u[F_idx + coord_idx] = fχm
#         integrator.u[F_idx + coord_idx + 3] = 0.
#         integrator.p.time_to_decay = rand(integrator.p.decay_dist)
#         integrator.p.last_decay_time = last_decay_time

#         set_t!(integrator, t_end′)
#         for j ∈ 1:n_times
#             step!(integrator, dτ, true)
#             cms[j] = integrator.u[F_idx + coord_idx + 3]
#         end

#         # solve χp
#         set_t!(integrator, t_end′)
#         integrator.u .= ut
#         update_u!(integrator.u, χp, n_states)
#         for i ∈ 1:n_excited
#             integrator.u[2n_states+i] = 0.
#         end
#         integrator.u[F_idx + coord_idx] = fχp
#         integrator.u[F_idx + coord_idx + 3] = 0.
#         integrator.p.time_to_decay = rand(integrator.p.decay_dist)
#         integrator.p.last_decay_time = last_decay_time

#         set_t!(integrator, t_end′)
#         for j ∈ 1:n_times
#             step!(integrator, dτ, true)
#             cps[j] = integrator.u[F_idx + coord_idx + 3]
#         end
        
#         # using the integrated version of the force
#         cm = cms[end]
#         cp = cps[end]
#         C = (1/4) * (μp^2 * cp - μm^2 * cm)

#         fτ = fτs[end]

#         Cs_integrated[i] = C
#         fτ_fts_integrated[i] = ft * fτ

#         for j ∈ 1:n_times

#             cm = cms[j]
#             cp = cps[j]

#             C = (1/4) * (μp^2 * cp - μm^2 * cm)
#             Cs[j] += C
            
#             fτ = fτs[j]
#             fτ_ft = fτ * ft
#             fτ_fts[j] += fτ_ft

#         end
        
#     end

#     Cs ./= n_avgs
#     fτ_fts ./= n_avgs
    
#     return Cs, fτ_fts, Cs_integrated, fτ_fts_integrated
# end
# export compute_diffusion