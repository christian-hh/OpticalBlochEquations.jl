QN_bounds = (
    label = "X",
    S = 1/2, 
    I = 1/2, 
    Λ = 0, 
    N = 1
)
basis = order_basis_by_m(enumerate_states(HundsCaseB_LinearMolecule, QN_bounds))

H_operator = :(
    BX * Rotation + 
    DX * RotationDistortion + 
    γX * SpinRotation + 
    bFX * Hyperfine_IS + 
    cX * (Hyperfine_Dipolar/3)
)

parameters = @params begin
    BX = 11633.6 * 1e6
    DX = 4.8078e-7 * 1e6
    γX = -9.2254 * 1e6
    bFX = (-762.976 - 9.412) * 1e6 # bFx = b + c/3
    cX = -28.236 * 1e6
end

X_state_ham = Hamiltonian(basis=basis, operator=H_operator, parameters=parameters)
evaluate!(X_state_ham)
QuantumStates.solve!(X_state_ham)
;