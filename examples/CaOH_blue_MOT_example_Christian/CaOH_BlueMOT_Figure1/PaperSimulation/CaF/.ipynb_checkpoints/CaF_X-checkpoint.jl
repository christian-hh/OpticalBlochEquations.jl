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
    BX = 10303.988 * 1e6
    DX = 0.014060 * 1e6
    γX = 39.65891 * 1e6
    bFX = 122.5569 * 1e6 
    cX = 40.1190 * 1e6
end

X_state_ham = Hamiltonian(basis=basis, operator=H_operator, parameters=parameters)
evaluate!(X_state_ham)
QuantumStates.solve!(X_state_ham)
;