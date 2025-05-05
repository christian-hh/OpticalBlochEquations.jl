# Define a Hund's case (b) basis for the Hamiltonian, using states from N=0 to N=3
QN_bounds = (
    label = "X",
    S = 1/2, 
    I = 1/2, 
    Λ = 0, 
    N = 1
)
X_state_basis = order_basis_by_m(enumerate_states(HundsCaseB_LinearMolecule, QN_bounds))

# Define the operator for the X state Hamiltonian of CaOH
X_state_operator = :(
    BX * Rotation + 
    DX * RotationDistortion + 
    γX * SpinRotation + 
    bFX * Hyperfine_IS + 
    cX * (Hyperfine_Dipolar/3)
)

X_state_parameters = QuantumStates.@params begin
    BX = 10023.0841 * 1e6
    DX = 1.154e-2 * 1e6
    γX = 34.7593 * 1e6
    bFX = 2.602 * 1e6
    cX = 2.053 * 1e6
end

X_state_ham = Hamiltonian(basis=X_state_basis, operator=X_state_operator, parameters=X_state_parameters)
evaluate!(X_state_ham)
QuantumStates.solve!(X_state_ham)