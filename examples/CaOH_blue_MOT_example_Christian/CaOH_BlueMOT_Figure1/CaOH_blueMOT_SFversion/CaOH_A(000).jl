QN_bounds = (
    label = "A",
    S = 1/2,
    I = 1/2,
    Λ = (-1,1),
    J = 1/2
)
A_state_basis = order_basis_by_m(enumerate_states(HundsCaseA_LinearMolecule, QN_bounds))

A_state_operator = :(
    T_A * DiagonalOperator +
    Be_A * Rotation + 
    Aso_A * SpinOrbit +
    q_A * ΛDoubling_q +
    p_A * ΛDoubling_p2q + q_A * (2ΛDoubling_p2q)
)

# Spectroscopic constants for CaOH, A state
A_state_parameters = QuantumStates.@params begin
    T_A = 15998.122 * c * 1e2
    Be_A = 0.3412200 * c * 1e2
    Aso_A = 66.8181 * c * 1e2
    p_A = -0.04287 * c * 1e2
    q_A = -0.3257e-3 * c * 1e2
end

A_state_ham_caseA = Hamiltonian(basis=A_state_basis, operator=A_state_operator, parameters=A_state_parameters)
evaluate!(A_state_ham_caseA)
QuantumStates.solve!(A_state_ham_caseA)

# Convert A state from Hund's case (b) to Hund's case (a)
QN_bounds = (
    label = "A",
    S = 1/2, 
    I = 1/2, 
    Λ = (-1,1), 
    N = 1,
    J = 1/2
)
A_state_caseB_basis = order_basis_by_m(enumerate_states(HundsCaseB_LinearMolecule, QN_bounds))
A_state_ham = convert_basis(A_state_ham_caseA, A_state_caseB_basis)