# -*- coding: utf-8 -*-
Zeeman_x(state, state′) = (Zeeman(state, state′,-1) - Zeeman(state, state′,1))/sqrt(2)
Zeeman_y(state, state′) = im*(Zeeman(state, state′,-1) + Zeeman(state, state′,1))/sqrt(2)
Zeeman_z(state, state′) = Zeeman(state, state′, 0)

Zeeman_x_mat = StructArray(operator_to_matrix_zero_padding2(Zeeman_x, X_states, A_states[1:4]) .* (2π*gS*_μB/Γ))
Zeeman_y_mat = StructArray(operator_to_matrix_zero_padding2(Zeeman_y, X_states, A_states[1:4]) .* (2π*gS*_μB/Γ))
Zeeman_z_mat = StructArray(operator_to_matrix_zero_padding2(Zeeman_z, X_states, A_states[1:4]) .* (2π*gS*_μB/Γ))
