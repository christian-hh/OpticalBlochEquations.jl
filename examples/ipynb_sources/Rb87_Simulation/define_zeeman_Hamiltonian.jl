zeeman(state, state′, p) = gS * zeeman_BdotS(state, state′, p) + zeeman_BdotL(state, state′, p)

zeeman_x(state, state′) = (zeeman(state, state′,-1) - zeeman(state, state′,1))/sqrt(2)
zeeman_y(state, state′) = im*(zeeman(state, state′,-1) + zeeman(state, state′,1))/sqrt(2)
zeeman_z(state, state′) = zeeman(state, state′, 0)

zeeman_x_mat = StructArray(operator_to_matrix(zeeman_x, states_to_use) .* (2π*_μB/Γ))
zeeman_y_mat = StructArray(operator_to_matrix(zeeman_y, states_to_use) .* (2π*_μB/Γ))
zeeman_z_mat = StructArray(operator_to_matrix(zeeman_z, states_to_use) .* (2π*_μB/Γ))