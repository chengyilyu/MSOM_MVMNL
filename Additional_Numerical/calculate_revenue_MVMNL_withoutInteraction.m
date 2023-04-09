function [revenue_MVMNL_withoutInteraction] = calculate_revenue_MVMNL_withoutInteraction (current_x)
    global I J L num_group_in_one_cluster utility_v0 revenue_matrix_r utility_matrix_v;
    
    revenue_MVMNL_withoutInteraction = 0;
    Gamma = 1;
    for j = 1 : J
        Gamma = Gamma*(1 + current_x(:,j)'*utility_matrix_v(:,j));
    end
    for j = 1 : J
        revenue_MVMNL_withoutInteraction = revenue_MVMNL_withoutInteraction + (current_x(:,j)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))/(( 1 + current_x(:,j)'*utility_matrix_v(:,j))* ( (utility_v0 - 1)/Gamma + 1));
    end

end