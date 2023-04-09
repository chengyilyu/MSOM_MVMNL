function [revenue_GroupwiseMNL] = calculate_revenue_GroupwiseMNL (current_x)
    global I J L num_group_in_one_cluster utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi;
    revenue_GroupwiseMNL = 0;
    for j = 1 : J
        revenue_GroupwiseMNL = revenue_GroupwiseMNL + (current_x(:,j)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))/(utility_v0 + current_x(:,j)'*utility_matrix_v(:,j));
    end
end