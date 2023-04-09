function [bestPolicy_GroupwiseMNL] = FindBestPolicy_GroupwiseMNL()
    global I J L num_group_in_one_cluster utility_v0 revenue_matrix_r utility_matrix_v;
    bestPolicy_GroupwiseMNL = zeros(I,J);
    bestRevenue_GroupwiseMNL = -inf;
    for j = 1 : J
        bestGroupRevenue = -inf;
        for i = 1 : I
            tempX = zeros(I,1);
            tempX(1:i) = 1;
            tempRevenue = (tempX'*(revenue_matrix_r(:,j).*utility_matrix_v(:,j)))/(utility_v0 + tempX'*utility_matrix_v(:,j));
            if tempRevenue > bestGroupRevenue
                bestGroupRevenue = tempRevenue;
                bestPolicy_GroupwiseMNL(:,j) = tempX;
            end
        end
    end
end