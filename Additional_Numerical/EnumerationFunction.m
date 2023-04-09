function [bestPolicy_GroupwiseMNL, bestPolicy_MVMNLwithoutInter, bestPolicy_MVMNLwithInter1, bestPolicy_MVMNLwithInter2, bestPolicy_MVMNLwithInter3,bestPolicy_MVMNLwithInter4,bestPolicy_MVMNLwithInter5] = EnumerationFunction()
    global I J L num_group_in_one_cluster utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi1 interation_para_phi2 interation_para_phi3 interation_para_phi4 interation_para_phi5;
    global bestPolicy_GroupwiseMNL bestPolicy_MVMNLwithoutInter I J bestPolicy_MVMNLwithInter1 bestPolicy_MVMNLwithInter2 bestPolicy_MVMNLwithInter3 bestPolicy_MVMNLwithInter4 bestPolicy_MVMNLwithInter5;
    global bestRevenue_GroupwiseMNL bestRevenue_MVMNLwithoutInter bestRevenue_MVMNLwithInter1 bestRevenue_MVMNLwithInter2 bestRevenue_MVMNLwithInter3 bestRevenue_MVMNLwithInter4 bestRevenue_MVMNLwithInter5;

    bestPolicy_GroupwiseMNL = zeros(I,J);
    bestPolicy_MVMNLwithInter1 = zeros(I,J); 
    bestPolicy_MVMNLwithInter2 = zeros(I,J); 
    bestPolicy_MVMNLwithInter3 = zeros(I,J); 
    bestPolicy_MVMNLwithInter4 = zeros(I,J); 
    bestPolicy_MVMNLwithInter5 = zeros(I,J); 
    bestPolicy_MVMNLwithoutInter = zeros(I,J);
    bestRevenue_GroupwiseMNL = -inf;
    bestRevenue_MVMNLwithInter1 = -inf;
    bestRevenue_MVMNLwithInter2 = -inf;
    bestRevenue_MVMNLwithInter3 = -inf;
    bestRevenue_MVMNLwithInter4 = -inf;
    bestRevenue_MVMNLwithInter5 = -inf;
    bestRevenue_MVMNLwithoutInter = -inf;

    recursion_to_find_the_optimal_policy (1,zeros(I,J));
end


function [ ] = recursion_to_find_the_optimal_policy(group_index, policy_x)
    global I J L num_group_in_one_cluster utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi1 interation_para_phi2 interation_para_phi3 interation_para_phi4 interation_para_phi5;
    global bestPolicy_GroupwiseMNL bestPolicy_MVMNLwithoutInter I J bestPolicy_MVMNLwithInter1 bestPolicy_MVMNLwithInter2 bestPolicy_MVMNLwithInter3 bestPolicy_MVMNLwithInter4 bestPolicy_MVMNLwithInter5;
    global bestRevenue_GroupwiseMNL bestRevenue_MVMNLwithoutInter bestRevenue_MVMNLwithInter1 bestRevenue_MVMNLwithInter2 bestRevenue_MVMNLwithInter3 bestRevenue_MVMNLwithInter4 bestRevenue_MVMNLwithInter5;
    
    if group_index == (J +1) 
        tempRevenue_GroupwiseMNL = calculate_revenue_GroupwiseMNL(policy_x);
        tempRevenue_MVMNLwithoutInter = calculate_revenue_MVMNL_withoutInteraction(policy_x);
        tempRevenue_MVMNLwithInter1 = calculate_revenue_MVMNL_withInteraction(policy_x,interation_para_phi1);
        tempRevenue_MVMNLwithInter2 = calculate_revenue_MVMNL_withInteraction(policy_x,interation_para_phi2);
        tempRevenue_MVMNLwithInter3 = calculate_revenue_MVMNL_withInteraction(policy_x,interation_para_phi3);
        tempRevenue_MVMNLwithInter4 = calculate_revenue_MVMNL_withInteraction(policy_x,interation_para_phi4);
        tempRevenue_MVMNLwithInter5 = calculate_revenue_MVMNL_withInteraction(policy_x,interation_para_phi5);
        if tempRevenue_GroupwiseMNL >= bestRevenue_GroupwiseMNL
            bestRevenue_GroupwiseMNL = tempRevenue_GroupwiseMNL;
            bestPolicy_GroupwiseMNL = policy_x;
        end
        if tempRevenue_MVMNLwithInter1 >= bestRevenue_MVMNLwithInter1
            bestRevenue_MVMNLwithInter1 = tempRevenue_MVMNLwithInter1;
            bestPolicy_MVMNLwithInter1 = policy_x;
        end
        if tempRevenue_MVMNLwithInter2 >= bestRevenue_MVMNLwithInter2
            bestRevenue_MVMNLwithInter2 = tempRevenue_MVMNLwithInter2;
            bestPolicy_MVMNLwithInter2 = policy_x;
        end
        if tempRevenue_MVMNLwithInter3 >= bestRevenue_MVMNLwithInter3
            bestRevenue_MVMNLwithInter3 = tempRevenue_MVMNLwithInter3;
            bestPolicy_MVMNLwithInter3 = policy_x;
        end
        if tempRevenue_MVMNLwithInter4 >= bestRevenue_MVMNLwithInter4
            bestRevenue_MVMNLwithInter4 = tempRevenue_MVMNLwithInter4;
            bestPolicy_MVMNLwithInter4 = policy_x;
        end
        if tempRevenue_MVMNLwithInter5 >= bestRevenue_MVMNLwithInter5
            bestRevenue_MVMNLwithInter5 = tempRevenue_MVMNLwithInter5;
            bestPolicy_MVMNLwithInter5 = policy_x;
        end
        if tempRevenue_MVMNLwithoutInter >= bestRevenue_MVMNLwithoutInter
            bestRevenue_MVMNLwithoutInter = tempRevenue_MVMNLwithoutInter;
            bestPolicy_MVMNLwithoutInter = policy_x;
        end
        return;
    end
    for i = 0 : I
        policy_x(1 : i, group_index) = 1;
        policy_x(i+1:end, group_index) = 0;
        recursion_to_find_the_optimal_policy (group_index+1, policy_x);
    end
end







