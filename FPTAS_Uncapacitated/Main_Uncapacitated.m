% FPTAS and Optimal iteration are only for clusters with two groups
% Heuristic algorithm can be applied to clusters with 3 or 4 groups
clear;
for sample = 1:10
    global I J L epsilon num_group_in_one_cluster cluster_iteration_limitation;
    global utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi;
    
    rand('seed',100*sum(clock));
    I = 21; %total number of items in each group
    J = 21;  %total number of groups
    L = 6;  %total number of bundles
    utility_v0 = J/2;
    num_group_in_one_cluster = 3;
    cluster_iteration_limitation = (I+1)^num_group_in_one_cluster;
    % load('uncapacitated_(3,3,1,2)_(3).mat');
%     revenue_matrix_r = sort(rand(I,J)*0.5,1,'descend'); % each column is non-increasing
%     revenue_matrix_r(1:floor(I/3),1:floor(I/3)) = 1.5 + rand(floor(I/3),floor(I/3))*0.5;
%     revenue_matrix_r(floor(I/3)+1:I,1:floor(I/3)) = rand(I - floor(I/3),floor(I/3))*1.5;
%     revenue_matrix_r = sort(revenue_matrix_r,1,'descend');
% 
% 
%     utility_matrix_v = 1 + rand(I,J);
%     interation_para_phi = ones(J,J);
%     for cluster_num = 1 : L
%         for index_1 = 1 : num_group_in_one_cluster - 1
%             for index_2 = (num_group_in_one_cluster*cluster_num - index_1 + 1) : num_group_in_one_cluster*cluster_num
%                 temp = rand*2;
%                 interation_para_phi(num_group_in_one_cluster*cluster_num - index_1, index_2) = temp;
%                 interation_para_phi(index_2, num_group_in_one_cluster*cluster_num - index_1) = temp;
%             end
%         end
%     end
    str = 'uncapacitated_(21,21,6,3)_('+string(sample)+').mat';
%     save(str,"interation_para_phi","revenue_matrix_r","utility_matrix_v");
    load(str);
    revenue_matrix_r = sort(revenue_matrix_r,1,'descend');

    fprintf('Unapacited case: \n');
    fprintf('There are %d groups, %d items in each group and %d clusters, each cluster has %d group. \n',J,I,L,num_group_in_one_cluster);
    % fprintf('The epsilon is: %.2f \n', epsilon);
    
    fprintf('------------------------%d \n',sample)
    %Optimal
    % [true_max_profit, true_optimal_x] = Optimal_uncapacitated ()
    % benchmarks
%     [best_PR_revenue, best_PR_X] = Pure_RO_uncapacitated ()
    [best_PairwiseMNL_revenue_MVMNL, best_PairwiseMNL_X] = GroupwiseMNL_RO_uncapacitated ()
    
%     epsilon = 0.8
    % FPTAS
    % [FPTAS_best_revenue, best_K, FPTAS_best_X] = FPTAS_uncapacitated () %maximum_revenue_index is the K which achieves the best revenue
    % SGS
%     [Heuristic_best_revenue_inf, Heuristic_maximum_revenue_index_inf, Heuristic_best_X_inf] = Heuristic_Golden_Section_uncapacitated( )
%     (Heuristic_best_revenue_inf - best_PR_revenue)/best_PR_revenue
end