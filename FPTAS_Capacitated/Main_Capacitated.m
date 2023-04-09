% FPTAS and Optimal iteration are only for clusters with two groups
% Heuristic algorithm can be applied to large cluster

clear;
for sample = 1 : 10
    global I J L epsilon num_group_in_one_cluster cluster_iteration_limitation;
    global utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi ;
    global lower_case_C upper_case_C order_matrix;
    
    rand('seed',100*sum(clock));
    I = 6; %total number of items in each group
    J = 6;  %total number of groups
    L = 1;  %total number of bundles
    utility_v0 = J/2;
    num_group_in_one_cluster = 3;
    lower_case_C = num_group_in_one_cluster; %local capacity
    upper_case_C = 30; %global capacity
    cluster_iteration_limitation = (I+1)^num_group_in_one_cluster;
    DivisionSampleNumber = 100;
    % load('capacitated_(7,7,2,3)_(3).mat');
%     revenue_matrix_r = sort(rand(I,J)*0.5,1,'descend'); % each column is non-increasing
%     revenue_matrix_r(1:floor(I/3),1:floor(I/3)) = sort(1.5 + rand(floor(I/3),floor(I/3))*0.5,1,'descend');
%     revenue_matrix_r(floor(I/3)+1:I,1:floor(I/3)) = sort(rand(I - floor(I/3),floor(I/3))*1.5,1,'descend');
%     revenue_matrix_r = sort(revenue_matrix_r,1,'descend');
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
    str = 'capacitated_(6,6,1,3)_('+string(sample)+').mat';
%     save(str,"interation_para_phi","revenue_matrix_r","utility_matrix_v");
    load(str);
    revenue_matrix_r = sort(revenue_matrix_r,1,'descend');

    fprintf('There are %d groups, %d items in each group and %d clusters, each cluster has %d group. \n',J,I,L,num_group_in_one_cluster);
    fprintf('The local capacity is %d, the global capacity is %d \n',lower_case_C,upper_case_C);
    fprintf('------------------------%d \n',sample)
    order_matrix = StaticMNL(I, J, revenue_matrix_r, utility_matrix_v);
    
    %Optimal
    % [true_max_profit, true_optimal_x] = Optimal_capacitated ()
    %Benchmark
    [best_PR_revenue, best_PR_X] = Pure_RO_capacitated ()
    [BestSampleRevenue_MVMNL, BestSampleX] = GroupwiseMNL_capacitated (DivisionSampleNumber)

%     epsilon = 0.8
    % [FPTAS_best_revenue2, best_K2, FPTAS_best_X2] = FPTAS_capacitated () %maximum_revenue_index is the K which achieves the best revenue
%     [Heuristic_best_revenue_GS2, maximum_revenue_index_GS2, Heuristic_best_X_GS2] = Heuristic_Golden_Section_capacitated ()
%     (Heuristic_best_revenue_GS2 - best_PR_revenue)/best_PR_revenue
end