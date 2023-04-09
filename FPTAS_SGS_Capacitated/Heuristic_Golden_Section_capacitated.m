function [Heuristic_best_revenue, maximum_revenue_index, Heuristic_best_X] = Heuristic_Golden_Section_capacitated ()
    global I J L epsilon num_group_in_one_cluster optimal_V epsilon_prime cluster_iteration_limitation;
    global utility_v0 revenue_matrix_r utility_matrix_v;
    global lower_case_C upper_case_C order_matrix;
    
    L_hat = J - (num_group_in_one_cluster-1)*L; %single group
    epsilon_prime = (1 + epsilon)^(1/(2*L_hat))-1; %(1+epsilon_prime)
    utility_matrix_v_single = utility_matrix_v(: , num_group_in_one_cluster*L+1 : J); %v_{i,j} for single groups

    %%%Heuristic
    %calculate the maximum M
    Heuristic_t_Start = tic;
    gamma_single_1 = 1+ sum(utility_matrix_v_single,1);
    gamma_bundle_1 = zeros(1,L);
    
    if L == 0
        M =  max(ceil(log(gamma_single_1)/log(1+epsilon_prime)));
        Max_K = sum(ceil(log(gamma_single_1)/log(1+epsilon_prime)));
    else
        for bundle_number = 1 : L
            [gamma_bundle_1(bundle_number),~] = calculate_gamma_Lambda_cluster (num_group_in_one_cluster*bundle_number - (num_group_in_one_cluster - 1), ones(I,num_group_in_one_cluster));
        end
        if J == 2*L
            M = max( ceil(log(gamma_bundle_1)/log(1+epsilon_prime)) );
            Max_K = sum(ceil(log(gamma_bundle_1)/log(1+epsilon_prime)));
        else
            M = max( max(ceil(log(gamma_single_1)/log(1+epsilon_prime)) ), max( ceil(log(gamma_bundle_1)/log(1+epsilon_prime)) ));
            Max_K = sum(ceil(log(gamma_single_1)/log(1+epsilon_prime)) ) + sum(ceil(log(gamma_bundle_1)/log(1+epsilon_prime)));
        end
    end

    % main 
    revenue = ones(1,Max_K)*(-inf);
    all_policy_x = {}; % all x for all K

    
    K_lower = L_hat;
    K_upper = Max_K;
    K_1 = 0;
    K_2 = 0;
    temp_K1_old = 0;
    temp_K2_old = 0;
    while(1)
        change_label = false;
        temp_K1 = floor(K_upper - 0.618*(K_upper - K_lower));
        temp_K2 = floor(K_lower + 0.618*(K_upper - K_lower));
        if ((temp_K1 ~= temp_K1_old) || (temp_K2 ~= temp_K2_old))
            K_1 = temp_K1;
            K_2 = temp_K2;
            temp_K1_old = temp_K1;
            temp_K2_old = temp_K2;
            change_label = true;
        end
        
        if revenue(K_1) == -inf
            [revenue(K_1),all_policy_x{K_1}] = calculate_revenue_for_K (K_1);
            temp_1 = 0;
            while (revenue(K_1) == -inf) ||((K_1 - temp_1 < L_hat) && (K_1 + temp_1 > Max_K ))
                temp_1 = temp_1 + 1;
                if K_1 - temp_1 >= L_hat + 1
                    [revenue(K_1 - temp_1),all_policy_x{K_1 - temp_1}] = calculate_revenue_for_K (K_1 - temp_1);
                    if revenue(K_1 - temp_1) ~= -inf 
                        K_1 = K_1 - temp_1;
                        break;
                    end
                end
                if K_1 + temp_1 <= L_hat + 1
                    [revenue(K_1 + temp_1),all_policy_x{K_1 + temp_1}] = calculate_revenue_for_K (K_1 + temp_1);
                    if revenue(K_1 + temp_1) ~= -inf 
                        K_1 = K_1 + temp_1;
                        break;
                    end
                end
            end
        end
        if revenue(K_1 - 1) == -inf
            [revenue(K_1 - 1),all_policy_x{K_1 - 1}] = calculate_revenue_for_K (K_1 - 1);
        end
        if revenue(K_1 + 1) == -inf
            [revenue(K_1 + 1),all_policy_x{K_1 + 1}] = calculate_revenue_for_K (K_1 + 1);
        end
        
        
        if revenue(K_2) == -inf
            [revenue(K_2),all_policy_x{K_2}] = calculate_revenue_for_K (K_2);
            temp_2 = 0;
            while (revenue(K_2) == -inf) ||((K_2 - temp_2 < L_hat) && (K_2 + temp_2 > Max_K ))
                temp_2 = temp_2 + 1;
                if K_2 - temp_2 >= L_hat
                    [revenue(K_2 - temp_2),all_policy_x{K_2 - temp_2}] = calculate_revenue_for_K (K_2 - temp_2);
                    if revenue(K_2 - temp_2) ~= -inf 
                        K_2 = K_2 - temp_2;
                        break;
                    end
                end
                if K_2 + temp_2 <= Max_K
                    [revenue(K_2 + temp_2),all_policy_x{K_2 + temp_2}] = calculate_revenue_for_K (K_2 + temp_2);
                    if revenue(K_2 + temp_2) ~= -inf 
                        K_2 = K_2 + temp_2;
                        break;
                    end
                end
            end
        end
        if revenue(K_2 - 1) == -inf
            [revenue(K_2 - 1),all_policy_x{K_2 - 1}] = calculate_revenue_for_K (K_2 - 1);
        end
        
        if revenue(K_2 + 1) == -inf
            [revenue(K_2 + 1),all_policy_x{K_2 + 1}] = calculate_revenue_for_K (K_2 + 1);
        end
        
        revenue_K1_around_vector = [revenue(K_1 - 1),revenue(K_1),revenue(K_1 + 1)];
        revenue_K2_around_vector = [revenue(K_2 - 1),revenue(K_2),revenue(K_2 + 1)];
               
        revenue_K1_around = mean (revenue_K1_around_vector(revenue_K1_around_vector ~= -inf) );
        revenue_K2_around = mean (revenue_K2_around_vector(revenue_K2_around_vector ~= -inf) );
       
        if revenue_K1_around >= revenue_K2_around
            K_upper = K_2;
        else
            K_lower = K_1;
        end
  
       if ((abs(K_upper - K_lower)) <= 1 || (change_label == false))
           break;
       end
       
%        Heuristic_t_temp = toc(Heuristic_t_Start);
%        [~,temp_maximum_revenue_index] = max(revenue);
%        temp_Heuristic_best_X = all_policy_x{temp_maximum_revenue_index};
%        temp_Heuristic_best_revenue = calculate_revenue(temp_Heuristic_best_X);
%        fprintf('Running time: %d minutes and %f seconds\n', floor(Heuristic_t_temp/60), rem(Heuristic_t_temp,60));
%        fprintf('Temp Heuristic best revenue = %f \n', temp_Heuristic_best_revenue);
%        if Heuristic_t_temp > 120*60
%            break;
%        end
    end
    [~,maximum_revenue_index] = max(revenue);
    Heuristic_best_X = all_policy_x{maximum_revenue_index};
    Heuristic_best_revenue = calculate_revenue(Heuristic_best_X);
    
    Heuristic_t_End = toc(Heuristic_t_Start);
    fprintf('Running time of the heuristic algorithm: %d minutes and %f seconds\n', floor(Heuristic_t_End/60), rem(Heuristic_t_End,60));
    
end

function[revenue_K,policy_x_K] = calculate_revenue_for_K (K)
    global I J L epsilon num_group_in_one_cluster optimal_V epsilon_prime jump_size cluster_iteration_limitation;
    global utility_v0 revenue_matrix_r utility_matrix_v;
    global lower_case_C upper_case_C order_matrix;
    
    L_hat = J - (num_group_in_one_cluster - 1)*L; %single group
    epsilon_prime = (1 + epsilon)^(1/(2*L_hat))-1; %(1+epsilon_prime)
   
    policy_x_K = [];
    
    %initialization the matrix of optimal V
    optimal_V = zeros(J+1,K+1,upper_case_C+1); %(j^th groups to J, log s_{j-1},1 to j-1 item number)
    Gamma = (1+epsilon_prime)^K;
    for j = J :-1: num_group_in_one_cluster*L+1
        for k = K-J+j : K 
            optimal_V(j,k+1,:) = -inf;
        end
        for k = 0 : j - (num_group_in_one_cluster-1)*L -2
            optimal_V(j,k+1,:) = -inf;
        end
    end

    for j = num_group_in_one_cluster*L - (num_group_in_one_cluster - 1) : -num_group_in_one_cluster : 1 
        for k = K-(L_hat - floor((j-1)/num_group_in_one_cluster) ) + 1 : K 
            optimal_V(j,k+1,:) = -inf;
        end
        for k = 0 : floor((j-1)/num_group_in_one_cluster) - 1
            optimal_V(j,k+1,:) = -inf;
        end
    end
    optimal_V(J+1,1:K,:) = -inf;

    matrix_x = {}; 
    optimal_x_history = {};

    %for single groups
    for j = J :-1: num_group_in_one_cluster*L + 1 
        for k = (j - (num_group_in_one_cluster-1)*L -1) : K-(J-j+1) % log s_{j-1)
            for occupied_c = 0 : upper_case_C
                temp_optimal_rev = -inf; %revenue from j to J
                for subpolicy_index = 1 : size(order_matrix{j},2)
                    for i = 0 : min(lower_case_C, upper_case_C - occupied_c)
                        temp_xj = zeros(1,I);
                        label = true;
                        for ordering_index = 1 : i
                            if order_matrix{j}(ordering_index,subpolicy_index) == -1
                                label = false;
                                break;
                            end
                            temp_xj(order_matrix{j}(ordering_index,subpolicy_index)) = 1;
                        end
                        if label == false
                            continue;
                        end  
                        temp_gamma = 1 + temp_xj*utility_matrix_v(:,j);
                        log_gamma_j_upper = ceil(log(temp_gamma)/log(1+epsilon_prime));
                        if log_gamma_j_upper == 0
                           log_gamma_j_upper = 1;
                        end
                        if ( 1 <= log_gamma_j_upper && log_gamma_j_upper <= K-k) % gamma_j_lower_bound < temp_gamma &&
                            gamma_j_upper_bound = (1+epsilon_prime)^(log_gamma_j_upper);
                            Vj_denominator =  (1+ (utility_v0-1)/Gamma)*gamma_j_upper_bound;
                            temp_rev = temp_xj * (revenue_matrix_r(:,j).* utility_matrix_v(:,j))/Vj_denominator + optimal_V(j+1,log_gamma_j_upper+k+1,(sum(temp_xj) + occupied_c +1));
                            if temp_rev > temp_optimal_rev 
                                temp_optimal_rev = temp_rev;
                                matrix_x{j,k+1,occupied_c+1} = temp_xj;
                                optimal_x_history{j,k+1,occupied_c+1} = [j+1,log_gamma_j_upper+k+1,(sum(temp_xj) + occupied_c +1)];
                            end %if
                        end %if
                    end %end for i
                end
                optimal_V(j,k+1,occupied_c+1) = temp_optimal_rev;
            end
        end
    end

    %for clusters
    for j = num_group_in_one_cluster*L - (num_group_in_one_cluster - 1) : -num_group_in_one_cluster : 1 
        for k = floor((j-1)/num_group_in_one_cluster)  : K-(L_hat - floor((j-1)/num_group_in_one_cluster) ) % log s_{j-1)
            for occupied_c = 0 : upper_case_C
                temp_optimal_rev = -inf; %revenue from j to J
                
                %initialization
                current_x_for_cluster = [ones(floor(I/2),num_group_in_one_cluster);zeros(I - floor(I/2),num_group_in_one_cluster)];
                [gamma_for_cluster_j,~] = calculate_gamma_Lambda_cluster (j, current_x_for_cluster);
                log_total_gamma_bundle_upper = ceil(log(gamma_for_cluster_j)/log(1+epsilon_prime));
                if (log_total_gamma_bundle_upper > K-k) || (sum(sum(current_x_for_cluster)) > (upper_case_C - occupied_c))
                    current_x_for_cluster = [ones(floor(I/4),num_group_in_one_cluster);zeros(I - floor(I/4),num_group_in_one_cluster)];
                    [gamma_for_cluster_j,~] = calculate_gamma_Lambda_cluster (j, current_x_for_cluster);
                    log_total_gamma_bundle_upper = ceil(log(gamma_for_cluster_j)/log(1+epsilon_prime));
                    if (log_total_gamma_bundle_upper > K-k) || (sum(sum(current_x_for_cluster)) > (upper_case_C - occupied_c))
                        current_x_for_cluster = zeros(I,num_group_in_one_cluster);
                    end
                end
%                 number_of_1_in_init = min(lower_case_C, floor(upper_case_C/J));
%                 while(1)                    
%                     current_x_for_cluster = [ones(number_of_1_in_init,num_group_in_one_cluster);zeros(I - number_of_1_in_init,num_group_in_one_cluster)];
%                     [gamma_for_cluster_j,~] = calculate_gamma_Lambda_cluster (j, current_x_for_cluster);
%                     log_total_gamma_bundle_upper = ceil(log(gamma_for_cluster_j)/log(1+epsilon_prime));
%                     if (log_total_gamma_bundle_upper > K-k) || (sum(sum(current_x_for_cluster)) > (upper_case_C - occupied_c))
%                         number_of_1_in_init = number_of_1_in_init - 1;
%                     else
%                         break;
%                     end
%                     if number_of_1_in_init < 0 
%                         break;
%                     end
%                 end
                
                current_index = 1;
                num_of_not_update = 0;
                num_of_iteration = 0;
                temp_local_optimal_rev_for_this_cluster = -inf;
                while(1)
                    %stopping condition
                    if num_of_not_update == num_group_in_one_cluster || num_of_iteration > cluster_iteration_limitation
%                             if num_of_not_update == num_group_in_one_cluster
%                                 fprintf('\n');
%                                 fprintf('The heuristic algorithm is stopped by the no update. The number of iterations is %d.\n', num_of_iteration);
%                             else
%                                 fprintf('\n');
%                                 fprintf('The heuristic algorithm is stopped by the iteration limit.\n');
%                             end
                        break;
                    end

                    [temp_local_optimal_rev_for_this_cluster,current_x_for_cluster, update_flag] = find_better_policy_in_group_j(j,current_index,current_x_for_cluster,temp_local_optimal_rev_for_this_cluster,Gamma,K,k,occupied_c);
                    if update_flag == 0 
                        num_of_not_update = num_of_not_update + 1;
                    else
                        num_of_not_update = 0;
                    end
                    num_of_iteration = num_of_iteration + 1;
                    current_index = mod(current_index + 1, num_group_in_one_cluster); 
                    if current_index == 0
                        current_index = num_group_in_one_cluster;
                    end
                end
                [gamma_for_cluster_j,~] = calculate_gamma_Lambda_cluster (j, current_x_for_cluster);
                log_total_gamma_bundle_upper = ceil(log(gamma_for_cluster_j)/log(1+epsilon_prime));
                if log_total_gamma_bundle_upper == 0
                    log_total_gamma_bundle_upper = 1;
                end
                if temp_local_optimal_rev_for_this_cluster > temp_optimal_rev 
                    temp_optimal_rev = temp_local_optimal_rev_for_this_cluster;
                    matrix_x{j,k+1,occupied_c+1} = current_x_for_cluster';
                    optimal_x_history{j,k+1,occupied_c+1} = [j+num_group_in_one_cluster,log_total_gamma_bundle_upper+k+1,(sum(sum(current_x_for_cluster))+occupied_c+1)];
                end

                optimal_V(j,k+1,occupied_c+1) = temp_optimal_rev;


            end
        end
    end
    revenue_K = max(optimal_V(1,1,:));
    if revenue_K ~= -inf
        optimal_x = zeros(J,I);
        last_j = 1;
        last_k = 1;
        [~,last_c] = max(optimal_V(1,1,:));
        while (last_j <= J)
            if size(matrix_x{last_j,last_k,last_c},1) > 1
                optimal_x(last_j:last_j+num_group_in_one_cluster - 1,:) = matrix_x{last_j,last_k,last_c};
            else
                optimal_x(last_j,:) = matrix_x{last_j,last_k,last_c};
            end
            ttt = optimal_x_history{last_j, last_k, last_c};
            last_j = ttt(1);
            last_k = ttt(2);
            last_c = ttt(3);
        end
        policy_x_K = optimal_x';
    end
%          clc;
%         fprintf('Now is processing K = %d.\n', K);

end


function  [temp_local_optimal_rev_for_this_cluster,current_x_for_cluster, update_flag] = find_better_policy_in_group_j(j_begin,current_index,current_x_for_cluster,temp_local_optimal_rev_for_this_cluster,Gamma,K,k,occupied_c)
    global I utility_v0 num_group_in_one_cluster optimal_V epsilon_prime order_matrix;
    global lower_case_C upper_case_C;
    
    temp_xj = current_x_for_cluster(:,current_index);
    temp_x_of_cluster = current_x_for_cluster;
    update_flag = false;
    
    for subpolicy_index = 1 : size(order_matrix{j_begin-1+current_index},2)
        for i = 0 : min(lower_case_C, (upper_case_C - occupied_c))
            temp_xj = zeros(1,I);
            label = true;
            for ordering_index = 1 : i
                if order_matrix{j_begin-1+current_index}(ordering_index,subpolicy_index) == -1
                    label = false;
                    break;
                end
                temp_xj(order_matrix{j_begin-1+current_index}(ordering_index,subpolicy_index)) = 1;
            end
            if label == false
                continue;
            end  
            temp_x_of_cluster(:,current_index) = temp_xj;
            [temp_gamma,temp_Lambda] = calculate_gamma_Lambda_cluster(j_begin,temp_x_of_cluster);
            log_total_gamma_cluster_upper = ceil(log(temp_gamma)/log(1+epsilon_prime));
            if log_total_gamma_cluster_upper == 0
                log_total_gamma_cluster_upper = 1;
            end
            if ((1 <= log_total_gamma_cluster_upper) && (log_total_gamma_cluster_upper <= (K-k)))
                if (sum(sum(temp_x_of_cluster)) <= (upper_case_C - occupied_c)) 
                    gamma_cluster_upper_bound = (1+epsilon_prime)^(log_total_gamma_cluster_upper);
                    temp_revenue = temp_Lambda/(gamma_cluster_upper_bound * (1+ (utility_v0-1)/Gamma)) + optimal_V(j_begin + num_group_in_one_cluster, log_total_gamma_cluster_upper+k+1,(sum(sum(temp_x_of_cluster))+occupied_c+1) );
                    if temp_revenue > temp_local_optimal_rev_for_this_cluster
                        temp_local_optimal_rev_for_this_cluster = temp_revenue;
                        current_x_for_cluster = temp_x_of_cluster;
                        update_flag = true;
                    end
                end
            end
        end
    end

end

function [revenue] = calculate_revenue (current_x)
    global I J L num_group_in_one_cluster utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi;
    vx_matrix = utility_matrix_v .* current_x;
    vx_matrix(end+1,:) = ones(1,J);
    
    revenue = 0;
    Lambda_of_cluster_vector = zeros(1,( J - (num_group_in_one_cluster - 1)*L));
    gamma_of_cluster_vector = zeros(1,( J - (num_group_in_one_cluster - 1)*L));
    
    for index_L = 1 : L
        j = num_group_in_one_cluster * (index_L - 1) + 1;
        j_end = num_group_in_one_cluster * index_L;
        current_x_for_cluster = current_x(:,j:j_end);
        
        if num_group_in_one_cluster == 2
            gamma_of_cluster_vector(index_L) = 1 + current_x_for_cluster(:,1)'*utility_matrix_v(:,j) + current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1) ...
                                + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1));
            Lambda_of_cluster_vector(index_L) = (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))...
                *(1 + interation_para_phi(j,j+1)*(current_x_for_cluster(:,2)' * utility_matrix_v(:,j+1))) ...
                + (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1))) ...
                *(1 + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)' * utility_matrix_v(:,j)));
        elseif num_group_in_one_cluster == 3
            gamma_of_cluster_vector(index_L) = 1 + current_x_for_cluster(:,1)'*utility_matrix_v(:,j) + current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1) + current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2) ...
                                + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                                + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2));
            Lambda_of_cluster_vector(index_L) = (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j))) + (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))...
                                + (current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2))) + ...
                                + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                                + interation_para_phi(j,j+1)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                                + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                                + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j));
        elseif num_group_in_one_cluster == 4
            gamma_of_cluster_vector(index_L) = 1 + current_x_for_cluster(:,1)'*utility_matrix_v(:,j) + current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1) + current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2) + current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3)...
                                + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                                + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+3)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                                * (current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3));

                Lambda_of_cluster_vector(index_L) = (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j))) + (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))...
                                + (current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2))) + (current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))...
                                + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                                + interation_para_phi(j,j+1)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                                + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                                + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...  
                                + interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))... 
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j)) ...
                                + interation_para_phi(j,j+3)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...  
                                + interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))... 
                                + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)'* (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                                + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))... 
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                                * (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                                * (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                                * (current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                                + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                                * (current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1));


        end
    end
    
    for index_j = num_group_in_one_cluster * L + 1 : J
        gamma_of_cluster_vector(index_j - (num_group_in_one_cluster - 1)*L) = 1 + current_x(:,index_j)'*utility_matrix_v(:,index_j);
        Lambda_of_cluster_vector(index_j - (num_group_in_one_cluster - 1)*L) = current_x(:,index_j)'* (revenue_matrix_r(:,index_j).* utility_matrix_v(:,index_j));
    end

    Gamma = prod(gamma_of_cluster_vector);

    for cluster_index = 1 : (L + J - num_group_in_one_cluster * L)
        revenue = revenue + Lambda_of_cluster_vector(cluster_index)/(gamma_of_cluster_vector(cluster_index)*((utility_v0-1)/Gamma + 1));
    end

end

function [gamma_for_cluster,Lambda_for_cluster] = calculate_gamma_Lambda_cluster (j, current_x_for_cluster)
    global I J num_group_in_one_cluster revenue_matrix_r utility_matrix_v interation_para_phi correlation_index;
    
    if num_group_in_one_cluster == 2
        gamma_for_cluster = 1 + current_x_for_cluster(:,1)'*utility_matrix_v(:,j) + current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1) ...
                            + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1));
        Lambda_for_cluster = (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))...
            *(1 + interation_para_phi(j,j+1)*(current_x_for_cluster(:,2)' * utility_matrix_v(:,j+1))) ...
            + (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1))) ...
            *(1 + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)' * utility_matrix_v(:,j)));
    elseif num_group_in_one_cluster == 3
        gamma_for_cluster = 1 + current_x_for_cluster(:,1)'*utility_matrix_v(:,j) + current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1) + current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2) ...
                            + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                            + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2));
        Lambda_for_cluster = (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j))) + (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))...
                            + (current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2))) + ...
                            + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                            + interation_para_phi(j,j+1)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                            + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                            + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j));
    elseif num_group_in_one_cluster == 4
        gamma_for_cluster = 1 + current_x_for_cluster(:,1)'*utility_matrix_v(:,j) + current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1) + current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2) + current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3)...
                            + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                            + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+3)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                            * (current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3));
        
            Lambda_for_cluster = (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j))) + (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))...
                            + (current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2))) + (current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))...
                            + interation_para_phi(j,j+1)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                            + interation_para_phi(j,j+1)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                            + interation_para_phi(j,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...
                            + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...  
                            + interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))... 
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j+1,j+2)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j)) ...
                            + interation_para_phi(j,j+3)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))...  
                            + interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))... 
                            + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)'* (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))...
                            + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+2,j+3)*(current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))... 
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                            * (current_x_for_cluster(:,1)' * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                            * (current_x_for_cluster(:,2)' * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                            * (current_x_for_cluster(:,3)' * (revenue_matrix_r(:,j+2).* utility_matrix_v(:,j+2)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1))*(current_x_for_cluster(:,4)'*utility_matrix_v(:,j+3))...
                            + interation_para_phi(j,j+1)*interation_para_phi(j,j+2)*interation_para_phi(j,j+3)*interation_para_phi(j+1,j+2)*interation_para_phi(j+1,j+3)*interation_para_phi(j+2,j+3)...
                            * (current_x_for_cluster(:,4)' * (revenue_matrix_r(:,j+3).* utility_matrix_v(:,j+3)))*(current_x_for_cluster(:,1)'*utility_matrix_v(:,j))*(current_x_for_cluster(:,3)'*utility_matrix_v(:,j+2))*(current_x_for_cluster(:,2)'*utility_matrix_v(:,j+1));
              
                            
    end
    
end