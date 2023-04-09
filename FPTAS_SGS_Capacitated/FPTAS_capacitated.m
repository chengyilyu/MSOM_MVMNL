function [FPTAS_best_revenue, maximum_revenue_index, FPTAS_best_X] = FPTAS_capacitated ()
    global I J L epsilon;
    global utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi;
    global lower_case_C upper_case_C order_matrix;


    L_hat = J - L;
    epsilon_prime = (1 + epsilon)^(1/(2*L_hat))-1; %(1+epsilon_prime)
    utility_matrix_v_bunlde = utility_matrix_v(: , 1 : 2*L); %v_{i,j} for bundles
    utility_matrix_v_single = utility_matrix_v(: , 2*L+1 : J); %v_{i,j} for single groups

    %calculate the maximum M
    FPTAS_t_Start = tic;
    gamma_single_1 = 1+ sum(utility_matrix_v_single,1);
    gamma_bundle_11 = zeros(1,L);
    sum_utility_matrix_v_bunlde = sum(utility_matrix_v_bunlde,1);
    
    if L == 0
        M =  max(ceil(log(gamma_single_1)/log(1+epsilon_prime)));
        Max_K = sum(ceil(log(gamma_single_1)/log(1+epsilon_prime)));
    else
        for bundle_number = 1 : L
            gamma_bundle_11(bundle_number) = 1 + sum_utility_matrix_v_bunlde(2*bundle_number-1) + sum_utility_matrix_v_bunlde(2*bundle_number)+...
                interation_para_phi(2*bundle_number-1,2*bundle_number)*sum_utility_matrix_v_bunlde(2*bundle_number-1)*sum_utility_matrix_v_bunlde(2*bundle_number);
        end
        if J == 2*L
            M = max( ceil(log(gamma_bundle_11)/log(1+epsilon_prime)) );
            Max_K = sum(ceil(log(gamma_bundle_11)/log(1+epsilon_prime)));
        else
            M = max( max(ceil(log(gamma_single_1)/log(1+epsilon_prime)) ), max( ceil(log(gamma_bundle_11)/log(1+epsilon_prime)) ));
            Max_K = sum(ceil(log(gamma_single_1)/log(1+epsilon_prime)) ) + sum(ceil(log(gamma_bundle_11)/log(1+epsilon_prime)));
        end
    end

    % main 
    revenue = zeros(1,L_hat*M )*(-inf);
    all_policy_x = {}; % all x for all K

    for K = L_hat : L_hat*M  %log \Gamma from J to JM
        %initialization the matrix of optimal V
        optimal_V = zeros(J+1,K+1,upper_case_C+1); %(j^th groups to J, log s_{j-1},1 to j-1 item number)
        Gamma = (1+epsilon_prime)^K;
        for j = J :-1: 2*L+1
            for k = K-J+j : K 
                optimal_V(j,k+1,:) = -inf;
            end
            for k = 0 : j - L -2
                optimal_V(j,k+1,:) = -inf;
            end
        end
        for j = 2*L-1 :-2: 1 
            for k = K-(L_hat - floor((j-1)/2) ) + 1 : K 
                optimal_V(j,k+1,:) = -inf;
            end
            for k = 0 : floor((j-1)/2) - 1
                optimal_V(j,k+1,:) = -inf;
            end
        end
        
        optimal_V(J+1,1:K,:) = -inf;

        matrix_x = {}; 
        optimal_x_history = {};

        %for single groups
        for j = J :-1: 2*L+1 
            for k = (j - L - 1) : K-(J-j+1) % log s_{j-1)
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

        %for bundles
        for j = 2*L-1 :-2: 1 
            for k = floor((j-1)/2) : K-(L_hat - floor((j-1)/2) ) % log s_{j-1)
                for occupied_c = 0 : upper_case_C
                    temp_optimal_rev = -inf; %revenue from j to J
                    for subpolicy_j1_index = 1 : size(order_matrix{j},2)
                        for subpolicy_j2_index = 1 : size(order_matrix{j+1},2)
                            for i1 = 0 : min(lower_case_C, upper_case_C - occupied_c)
                                temp_x_j1 = zeros(1,I);
                                label_j1 = true;
                                for ordering_index_1 = 1 : i1
                                    if order_matrix{j}(ordering_index_1,subpolicy_j1_index) == -1
                                        label_j1 = false;
                                        break;
                                    end
                                    temp_x_j1(order_matrix{j}(ordering_index_1,subpolicy_j1_index)) = 1;
                                end
                                if label_j1 == false
                                    break;
                                end
                                
                                for i2 = 0 : min(lower_case_C, upper_case_C - occupied_c - i1)
                                    temp_x_j2 = zeros(1,I);
                                    label_j2 = true;
                                    for ordering_index_2 = 1 : i2
                                        if order_matrix{j+1}(ordering_index_2,subpolicy_j2_index) == -1
                                            label_j2 = false;
                                            break;
                                        end
                                        temp_x_j2(order_matrix{j+1}(ordering_index_2,subpolicy_j2_index)) = 1;
                                    end
                                    if label_j2 == false
                                        break;
                                    end
                                    
                                    temp_gamma = 1 + temp_x_j1*utility_matrix_v(:,j) + temp_x_j2*utility_matrix_v(:,j+1) ...
                                            + interation_para_phi(j,j+1)*(temp_x_j1*utility_matrix_v(:,j))*(temp_x_j2*utility_matrix_v(:,j+1));
                                    log_total_gamma_bundle_upper = ceil(log(temp_gamma)/log(1+epsilon_prime)); 
                                    if log_total_gamma_bundle_upper == 0
                                        log_total_gamma_bundle_upper = 1;
                                    end
                                    if (1 <= log_total_gamma_bundle_upper && log_total_gamma_bundle_upper <= K-k) %gamma_bundle_lower_bound < temp_gamma && 
                                        gamma_bundle_upper_bound = (1+epsilon_prime)^(log_total_gamma_bundle_upper);
                                        Vj_denominator =  (1+ (utility_v0-1)/Gamma)*gamma_bundle_upper_bound;
                                        temp_rev_numerator_part1 = (temp_x_j1 * (revenue_matrix_r(:,j).* utility_matrix_v(:,j)))...
                                            *(1 + interation_para_phi(j,j+1)*(temp_x_j2 * utility_matrix_v(:,j+1)));
                                        temp_rev_numerator_part2 = (temp_x_j2 * (revenue_matrix_r(:,j+1).* utility_matrix_v(:,j+1))) ...
                                            *(1 + interation_para_phi(j,j+1)*(temp_x_j1 * utility_matrix_v(:,j)));
                                        temp_rev_numerator = temp_rev_numerator_part1 + temp_rev_numerator_part2;
                                        temp_rev = temp_rev_numerator / Vj_denominator + optimal_V(j + 2,log_total_gamma_bundle_upper + k+1, (sum(temp_x_j1)+sum(temp_x_j2) + occupied_c +1));
                                        if temp_rev > temp_optimal_rev 
                                            temp_optimal_rev = temp_rev;
                                            matrix_x{j,k+1,occupied_c+1} = [temp_x_j1; temp_x_j2];
                                            optimal_x_history{j,k+1,occupied_c+1} = [j + 2,log_total_gamma_bundle_upper + k+1, (sum(temp_x_j1)+sum(temp_x_j2) + occupied_c +1)];
                                        end %if
                                    end
                                end
                            end
                        end
                    end
                    optimal_V(j,k+1,occupied_c+1) = temp_optimal_rev;
                end
            end
        end
        revenue(K) = max(optimal_V(1,1,:));
        if revenue(K) ~= -inf
            optimal_x = zeros(J,I);
            last_j = 1;
            last_k = 1;
            [~,last_c] = max(optimal_V(1,1,:));
            while (last_j <= J)
                if size(matrix_x{last_j,last_k,last_c},1) > 1
                    optimal_x(last_j:last_j+1,:) = matrix_x{last_j,last_k,last_c};
                else
                    optimal_x(last_j,:) = matrix_x{last_j,last_k,last_c};
                end
                ttt = optimal_x_history{last_j, last_k, last_c};
                last_j = ttt(1);
                last_k = ttt(2);
                last_c = ttt(3);
            end
            all_policy_x{K} = optimal_x';
        end
%         clc;
%          fprintf('Now is processing K = %d of %d.\n', K, M*J);
    end

    [~, maximum_revenue_index] = max(revenue);
    FPTAS_best_X = all_policy_x{maximum_revenue_index};
    FPTAS_best_revenue = calculate_revenue(FPTAS_best_X);

    FPTAS_t_End = toc(FPTAS_t_Start);
    fprintf('\n');
    fprintf('Running time of the FPTAS: %d minutes and %f seconds\n', floor(FPTAS_t_End/60), rem(FPTAS_t_End,60));
    
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