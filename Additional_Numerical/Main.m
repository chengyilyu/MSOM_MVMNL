clear;
global I J L num_group_in_one_cluster utility_v0 revenue_matrix_r utility_matrix_v interation_para_phi1 interation_para_phi2 interation_para_phi3 interation_para_phi4 interation_para_phi5;

rand('seed',100*sum(clock));
I = 10; %total number of items in each group
J = 5;  %total number of groups
L = 1;  %total number of bundles
num_group_in_one_cluster = J;
interationTerm1 = 0.2;
interationTerm2 = 0.4;
interationTerm3 = 0.6; 
interationTerm4 = 0.8;
interationTerm5 = 0.9;
utility_v0 = 1.5;
% load('uncapacitated_(3,3,1,2)_(3).mat');

%%%
% revenue_matrix_r = sort(rand(I,J)*0.5,1,'descend'); % each column is non-increasing
% revenue_matrix_r(:,1:floor(I/3)) = rand(I,floor(I/3))*2;
% revenue_matrix_r = sort(revenue_matrix_r,1,'descend');

% %% Random
% revenue_matrix_r = sort(rand(I,J)*2,1,'descend');

% %%% Similar Price
% revenue_matrix_r(1:I,1) = sort(rand(I,1)*2,1,'descend');
% revenue_matrix_r(1:I,2) = max(revenue_matrix_r(1:I,1) - 0.05,0);
% revenue_matrix_r(1:I,3) = max(revenue_matrix_r(1:I,2) - 0.05,0);
% revenue_matrix_r(1:I,4) = max(revenue_matrix_r(1:I,3) - 0.05,0);
% revenue_matrix_r(1:I,5) = max(revenue_matrix_r(1:I,4) - 0.05,0);

basePrice = 0;
baseVec = [1 : -0.1 : 0];
baseVec = baseVec(1:I);
revenue_matrix_r(1:I,1) = baseVec + ones(1,I)*basePrice;
revenue_matrix_r(1:I,2) = baseVec + ones(1,I)*(basePrice+0.1);
revenue_matrix_r(1:I,3) = baseVec + ones(1,I)*(basePrice+0.2);
revenue_matrix_r(1:I,4) = baseVec + ones(1,I)*(basePrice+0.3);
revenue_matrix_r(1:I,5) = baseVec + ones(1,I)*(basePrice+0.4);

utility_matrix_v = 1 + rand(I,J);

interation_para_phi1 = ones(J,J);
for cluster_num = 1 : L
    for index_1 = 1 : num_group_in_one_cluster - 1
        for index_2 = (num_group_in_one_cluster*cluster_num - index_1 + 1) : num_group_in_one_cluster*cluster_num
            interation_para_phi1(num_group_in_one_cluster*cluster_num - index_1, index_2) = interationTerm1;
            interation_para_phi1(index_2, num_group_in_one_cluster*cluster_num - index_1) = interationTerm1;
        end
    end
end

interation_para_phi2 = ones(J,J);
for cluster_num = 1 : L
    for index_1 = 1 : num_group_in_one_cluster - 1
        for index_2 = (num_group_in_one_cluster*cluster_num - index_1 + 1) : num_group_in_one_cluster*cluster_num
            interation_para_phi2(num_group_in_one_cluster*cluster_num - index_1, index_2) = interationTerm2;
            interation_para_phi2(index_2, num_group_in_one_cluster*cluster_num - index_1) = interationTerm2;
        end
    end
end

interation_para_phi3 = ones(J,J);
for cluster_num = 1 : L
    for index_1 = 1 : num_group_in_one_cluster - 1
        for index_2 = (num_group_in_one_cluster*cluster_num - index_1 + 1) : num_group_in_one_cluster*cluster_num
            interation_para_phi3(num_group_in_one_cluster*cluster_num - index_1, index_2) = interationTerm3;
            interation_para_phi3(index_2, num_group_in_one_cluster*cluster_num - index_1) = interationTerm3;
        end
    end
end

interation_para_phi4 = ones(J,J);
for cluster_num = 1 : L
    for index_1 = 1 : num_group_in_one_cluster - 1
        for index_2 = (num_group_in_one_cluster*cluster_num - index_1 + 1) : num_group_in_one_cluster*cluster_num
            interation_para_phi4(num_group_in_one_cluster*cluster_num - index_1, index_2) = interationTerm4;
            interation_para_phi4(index_2, num_group_in_one_cluster*cluster_num - index_1) = interationTerm4;
        end
    end
end

interation_para_phi5 = ones(J,J);
for cluster_num = 1 : L
    for index_1 = 1 : num_group_in_one_cluster - 1
        for index_2 = (num_group_in_one_cluster*cluster_num - index_1 + 1) : num_group_in_one_cluster*cluster_num
            interation_para_phi5(num_group_in_one_cluster*cluster_num - index_1, index_2) = interationTerm5;
            interation_para_phi5(index_2, num_group_in_one_cluster*cluster_num - index_1) = interationTerm5;
        end
    end
end

% save('uncapacitated_(5,5)_(00)', "revenue_matrix_r","utility_matrix_v");
load('uncapacitated_(10,5)_(00)');

fprintf('There are %d groups, %d items in each group and %d clusters, each cluster has %d group. \n',J,I,L,num_group_in_one_cluster);


[bestPolicy_GroupwiseMNL, bestPolicy_MVMNLwithoutInter, bestPolicy_MVMNLwithInter1, bestPolicy_MVMNLwithInter2, ...
    bestPolicy_MVMNLwithInter3, bestPolicy_MVMNLwithInter4, bestPolicy_MVMNLwithInter5] = EnumerationFunction();

bestRevenueGroupwiseMNL_UnderMVMNLwithInter1 = calculate_revenue_MVMNL_withInteraction (bestPolicy_GroupwiseMNL, interation_para_phi1);
bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter1 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithoutInter, interation_para_phi1);
bestRevenueVMNLwithInter_UnderMVMNLwithInter1 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithInter1,interation_para_phi1);


bestRevenueGroupwiseMNL_UnderMVMNLwithInter2 = calculate_revenue_MVMNL_withInteraction (bestPolicy_GroupwiseMNL, interation_para_phi2);
bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter2 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithoutInter, interation_para_phi2);
bestRevenueVMNLwithInter_UnderMVMNLwithInter2 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithInter2,interation_para_phi2);


bestRevenueGroupwiseMNL_UnderMVMNLwithInter3 = calculate_revenue_MVMNL_withInteraction (bestPolicy_GroupwiseMNL, interation_para_phi3);
bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter3 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithoutInter, interation_para_phi3);
bestRevenueVMNLwithInter_UnderMVMNLwithInter3 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithInter3,interation_para_phi3);

bestRevenueGroupwiseMNL_UnderMVMNLwithInter4 = calculate_revenue_MVMNL_withInteraction (bestPolicy_GroupwiseMNL, interation_para_phi4);
bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter4 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithoutInter, interation_para_phi4);
bestRevenueVMNLwithInter_UnderMVMNLwithInter4 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithInter4,interation_para_phi4);

bestRevenueGroupwiseMNL_UnderMVMNLwithInter5 = calculate_revenue_MVMNL_withInteraction (bestPolicy_GroupwiseMNL, interation_para_phi5);
bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter5 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithoutInter, interation_para_phi5);
bestRevenueVMNLwithInter_UnderMVMNLwithInter5 = calculate_revenue_MVMNL_withInteraction (bestPolicy_MVMNLwithInter5,interation_para_phi5);

fprintf('---------------------------------------------\n')
fprintf('The interation term is: %.1f \n', interationTerm1);
fprintf('---\n')
fprintf('The best policy under GroupwiseMNL is: \n');
bestPolicy_GroupwiseMNL
sum(bestPolicy_GroupwiseMNL,1)
fprintf('The corresponding revenue under GroupwiseMNL is: %.5f\n', bestRevenueGroupwiseMNL_UnderMVMNLwithInter1);
fprintf('---\n')
fprintf('The best policy under MVMNL without Interaction is: \n');
bestPolicy_MVMNLwithoutInter
sum(bestPolicy_MVMNLwithoutInter,1)
fprintf('The corresponding revenue under MVMNL with interaction is: %.5f\n', bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter1);
fprintf('---\n')
fprintf('The best policy under MVMNL with Interaction is: \n');
bestPolicy_MVMNLwithInter1
sum(bestPolicy_MVMNLwithInter1,1)
fprintf('The corresponding revenue under MVMNL with interaction is: %.5f\n', bestRevenueVMNLwithInter_UnderMVMNLwithInter1);
fprintf('\n')
fprintf('---------------------------------------------\n')
fprintf('The interation term is: %.1f \n', interationTerm2);
fprintf('---\n')
fprintf('Groupwise MNL revenue under GroupwiseMNL is: %.5f\n', bestRevenueGroupwiseMNL_UnderMVMNLwithInter2);
fprintf('MVMNL without Interaction revenue under MVMNL with interaction is: %.5f\n', bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter2);
fprintf('---\n')
fprintf('The best policy under MVMNL with Interaction is: \n');
bestPolicy_MVMNLwithInter2
sum(bestPolicy_MVMNLwithInter2,1)
fprintf('The corresponding revenue under MVMNL with interaction is: %.5f\n', bestRevenueVMNLwithInter_UnderMVMNLwithInter2);
fprintf('\n')
fprintf('---------------------------------------------\n')
fprintf('The interation term is: %.1f \n', interationTerm3);
fprintf('---\n')
fprintf('Groupwise MNL revenue under GroupwiseMNL is: %.5f\n', bestRevenueGroupwiseMNL_UnderMVMNLwithInter3);
fprintf('MVMNL without Interaction revenue under MVMNL with interaction is: %.5f\n', bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter3);
fprintf('---\n')
fprintf('The best policy under MVMNL with Interaction is: \n');
bestPolicy_MVMNLwithInter3
sum(bestPolicy_MVMNLwithInter3,1)
fprintf('The corresponding revenue under MVMNL with interaction is: %.5f\n', bestRevenueVMNLwithInter_UnderMVMNLwithInter3);

fprintf('\n')
fprintf('---------------------------------------------\n')
fprintf('The interation term is: %.1f \n', interationTerm4);
fprintf('---\n')
fprintf('Groupwise MNL revenue under GroupwiseMNL is: %.5f\n', bestRevenueGroupwiseMNL_UnderMVMNLwithInter4);
fprintf('MVMNL without Interaction revenue under MVMNL with interaction is: %.5f\n', bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter4);
fprintf('---\n')
fprintf('The best policy under MVMNL with Interaction is: \n');
bestPolicy_MVMNLwithInter4
sum(bestPolicy_MVMNLwithInter4,1)
fprintf('The corresponding revenue under MVMNL with interaction is: %.5f\n', bestRevenueVMNLwithInter_UnderMVMNLwithInter4);


fprintf('\n')
fprintf('---------------------------------------------\n')
fprintf('The interation term is: %.1f \n', interationTerm5);
fprintf('---\n')
fprintf('Groupwise MNL revenue under GroupwiseMNL is: %.5f\n', bestRevenueGroupwiseMNL_UnderMVMNLwithInter5);
fprintf('MVMNL without Interaction revenue under MVMNL with interaction is: %.5f\n', bestRevenueMVMNLwithoutInter_UnderMVMNLwithInter5);
fprintf('---\n')
fprintf('The best policy under MVMNL with Interaction is: \n');
bestPolicy_MVMNLwithInter5
sum(bestPolicy_MVMNLwithInter5,1)
fprintf('The corresponding revenue under MVMNL with interaction is: %.5f\n', bestRevenueVMNLwithInter_UnderMVMNLwithInter5);
