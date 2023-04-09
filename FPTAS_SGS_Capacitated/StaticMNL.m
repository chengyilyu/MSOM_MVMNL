function [order_matrix] = StaticMNL(I, J, revenue_matrix_r, utility_matrix_v)

order_matrix = {};

% lower_case_C %local capacity
% upper_case_C %global capacity
% I = 3; %total number of items in each group
% J = 3;  %total number of groups
% revenue_matrix_r = sort(rand(I,J),1,'descend'); % each column is non-increasing
% utility_matrix_v = rand(I,J);

for j = 1 : J
    temp_v = utility_matrix_v(:,j);
    temp_r = revenue_matrix_r(:,j);
    temp_v(I+1) = 0;
    temp_r(I+1) = 0;
    order_matrix_of_j = [];
    
    interaction_point = [];
    interaction_line_num = [];
    for i1 = 1 : I
        for i2 = i1+1 : I+1
            cross_value = (temp_v(i1) * temp_r(i1) - temp_v(i2) * temp_r(i2)) / (temp_v(i1) - temp_v(i2));
            interaction_point(end+1) = cross_value;
            interaction_line_num(1,end+1) = i1;
            interaction_line_num(2,end) = i2;
        end
    end
    
    [interaction_point, index] = sort(interaction_point);
    interaction_line_num = interaction_line_num(:,index);
    
    
    [~,ordering_of_I] = sort(temp_v, 'descend');
    [~,index_of_v] = sort(ordering_of_I);
    
    order_matrix_of_j = [];
    ordering_of_I(end) = -1;
    order_matrix_of_j(:,1) = ordering_of_I;
    for index = 1 : size(interaction_point,2)-1
        i1 = interaction_line_num(1,index);
        i2 = interaction_line_num(2,index);
        ordering_of_I(index_of_v(i1)) = i2;
        ordering_of_I(index_of_v(i2)) = i1;
        tmp = index_of_v(i1);
        index_of_v(i1) = index_of_v(i2);
        index_of_v(i2) = tmp;
        
        for below_0 = index_of_v(I+1) : I+1
            ordering_of_I(below_0) = -1;
        end
        
        order_matrix_of_j(:,end+1) = ordering_of_I;
    end
    
    order_matrix{j} = order_matrix_of_j;
end
end