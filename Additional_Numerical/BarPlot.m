%% v_0 = 1.5
x = 0.1 : 0.1 : 0.9;
y = [[3,4,4,5,6]
[3,4,5,5,6]
[4,4,5,5,6]
[4,4,5,5,6]
[4,4,5,5,6]
[4,4,5,5,5]
[4,4,5,5,5]
[4,4,4,4,5]
[4,4,4,4,4]

];


b = bar(x,y);
legend("Group 1", "Group 2", "Group 3", "Group 4", "Group 5");
ylim([0 11])
ylabel("Size of Local Assortment");
xlabel("$$\phi$$",'Interpreter','latex')

b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660 0.6740 0.1880];
b(5).FaceColor = [0.3010 0.7450 0.9330];
title("$$v_0 = 1.5$$",'Interpreter','latex')

%% v_0 = 15
x = 0.1 : 0.1 : 0.9;
y = [[5,6,6,7,8];
[5,5,6,7,7];
[5,5,6,6,7];
[4,5,5,6,6];
[4,5,5,6,6];
[4,5,5,5,6];
[4,4,5,5,5];
[4,4,4,5,5];
[4,4,4,4,4]
];


b = bar(x,y);
legend("Group 1", "Group 2", "Group 3", "Group 4", "Group 5");
ylim([0 11])
ylabel("Size of Local Assortment");
xlabel("$$\phi$$",'Interpreter','latex')

b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660 0.6740 0.1880];
b(5).FaceColor = [0.3010 0.7450 0.9330];
title("$$v_0 = 15$$",'Interpreter','latex')


%% v_0 = 150
x = 0.1 : 0.1 : 0.9;
y = [9,10,10,10,10;
    8,9,9,10,10;
    7,8,8,9,10;
    7,7,7,8,9;
    6,6,7,7,8;
    6,6,6,7,7;
    5,5,6,6,6;
    5,5,5,5,6;
    4,4,4,5,5
];


b = bar(x,y);
legend("Group 1", "Group 2", "Group 3", "Group 4", "Group 5");
ylim([0 11])
ylabel("Size of Local Assortment");
xlabel("$$\phi$$",'Interpreter','latex')

b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0.9290 0.6940 0.1250];
b(4).FaceColor = [0.4660 0.6740 0.1880];
b(5).FaceColor = [0.3010 0.7450 0.9330];
title("$$v_0 = 150$$",'Interpreter','latex')



