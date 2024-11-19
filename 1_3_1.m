close all; clear; clc;
NACA = '0006';
alpha1 = 10; % AOA [deg]

%% Experimental Data
Exp0012Data = load("Dataset 9E6 Reynolds.csv"); % NACA 0006 (section aoa [deg], section C_L)
Exp2412Data = load("NACA_2412.csv"); % NACA 0006 (section aoa [deg], section C_L)
Exp4412Data = load("NACA_4412.txt"); % NACA 0012 (section aoa [deg], section C_L)



 %% 1.1.1 Plotting for comparing airfoils
 [a1, c_l1, x1,y_camber1,x_u1,y_u1,x_L1,y_L1] = plotAirfoil('0012',100);
 [a2, c_l2, x2,y_camber2,x_u2,y_u2,x_L2,y_L2] = plotAirfoil('2421',100);
 [a3, c_l3, x3,y_camber3,x_u3,y_u3,x_L3,y_L3] = plotAirfoil('4412',100);
figure('Position', [40 60 700 400]); hold on; grid on;
title("Comparing NACA Airfoils")
xlabel("Distance along Chord");
ylabel("Thickness");
plot(x_u1,y_u1,'r.','LineWidth',2);
plot(x_L1,y_L1,'r.','LineWidth',2);
plot(x_u2,y_u2,'b.','LineWidth',2);
plot(x_L2,y_L2,'b.','LineWidth',2);
plot(x_u3,y_u3,'m.','LineWidth',2);
plot(x_L3,y_L3,'m.','LineWidth',2);
plot(x2, y_camber2,'b');
ylim([-0.2 0.3])
xlim([-0.001 1.001])
legend("NACA 0012","","NACA 2421","","NACA 4412","","NACA 2421 Camber")

%% Finding true steady state value for c_l convergence
[~,~,~,~,x_u,y_u,x_L,y_L] = plotAirfoil(NACA,2000);
x_B = [flip(x_L(2:end)),x_u(1:end)];
y_B = [flip(y_L(2:end)),y_u(1:end)];
C_l_true = Vortex_Panel(x_B,y_B,alpha1);

% figure(); hold on;
% x_B_new = x_B;
% y_B_new = y_B;
% c = linspace(1,10,length(x_B_new));
% scatter(x_B_new, y_B_new, [], c);
% ylim([-0.2 0.3])
% xlim([-0.001 1.001])

%% 1.2.1 C_l convergence
j = 1;
for i=4:1:100
    [~,~,~,~,x_u,y_u,x_L,y_L] = plotAirfoil(NACA,i);
    x_B = [flip(x_L(2:end)),x_u(1:end)];
    y_B = [flip(y_L(2:end)),y_u(1:end)];
    c_l(i) = Vortex_Panel(x_B,y_B,alpha1);
    x_loc(i) = i;
    C_l_goal(i) = 100*abs(c_l(i)-(C_l_true))/C_l_true;
end
k = find(C_l_goal(5:end)<1,1)+5;

figure();
hold on; grid on; grid minor;
plot(x_loc, c_l,'LineWidth',2);
title("NACA 0006: Convergence of c_l");
xlabel("Number of Total Panels on Airfoil");
ylabel("Sectional Coefficient of Lift, c_l");


%% 1.2.1 Angle of attack
i=1;
for a=-8:1:14
    [a2,cl2,~,~,x_u2,y_u2,x_L2,y_L2] = plotAirfoil('0012',100);
    x_B2 = [flip(x_L2(2:end)),x_u2(1:end)];
    y_B2 = [flip(y_L2(2:end)),y_u2(1:end)];
    c_l2a(i) = Vortex_Panel(x_B2,y_B2,a);
    [a12,cl12,~,~,x_u12,y_u12,x_L12,y_L12] = plotAirfoil('2412',100);
    x_B12 = [flip(x_L12(2:end)),x_u12(1:end)];
    y_B12 = [flip(y_L12(2:end)),y_u12(1:end)];
    c_l12a(i) = Vortex_Panel(x_B12,y_B12,a);
    [a412,cl412,~,~,x_u412,y_u412,x_L412,y_L412] = plotAirfoil('4412',100);
    x_B412 = [flip(x_L412(2:end)),x_u412(1:end)];
    y_B412 = [flip(y_L412(2:end)),y_u412(1:end)];
    c_l412a(i) = Vortex_Panel(x_B412,y_B412,a);





    i = i + 1;
end

figure();
hold on; grid on; grid minor;
plot(-8:1:14, c_l2a,'b','LineWidth',1);
plot(-8:1:14, c_l12a,'r','LineWidth',1);
plot(-8:1:14, c_l412a,'g','LineWidth',1);

plot(Exp0012Data(:,1), Exp0012Data(:,2),'b--','LineWidth',1);
plot(Exp2412Data(:,1), Exp2412Data(:,2),'r--','LineWidth',1);
plot(Exp4412Data(:,1), Exp4412Data(:,2),'g--','LineWidth',1);

plot(a2, cl2,'b.-','LineWidth',1);
plot(a12, cl12,'r.-','LineWidth',1);
plot(a412, cl412,'g.-','LineWidth',1);
title("NACA Airfoil c_l");
xlabel("Angle of Attack [deg]");
ylabel("Sectional Coefficient of Lift, c_l");
legend("NACA 0012 Vortex Panel","NACA 2412 Vortex Panel","NACA 4412 Vortex Panel", "NACA 0012 Experimental", "NACA 2412 Experimental", "NACA 4412 Experimental", "NACA 0012 TAT", "NACA 2412 TAT", "NACA 4412 TAT","Location","southeast")
xlim([-8, 10]);
