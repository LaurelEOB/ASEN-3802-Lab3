%% ASEN 3802 Lab 3 Part 3 Convergence, Drag, and Twist Study
%   Task 1: Objective was to preform a convergence study to determine the
%   appropriate number of number of terms to use in the code to gain an
%   acceptable level of accuracy. Found how many terms are required to
%   obtain a solution for each with one percent relative error
%   Task 2: Objective was to analyze the Total Drag Coefficient and compare
%   the two components of drag (profile and induced). Profile drag was
%   estimated from prior experimental results
%   Task 3: Objective was to see what impact twist has on the total lift to
%   drag ratio and spac efficiency factor e
%       Authors: Laurel O'Brien, Gabriel Burdan, Brody Ambroggio, William Wallingford
%       Collaborators: Alyxis Ellington, Samantha Sheppard, Professor Hoke
%       Date Last Revised: 12-12-2024
close all; clc; clear;

%% Part 3, Task 1
% Constants for the C180 wing
aoa1 = 4;
b = 35*12 + 10; % span [in]
c_r = 5*12 + 4; %chord at the root [in]
c_t = 3*12 + 7; % chord at the tips [in]
geo_r1 = (2+aoa1); % geometric angle of attack (geometric twist + alpha) at the roots [deg]
geo_t1 = 0+aoa1; % geometric angle of attack (geometric twist + alpha) at the tips [deg]
NACA_r = '2412';
NACA_t = '0012';
aero_r = -2.0781; % zero-lift angle of attack at the root [deg]
aero_t = 0; % zero-lift angle of attack at the tips [deg]
a0_r = 0.1187*180/pi; % cross-sectional lift slope at the root [per radian]
a0_t = 0.1192*180/pi; % cross-sectional lift slope at the tips [per radian]

numConv = 20;

% Values that the coefficients converge to
[~,C_l_true,C_Di_true] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t1,geo_r1,1000);

% Convergence study
for N=1:numConv
    [e(N),c_L(N),c_Di(N)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t1,geo_r1,N);

    % Percent error
    C_l_goal(N) = 100*abs(c_L(N)-(C_l_true))/C_l_true;
    C_Di_goal(N) = 100*abs(c_Di(N)-(C_Di_true))/C_Di_true;
end 

% Finding when percent error is below 1%
k_CL= find(C_l_goal(1:end)<1,1); % Should be 4
k_CDi= find(C_Di_goal(1:end)<1,1); % Should be 8

% Plotting
figure('Position', [40 60 600 300]); hold on; grid on; grid minor;
title("Convergence of C_L");
plot(1:numConv, c_L,"LineWidth",2);
text(k_CL+2, c_L(2), "C_L converges to 1% after "+k_CL+" odd terms");
xlabel("Number of Odd Terms");
ylabel("Lift Coefficient, C_L");
yline(c_L(k_CL),":",'LineWidth',1);
xa1 = xline(k_CL,'-',[num2str(k_CL) ' Terms']);
xa1.LabelVerticalAlignment = "middle";
xa1.LabelHorizontalAlignment = 'center';
legend("","1% of final value");

figure('Position', [40 460 600 300]); hold on; grid on; grid minor;
title("Convergence of C_D_i");
plot(1:numConv, c_Di,"LineWidth",2);
text(k_CDi+2, c_Di(2), "C_D_i converges to 1% after "+k_CDi+" odd terms");
xlabel("Number of Odd Terms");
ylabel("Induced Drag Coefficient, C_D_i");
yline(c_Di(k_CDi),":",'LineWidth',1);
xa2 = xline(k_CDi,'-',[num2str(k_CDi) ' Terms']);
xa2.LabelVerticalAlignment = "middle";
xa2.LabelHorizontalAlignment = 'center';
legend("","1% of final value");




%% Part 3, Task 2

% Data from Abbot's Theory Of Wing Sections
Exp2412Data1 = load("NACA 2412 c_d.csv"); % NACA 2412 (section lift coefficient, section drag coefficient)
Exp2412Data2 = load("NACA 2412 with neg AoA.txt"); % NACA 2412 (section aoa [deg], section C_L)

% Getting C_Di data into a function of aoa
p_a = polyfit(Exp2412Data2(:,1), Exp2412Data2(:,2),4);
aoa_a = Exp2412Data2(1,1):0.1:Exp2412Data2(end,1)+0.1;
c_L_a = polyval(p_a,aoa_a);
p_b = polyfit(Exp2412Data1(:,1), Exp2412Data1(:,2),7);
c_D_b = polyval(p_b,c_L_a); % Profile Drag Coefficient

%Temporary Plotting
figure(); hold on;
xlabel("C_L")
yyaxis right
plot(Exp2412Data1(:,1),Exp2412Data1(:,2));
plot(c_L_a, c_D_b,'r');
ylabel("C_D");
yyaxis left
plot(Exp2412Data2(:,1), Exp2412Data2(:,2));
plot(aoa_a, c_L_a,'r');
ylabel("aoa");

% Goes through and finds total drag coefficients
for i = 1:length(aoa_a)
    geo_r1 = (2+aoa_a(i)); % geometric angle of attack (geometric twist + alpha) at the roots [deg]
    geo_t1 = 0+aoa_a(i); % geometric angle of attack (geometric twist + alpha) at the tips [deg]
    [~,c_L2(i),c_Di2(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t1,geo_r1,20);
    c_D(i) = c_Di2(i) + c_D_b(i);

end

% Important values
[maxCLCD, maxCLCD_I] = max(c_L2./c_D);
maxCLCD_aoa = aoa_a(maxCLCD_I);
C_L_cruise = 0.27085;
actual_aoa = interp1(c_L2,aoa_a,C_L_cruise);

% Interpolating important values
s.c_L = {interp1(aoa_a,c_L2,4),interp1(aoa_a,c_L2,maxCLCD_aoa),interp1(aoa_a,c_L2,actual_aoa)};
s.c_D = {interp1(aoa_a,c_D,4),interp1(aoa_a,c_D,maxCLCD_aoa),interp1(aoa_a,c_D,actual_aoa)};
s.LD = {interp1(aoa_a,c_L2./c_D,4),interp1(aoa_a,c_L2./c_D,maxCLCD_aoa),interp1(aoa_a,c_L2./c_D,actual_aoa)};

% Plot 1
figure(); hold on; grid on; grid minor;
plot(aoa_a,c_D,"LineWidth",2.5)
plot(aoa_a, c_Di2,"LineWidth",2.5)
plot(aoa_a,c_D_b,"LineWidth",2.5)
xlabel("Angle of Attack (\alpha)");
ylabel("Drag Cofficient");
xl = xline(4,'-',['4' char(176) ' AoA']);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
x2 = xline(maxCLCD_aoa,'-',[num2str(maxCLCD_aoa,4) char(176) ' AoA (Max L/D)']);
x2.LabelVerticalAlignment = 'middle';
x2.LabelHorizontalAlignment = 'center';
x3 = xline(actual_aoa,'-',[num2str(actual_aoa) char(176) ' AoA (C180 Cruise)']);
x3.LabelVerticalAlignment = 'middle';
x3.LabelHorizontalAlignment = 'center';
title("Components of Drag at Different AoA");
xlim([-5,20])
ylim([-0.01 0.2])
yline(0,":",'LineWidth',0.5);
xline(0,":",'LineWidth',0.5);
legend("Total Drag Coefficient","Induced Drag Coefficient","Profile Drag Coefficient","","",'Location','north')

% Plot 2
figure(); hold on; grid on; grid minor;
plot(aoa_a,c_L2./c_D,"LineWidth",2.5)
xlabel("Angle of Attack (\alpha)");
ylabel("C_L / C_D");
xl = xline(4,'-',['4' char(176) ' AoA']);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'center';
x2 = xline(maxCLCD_aoa,'-',[num2str(maxCLCD_aoa,4) char(176) ' AoA (Max L/D)']);
x2.LabelVerticalAlignment = 'bottom';
x2.LabelHorizontalAlignment = 'center';
x3 = xline(actual_aoa,'-',[num2str(actual_aoa) char(176) ' AoA (C180 Cruise)']);
x3.LabelVerticalAlignment = "bottom";
x3.LabelHorizontalAlignment = 'center';
title("Ratio of Lift to Drag as a Function of AoA");
xlim([-5,20])
yline(0,":",'LineWidth',0.5);
xline(0,":",'LineWidth',0.5);



%% Part 3, Task 3
% Varying twist angle
twist_angle = -5:0.1:20;
for i = 1:length(twist_angle)
    geo_r1 = 2+twist_angle(i); % geometric angle of attack (geometric twist + alpha) at the roots [deg]
    geo_t1 = 2-twist_angle(i); % geometric angle of attack (geometric twist + alpha) at the tips [deg]
    [e3(i),c_L3(i),c_Di3(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t1,geo_r1,20);
    c_D3(i) = c_Di3(i) + c_D_b(19);

end

% Plotting
figure(); hold on; grid on; grid minor;
plot(twist_angle,c_L3./c_D3,"LineWidth",2.5)
xlabel("Twist Angle (Washout)");
ylabel("C_L / C_D");
title("Effect of Modifying Twist on C_L/C_D, \alpha=2"+char(176));
figure(); hold on; grid on; grid minor;
plot(twist_angle,e3,"LineWidth",2.5)
xlabel("Twist Angle (Washout)");
ylabel("Span Efficiency Factor (e)");
title("Effect of Modifying Twist on e, \alpha=2"+char(176));
