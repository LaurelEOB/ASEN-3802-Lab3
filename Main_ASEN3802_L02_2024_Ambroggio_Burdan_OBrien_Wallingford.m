%% ASEN 3802 Lab 3 Part 2 Prandtl Lifting Line Code
%   Objective was to take the Prandtl lifting line code and use it to
%   recreate figure 5.20 in Anderson (figure of different AR, taper, and delta). 
%       Authors: Laurel O'Brien, Gabriel Burdan, Brody Ambroggio, William Wallingford
%       Collaborators: Alyxis Ellington, Samantha Sheppard, Professor Hoke
%       Date Last Revised: 12-5-2024
close all; clc; clear;

taper_ratio = 0:0.025:1;
c_r = 3;
aoa = 5;

%% Creating figure 5.20 from Anderson using different AR
for AR=4:2:10 % Varying aspect ratio
    for i=1:length(taper_ratio) % Varying taper ratio

        c_t = taper_ratio(i) * c_r;
        b = AR*(c_r+c_t)*0.5;

        [e,c_L,c_Di] = PLLT(b,2*pi,2*pi,c_t,c_r,   0,    0,      aoa,   aoa,50);

        inducedDragFactor(i, (AR/2)-1) = (c_Di*pi*AR/((c_L)^2)) - 1;
    end
end

%% Graphing 
hold on; grid on; grid minor;
plot(taper_ratio, inducedDragFactor, "LineWidth",2)
title("")
xlabel("Taper ratio, c_t/c_r");
ylabel("Induced drag factor, \delta");
legend("AR = 4", "AR = 6", "AR = 8", "AR = 10");
title("Recreation of Figure 5.20 in Anderson")

