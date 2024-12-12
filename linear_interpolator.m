%% linear interpolator
clc
clear 
close all

x1 = -10; % first value
x2 = 0; % second value
y1 = 1.680; % first x val
y2 = 1.729; % second x val
x = -0.67723; % desired location of value
y = y1 + (x-x1)*(y2-y1)/(x2-x1); % output value
disp(y);