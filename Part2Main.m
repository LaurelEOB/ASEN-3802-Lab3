close all; clc; clear;

taper_ratio = 0:0.025:1;
c_r = 3;
aoa = 5;




for AR=4:2:10
    for i=1:length(taper_ratio)
        c_t = taper_ratio(i) * c_r;
        b = AR*(c_r+c_t)*0.5;

        %              PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
        [e,c_L,c_Di] = PLLT(b,2*pi,2*pi,c_t,c_r,   0,    0,      aoa,   aoa,50);

        inducedDragFactor(i, (AR/2)-1) = (c_Di*pi*AR/((c_L)^2)) - 1;
    end
end
hold on; grid on; grid minor;
plot(taper_ratio, inducedDragFactor)
xlabel("Taper ratio, c_t/c_r");
ylabel("Induced drag factor, \delta");



%             PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%[e,c_L_AR10,Cd] = PLLT(35,2*pi,2*pi,  2,  6,   0,    0,      5,   5,   50);


function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    
    % e: Span efficiency factor (compute)
    % C_L; Coefficient of lift (compute)
    % C_Di; Induced Coefficient of drag (compute)
    % b: span [ft]
    % a0_t: cross-sectional lift slope at the tips [per radian]
    % a0_r: cross-sectional lift slope at the root [per radian]
    % c_t: chord at the tips [ft]
    % c_r: chord at the root [ft]
    % aero_t: zero-lift angle of attack at the tips [deg]
    % aero_r: zero-lift angle of attack at the root [deg]
    % geo_t: geometric angle of attack (geometric twist + alpha) at the tips [deg]
    % geo_r: geometric angle of attack (geometric twist + alpha) at the roots [deg]
    % N: number of odd terms to include in the series expansion for circulation

    aero_t=aero_t*pi/180;
    aero_r=aero_r*pi/180;
    geo_r=geo_r*pi/180;
    geo_t=geo_t*pi/180;
    s = b*(c_r+c_t)/2;
    AR = (b^2)/s;

    % 
    % y = linspace(0, (b/2)-0.0001, N+1);
    % %y = y(1:end-1);
    % theta1 = flip(acos(2*y/b));

    theta = (1:N) * pi / (2*N);
    
    a0 = a0_r + (a0_t-a0_r).*(cos(theta));
    chord = c_r + (c_t-c_r).*(cos(theta));
    aero = aero_r + (aero_t-aero_r).*(cos(theta));
    geo = geo_r + (geo_t-geo_r).*(cos(theta));

    matrix1 = zeros(N,1);
    matrix2 = zeros(N);

    for n=1:N

        matrix1(n,1) = geo(n)-aero(n);

        for m=1:N
            mult = 2*m-1;
            matrix2(n,m) = ( (4*b/a0(n)/chord(n)) * sin(mult*theta(n)) ) + (mult*sin(mult*theta(n))/sin(theta(n)));
        end
    end 

    FoCoeff = matrix2\matrix1;

    for n=2:N
        deltaSum(n,1) = n * (FoCoeff(n)/FoCoeff(1))^2; 
    end
    delta = sum(deltaSum);
    
    
    c_L = FoCoeff(1)*pi*AR;
    e = 1 / (1+delta);
    c_Di = (c_L^2)/pi/e/AR;

end
