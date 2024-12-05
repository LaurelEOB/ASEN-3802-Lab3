%% PPLT Function
function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    %   PLLT is a function to solve the fundamental equation of Prandtl Lifting
    %   Line Theory for finite wings with thick airfoils. 
    %       Authors: Laurel O'Brien, Gabriel Burdan, Brody Ambroggio, William Wallingford
    %       Collaborators: Alyxis Ellington, Samantha Sheppard, Professor Hoke
    %       Date Last Revised: 12-5-2024

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

    % Converting to rad
    aero_t=aero_t*pi/180;
    aero_r=aero_r*pi/180;
    geo_r=geo_r*pi/180;
    geo_t=geo_t*pi/180;
    % Constants
    s = b*(c_r+c_t)/2;
    AR = (b^2)/s;

    % Theta coordinates along wing
    theta = (1:N) * pi / (2*N);
    
    % Linearly varing properties
    a0 = a0_r + (a0_t-a0_r).*(cos(theta));
    chord = c_r + (c_t-c_r).*(cos(theta));
    aero = aero_r + (aero_t-aero_r).*(cos(theta));
    geo = geo_r + (geo_t-geo_r).*(cos(theta));

    % Preallocating 
    matrix1 = zeros(N,1);
    deltaSum = zeros(N,1);
    matrix2 = zeros(N);

    % Finding linear system of equations
    for n=1:N

        matrix1(n,1) = geo(n)-aero(n);

        for m=1:N
            oddTerms = 2*m-1; % Only using odd terms
            matrix2(n,m) = ( (4*b/a0(n)/chord(n)) * sin(oddTerms*theta(n)) ) + (oddTerms*sin(oddTerms*theta(n))/sin(theta(n)));
        end
    end 

    % Fourier Coefficients An
    FoCoeff = matrix2\matrix1;

    % Delta calculation
    for n=2:N
        oddTerms = 2*n-1;
        deltaSum(n,1) = oddTerms * (FoCoeff(n)/FoCoeff(1))^2; 
    end
    delta = sum(deltaSum);
    
    % Final values
    c_L = FoCoeff(1)*pi*AR;
    e = 1 / (1+delta);
    c_Di = (c_L^2)/pi/e/AR;

end