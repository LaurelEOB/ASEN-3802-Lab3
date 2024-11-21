function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_r,N)
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


end
