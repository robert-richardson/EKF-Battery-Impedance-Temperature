function H_x = x_Jacobian(s)
% This function calculates H_x, Jacobian matrix of partial derivatives of f
% with respect to x for each time step.

% Copyright (c) 2016 by Robert Richardson, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

persistent a3 a2 a1 R T_inf
persistent firstRun

if isempty(firstRun)
    R = s.r_o;
    T_inf =  s.T_inf;
    a1 = s.a1;
    a2 = s.a2;
    a3 = s.a3;    
    firstRun = 1;
end

Tbar_ss = s.x(1);
gam_ss = s.x(2);
C21 = s.C(2,1);
C22 = s.C(2,2);
D22 = s.D(2,2);

% H_x(1) =>d(Z)/d(x1)
H_x(1) = -(a2 + 6*Tbar_ss*a1 - 4*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf) - ...
    4*C21*Tbar_ss*a1 + 4*C21*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf) + ...
    (15*R*gam_ss*a1)/8 - (15*C21*R*gam_ss*a1)/8)/(a3 + Tbar_ss*a2 + ...
    3*Tbar_ss^2*a1 + 2*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf)^2 + ...
    (15*R^2*gam_ss^2*a1)/32 - 4*Tbar_ss*a1*(C22*gam_ss + ...
    C21*Tbar_ss + D22*T_inf) + (15*R*Tbar_ss*gam_ss*a1)/8 - ...
    (15*R*gam_ss*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf))/8)^2;
% H_x(1) = -(a2 + 2*Tbar_ss*a1)/(a1*Tbar_ss^2 + a2*Tbar_ss + a3)^2; % linear approx

% H_x(1) =>d(Z)/d(x2)
H_x(2) = (4*C22*Tbar_ss*a1 - 4*C22*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf) -...
    (15*R*Tbar_ss*a1)/8 + (15*R*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf))/8 -...
    (15*R^2*gam_ss*a1)/16 + (15*C22*R*gam_ss*a1)/8)/...
    (a3 + Tbar_ss*a2 + 3*Tbar_ss^2*a1 + ...
    2*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf)^2 + ...
    (15*R^2*gam_ss^2*a1)/32 - 4*Tbar_ss*a1*(C22*gam_ss + ...
    C21*Tbar_ss + D22*T_inf) + (15*R*Tbar_ss*gam_ss*a1)/8 - ...
    (15*R*gam_ss*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf))/8)^2;
% H_x(2)=0; % linear approx
end
