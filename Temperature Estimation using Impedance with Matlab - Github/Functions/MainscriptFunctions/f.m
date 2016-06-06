function y = f(s)
% This function calculates f(x, h), the non-linear function relating the
% state vector to the measurement (eq. 9)

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


admittance = a3 + Tbar_ss*a2 + 3*Tbar_ss^2*a1 +...
    2*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf)^2 +...
    (15*R^2*gam_ss^2*a1)/32 -...
    4*Tbar_ss*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf) +...
    (15*R*Tbar_ss*gam_ss*a1)/8 -...
    (15*R*gam_ss*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf))/8;

% admittance = a3 + Tbar_ss*a2 + (Tbar_ss^2)*a1; % (linear approx.)    
y = 1/admittance;
end