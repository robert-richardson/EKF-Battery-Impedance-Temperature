function [A, B, C, D] = param_update(s)
% This function updates the values of the state matrices at each time step
% using the current estimate of the convection coefficient, h.

% Copyright (c) 2016 by Robert Richardson, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

persistent al kt R Vb deltat
persistent firstRun

if  isempty(firstRun)
    al = s.al; kt = s.kt; R = s.r_o; Vb = s.Vb;
    deltat = s.delta_t;
    firstRun = 1;
end

h = s.h;

A = eye(2) + deltat*[-48*al*h/(R*(24*kt+R*h)),      -15*al*h / (24*kt + R*h) ;                  ...
    -320*al*h /((R^2)*(24*kt+R*h)),                 -120*al*(4*kt + R*h) / ((R^2)*(24*kt+R*h))];
B = deltat*[al/(kt*Vb),                             48*al*h / (R*(24*kt + R*h));                ...
    0,                                              320*al*h / ((R^2)*(24*kt+R*h))];
C = [(24*kt - 3*R*h)/(24*kt + R*h),                 -(120*R*kt+15*(R^2)*h)/(8*(24*kt + R*h));   ...
    24*kt/(24*kt + R*h),                            15*R*kt/(48*kt + 2*R*h)];
D = [0 ,                                            4*R*h / (24*kt + R*h);                      ...
    0 ,                                             R*h / (24*kt + R*h)];

end