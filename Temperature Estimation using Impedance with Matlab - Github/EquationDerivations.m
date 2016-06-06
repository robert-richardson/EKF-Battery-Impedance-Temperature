%%
%-------------------------------------------------------------------------%
% Derivation of equations
%-------------------------------------------------------------------------%
%
%{
This code steps through the derivations of the equations used in the
mainscript functions, using Matlab Symbolic Math.

Note that a published html document showing this code and its results is
located in the 'html' sub-folder within this package.

Usage:
(i)   Run the program
(ii)  The program pauses at different steps in the derivation. Press any 
      key to continue to the next derivation.

I would ask that you cite this paper as Richardson, Robert R., and 
David A. Howey. "Sensorless battery internal temperature estimation using a
kalman filter with impedance measurement." Sustainable Energy, IEEE 
Transactions on 6.4 (2015): 1190-1199. if you want to use this code for 
your own research. For further details on the work of the Energy Power 
Group at Oxford, please see epg.eng.ox.ac.uk.

Copyright (c) 2016 by Robert Richardson, David Howey
and The Chancellor, Masters and Scholars of the University of Oxford.
See the licence file LICENCE.txt for more information.

%}


%%
%-------------------------------------------------------------------------%
% Clear environemnt and define symbolic math
%-------------------------------------------------------------------------%
%
% Clear environment
clear; close all; clc

% Define initial symbols
syms a b d Ts Tbar gam R

a = 4*Ts - 3*Tbar - 15*R*gam/8;
b = -18*Ts + 18*Tbar + 15*R*gam/2;
d = 15*Ts - 15*Tbar - 45*R*gam/8;


%% 
%-------------------------------------------------------------------------%
% Derivation of admittance
%-------------------------------------------------------------------------%
% (as used in function 'f.m')

syms a3 a2 a1
admittance = a1*a^2 + a1*a*b + (2*a1*a*d)/3 + a2*a + (a1*b^2)/3 + (a1*b*d)/2 +...
    (a2*b)/2 + (a1*d^2)/5 + (a2*d)/3 + a3;
admittance = simplify(admittance);

% Display to screen
fprintf('admittance (intermediate step):\n\n')
disp(admittance)

% admittance = (15*a1*R^2*gam^2)/32 + (15*a1*R*Tbar*gam)/8 -...
%     (15*a1*R*Ts*gam)/8 + 3*a1*Tbar^2 -...
%     4*a1*Tbar*Ts + a2*Tbar + 2*a1*Ts^2 + a3

disp('press any key to continue...')
pause


% Now substitute for Ts from eq. ()
% where:
% Ts = (C21*Tbar + C22*gam + D22*Tinf)
% where:
% Cij and Dij are the (i,j) elements of matrices C and D

syms Tinf C21 C22 D22
admittance = (15*a1*R^2*gam^2)/32 + (15*a1*R*Tbar*gam)/8 -...
    (15*a1*R*(C21*Tbar + C22*gam + D22*Tinf)*gam)/8 + 3*a1*Tbar^2 -...
    4*a1*Tbar*(C21*Tbar + C22*gam + D22*Tinf) +...
    a2*Tbar + 2*a1*(C21*Tbar + C22*gam + D22*Tinf)^2 + a3;

% Display to screen
fprintf('\n\nadmittance (as used in function f.m):\n\n')
disp(admittance)

% admittance = a3 + Tbar*a2 + 3*Tbar^2*a1 +...
%     2*a1*(C22*gam + C21*Tbar + D22*Tinf)^2 +...
%     (15*R^2*a1*gam^2)/32 -...
%     4*Tbar*a1*(C22*gam + C21*Tbar + D22*Tinf) +...
%     (15*R*Tbar*a1*gam)/8 -...
%     (15*R*a1*gam*(C22*gam + C21*Tbar + D22*Tinf))/8;

disp('press any key to continue...')
pause


%% 
%-------------------------------------------------------------------------%
% Derivation of state Jacobian
%-------------------------------------------------------------------------%
% (as used in function 'x_Jacobian.m')
% 
% H_x(1) = d(Z)/d(x1)
% H_x(2) = d(Z)/d(x2)

Z = 1/admittance;

H_x_1 = diff(Z,Tbar);
fprintf('\n\n H_x(1) (as used in function x_Jacobian.m):\n\n')
disp(H_x_1)

H_x_2 = diff(Z,gam);
fprintf('\n\nH_x(2) (as used in function x_Jacobian.m):\n\n')
disp(H_x_2)


% H_x(1) = -(a2 + 6*Tbar_ss*a1 - 4*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf) - ...
%     4*C21*Tbar_ss*a1 + 4*C21*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf) + ...
%     (15*R*gam_ss*a1)/8 - (15*C21*R*gam_ss*a1)/8)/(a3 + Tbar_ss*a2 + ...
%     3*Tbar_ss^2*a1 + 2*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf)^2 + ...
%     (15*R^2*gam_ss^2*a1)/32 - 4*Tbar_ss*a1*(C22*gam_ss + ...
%     C21*Tbar_ss + D22*T_inf) + (15*R*Tbar_ss*gam_ss*a1)/8 - ...
%     (15*R*gam_ss*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf))/8)^2;
% 
% H_x(2) = (4*C22*Tbar_ss*a1 - 4*C22*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf) -...
%     (15*R*Tbar_ss*a1)/8 + (15*R*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf))/8 -...
%     (15*R^2*gam_ss*a1)/16 + (15*C22*R*gam_ss*a1)/8)/...
%     (a3 + Tbar_ss*a2 + 3*Tbar_ss^2*a1 + ...
%     2*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf)^2 + ...
%     (15*R^2*gam_ss^2*a1)/32 - 4*Tbar_ss*a1*(C22*gam_ss + ...
%     C21*Tbar_ss + D22*T_inf) + (15*R*Tbar_ss*gam_ss*a1)/8 - ...
%     (15*R*gam_ss*a1*(C22*gam_ss + C21*Tbar_ss + D22*T_inf))/8)^2;

disp('press any key to continue...')
pause

%% 
%-------------------------------------------------------------------------%
% Derivation of parameter Jacobian
%-------------------------------------------------------------------------%
% (as used in function 'h_Jacobian.m')
% d(Z)/d(h)


% Manually substitute for C(2,1), C(2,2) and D(2,2) into admittance to give
% admittance in terms of the base parameters
% where:
% C21 = (24*kt/(24*kt + R*h));
% C22 = (15*R*kt/(48*kt + 2*R*h));
% D22 = (R*h / (24*kt + R*h));
syms alpha kt h

admittance_base_parameters = (15*a1*R^2*gam^2)/32 + (15*a1*R*Tbar*gam)/8 -...
    (15*a1*R*((24*kt/(24*kt + R*h))*Tbar + (15*R*kt/(48*kt + 2*R*h))*gam + (R*h / (24*kt + R*h))*Tinf)*gam)/8 + 3*a1*Tbar^2 -...
    4*a1*Tbar*((24*kt/(24*kt + R*h))*Tbar + (15*R*kt/(48*kt + 2*R*h))*gam + (R*h / (24*kt + R*h))*Tinf) +...
    a2*Tbar + 2*a1*((24*kt/(24*kt + R*h))*Tbar + (15*R*kt/(48*kt + 2*R*h))*gam + (R*h / (24*kt + R*h))*Tinf)^2 + a3;

Z_base_parameters = 1/admittance_base_parameters;
H_h = diff(Z_base_parameters,h);

fprintf('\n\nH_h (as used in function h_Jacobian.m):\n\n')
disp(H_h)

% H_h = -(4*Tbar*a1*((24*R*Tbar*kt)/(24*kt + R*h)^2 -                 ...
%     (R*T_inf)/(24*kt + R*h) + (R^2*T_inf*h)/(24*kt + R*h)^2 +       ...
%     (30*R^2*gam*kt)/(48*kt + 2*R*h)^2) -                            ...
%     4*a1*((24*Tbar*kt)/(24*kt + R*h) + (R*T_inf*h)/(24*kt + R*h) +  ...
%     (15*R*gam*kt)/(48*kt + 2*R*h))*((24*R*Tbar*kt)/(24*kt + R*h)^2 -...
%     (R*T_inf)/(24*kt + R*h) + (R^2*T_inf*h)/(24*kt + R*h)^2 +       ...
%     (30*R^2*gam*kt)/(48*kt + 2*R*h)^2) +                            ...
%     (15*R*gam*a1*((24*R*Tbar*kt)/(24*kt + R*h)^2 -                  ...
%     (R*T_inf)/(24*kt + R*h) + (R^2*T_inf*h)/(24*kt + R*h)^2 +       ...
%     (30*R^2*gam*kt)/(48*kt + 2*R*h)^2))/8)/                         ...
%     (a3 + Tbar*a2 + 3*Tbar^2*a1 + 2*a1*((24*Tbar*kt)/(24*kt + R*h) +...
%     (R*T_inf*h)/(24*kt + R*h) + (15*R*gam*kt)/(48*kt + 2*R*h))^2 +  ...
%     (15*R^2*gam^2*a1)/32 - 4*Tbar*a1*((24*Tbar*kt)/(24*kt + R*h) +  ...
%     (R*T_inf*h)/(24*kt + R*h) + (15*R*gam*kt)/(48*kt + 2*R*h)) +    ...
%     (15*R*Tbar*gam*a1)/8 - (15*R*gam*a1*((24*Tbar*kt)/(24*kt + R*h)+...
%     (R*T_inf*h)/(24*kt + R*h) + (15*R*gam*kt)/(48*kt + 2*R*h)))/8)^2;



