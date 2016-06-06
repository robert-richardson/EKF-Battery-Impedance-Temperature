function s = func_DEKF(s)
% This function carries out the DEKF prediction update and measurement correction
% steps for each time step.

% Copyright (c) 2016 by Robert Richardson, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

if isnan(s.x)
   % initialize state estimate from first observation
   if diff(size(s.H_x))
      error('Observation matrix must be square and invertible for state autointialization.');
   end
   s.x = inv(s.H_x)*s.z;
   s.P_x = inv(s.H_x)*s.R_n*inv(s.H_x'); 
end
   
% Prediction update parameter:
s.P_h = s.P_h + s.R_e;

% (Update state matrices):
[s.A, s.B, s.C, s.D] = param_update(s);

% Prediction update state:
s.x = s.A*s.x + s.B*s.u;
s.P_x = s.A * s.P_x * s.A' + s.R_v;

% Apply measurement update if a measurement occurs on this time step
if s.z ~= -1;
   % Linearise measurements
   s.H_x = x_Jacobian(s);
   
   % Measurement update state:
   K_x = s.P_x*s.H_x'*inv(s.H_x*s.P_x*s.H_x'+s.R_n);
   s.x = s.x + K_x*(s.z-f(s));
   s.P_x = s.P_x - K_x*s.H_x*s.P_x;
   
   % (Update parameter matrices):
   s.H_h = h_Jacobian(s);
   
   % Measurement update parameter:
   K_h = s.P_h*s.H_h'*inv(s.H_h*s.P_h*s.H_h'+s.R_r);
   s.h = s.h + K_h*(s.z-f(s));
   s.P_h = s.P_h - K_h*s.H_h*s.P_h;
end

end