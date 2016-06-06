function s = func_EKF(s)
% This function carries out the EKF prediction update and measurement correction
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

% Prediction update state:
s.x = s.A*s.x + s.B*s.u;
s.P_x = s.A * s.P_x * s.A' + s.R_v;

% Apply measurement update if a measurement occurs on this time step
if s.z ~= -1;
   % Linearise measurements
   s.H_x = x_Jacobian(s);
   
   % Measurement update state:
   K = s.P_x*s.H_x'*inv(s.H_x*s.P_x*s.H_x'+s.R_n);
   s.x = s.x + K*(s.z-f(s));
   s.P_x = s.P_x - K*s.H_x*s.P_x;
end

end