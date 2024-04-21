clear; close all; clc;

% Rits Equation yn = x(1-x)(C1 + C2x + ... + Cnx^n-1)
% Substitute y2 and use base equation twice

syms x y c1 c2

y = x*(1-x)*(c1+c2*x); % First equation
y1 = diff(y);

f = y1^2 - y^2 - 2 * x * y; % Substitute back and use equation 2nd time

I = int(f,0,1); % Boundary condition                       
C1 = diff(I,c1);
C2 = diff(I,c2);

sol = solve([C1, C2], [c1, c2]);

c1 = sol.c1;
c2 = sol.c2;

x_val = 3/4;
y_val = (sin(x_val)/sin(1)) - x_val;
y1_val = (5/18)*x_val*(1-x_val);
y2_val = x_val*(1-x_val)*(c1+c2*x_val);

output_vector = [x_val, y_val, y1_val, y2_val];
fprintf('Output Vector: [%f, %f, %f, %f]\n', output_vector);
