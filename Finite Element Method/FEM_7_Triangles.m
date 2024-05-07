clear; close all; clc;

%% Set up the Dimensions/Coordinates of the problem

rows = 7;
cols = 3;

C = cell(rows, cols);

C{1,1} = [0, 0];
C{1,2} = [0.5, 0];
C{1,3} = [0.5, 0.25];

C{2,1} = [0.5, 0];
C{2,2} = [1, 0];
C{2,3} = [0.5, 0.25];

C{3,1} = [1, 0];
C{3,2} = [1, 0.25];
C{3,3} = [0.5, 0.25];

C{4,1} = [1, 0.25];
C{4,2} = [1, 0.5];
C{4,3} = [0.5, 0.25];

C{5,1} = [1, 0.5];
C{5,2} = [0.5, 0.5];
C{5,3} = [0.5, 0.25];

C{6,1} = [0.5, 0.5];
C{6,2} = [0, 0.5];
C{6,3} = [0.5, 0.25];

C{7,1} = [0, 0.5];
C{7,2} = [0, 0];
C{7,3} = [0.5, 0.25];

%%

E1 = FEM(C, 1);
E2 = FEM(C, 2);
E3 = FEM(C, 3);
E4 = FEM(C, 4);
E5 = FEM(C, 5);
E6 = FEM(C, 6);
E7 = FEM(C, 7);

Global_Element_Matrix = zeros(8,8);

Global_Element_Matrix(1,1) = E1(1,1) + E7(2,2);
Global_Element_Matrix(1,2) = E2(1,2);
Global_Element_Matrix(1,7) = E7(2,1);
Global_Element_Matrix(1,8) = E1(1,3) + E7(2,3);

Global_Element_Matrix(2,1) = E1(1,2);
Global_Element_Matrix(2,2) = E1(2,2) + E2(1,1);
Global_Element_Matrix(2,3) = E2(1,2);
Global_Element_Matrix(2,8) = E1(2,3) + E2(1,3);

Global_Element_Matrix(3,2) = E2(1,2);
Global_Element_Matrix(3,3) = E2(2,2) + E3(1,1);
Global_Element_Matrix(3,4) = E3(1,2);
Global_Element_Matrix(3,8) = E2(2,3) + E3(1,3);

Global_Element_Matrix(4,3) = E3(1,2);
Global_Element_Matrix(4,4) = E3(2,2) + E4(1,1);
Global_Element_Matrix(4,5) = E4(1,2);
Global_Element_Matrix(4,8) = E3(2,3) + E4(1,3);

Global_Element_Matrix(5,4) = E4(2,1);
Global_Element_Matrix(5,5) = E4(2,2) + E5(1,1);
Global_Element_Matrix(5,6) = E5(1,2);
Global_Element_Matrix(5,8) = E4(2,3) + E5(1,3);

Global_Element_Matrix(6,5) = E5(2,1);
Global_Element_Matrix(6,6) = E5(2,2) + E6(1,1);
Global_Element_Matrix(6,7) = E6(1,2);
Global_Element_Matrix(6,8) = E5(2,3) + E6(1,3);

Global_Element_Matrix(7,1) = E7(2,1);
Global_Element_Matrix(7,6) = E6(2,1);
Global_Element_Matrix(7,7) = E6(2,2) + E7(1,1);
Global_Element_Matrix(7,8) = E6(2,3) + E7(1,3);

Global_Element_Matrix(8,1) = E1(3,1) + E7(3,2);
Global_Element_Matrix(8,2) = E1(3,2) + E2(3,1);
Global_Element_Matrix(8,3) = E2(3,2) + E3(3,1);
Global_Element_Matrix(8,4) = E3(3,2) + E4(3,1);
Global_Element_Matrix(8,5) = E4(3,2) + E5(3,1);
Global_Element_Matrix(8,6) = E5(3,2) + E6(3,1);
Global_Element_Matrix(8,7) = E6(3,2) + E7(3,1);
Global_Element_Matrix(8,8) = E1(3,3) + E2(3,3) + E3(3,3) + E4(3,3) + E5(3,3) + E6(3,3) + E7(3,3);

% Calculated Phi gives the same answer as using symbolic method
% Calculated_Phi_8 = (-1/Global_Element_Matrix(8,8)) * ( (Global_Element_Matrix(2,8)*50) + (Global_Element_Matrix(3,8)*75) + (Global_Element_Matrix(4,8) * 100) ...
%     + (Global_Element_Matrix(5,8)*75) + (Global_Element_Matrix(6,8)*50))

syms Phi_8

Phi_Matrix = [0; 0.3; 50; 100; 50; 0.3; 0; Phi_8]; 

L = 0.5 * transpose(Phi_Matrix) * Global_Element_Matrix * Phi_Matrix;

G8 = diff(L,Phi_8)

sol = solve(G8,Phi_8)


%%

function Matrix_Equation = FEM(C, Element)

x1 = C{Element,1}(1);
y1 = C{Element,1}(2);

x2 = C{Element,2}(1);
y2 = C{Element,2}(2);

x3 = C{Element,3}(1);
y3 = C{Element,3}(2);

Area = (0.5)*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));

c12 = 1/(4*Area)*((y2-y3)*(y3-y1)+(x3-x2)*(x1-x3));
c13 = 1/(4*Area)*((y2-y3)*(y1-y2)+(x3-x2)*(x2-x1));
c23 = 1/(4*Area)*((y3-y1)*(y1-y2)+(x1-x3)*(x2-x1));
c11 = 1/(4*Area)*((y2-y3)^2+(x3-x2)^2);
c22 = 1/(4*Area)*((y3-y1)^2+(x1-x3)^2);
c33 = 1/(4*Area)*((y1-y2)^2+(x2-x1)^2);

c21 = c12;
c31 = c13;
c32 = c23;

Matrix_Equation = [c11,c12,c13; c21,c22,c23; c31,c32,c33];

end


