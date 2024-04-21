clear; close all; clc;

Epsilon_Zero = 8.854 * 10^-12;
Mu_Zero = 4*pi * 10^-7;
N_Zero = sqrt(Mu_Zero / Epsilon_Zero);
C = 3 * 10^8;

Frequency = 1 * 10^9;

Lambda = C / Frequency;
B_Zero = 2*pi / Lambda;

Vg = 1;

Antenna_Lambda = 1;

Length = Antenna_Lambda * Lambda;
Wire_Radius = 0.0001 * Lambda;

Section_Count = 100;

Section_Length = Length / Section_Count;

Z = zeros(1, Section_Count);

for Idx = 1:1:Section_Count
    Z(1, Idx) = (-Length / 2) + (Section_Length / 2) + (Idx - 1) * Section_Length;
end

Z(Section_Count + 1) = (Length/2) - 2 * Section_Length;

R_Distance_Matrix = zeros(Section_Count + 1, Section_Count + 1);
A = zeros(Section_Count + 1, Section_Count + 1);

for R_X = 1:1:Section_Count+1
    for R_Y = 1:1:Section_Count+1
        R_Distance_Matrix(R_X, R_Y) = sqrt(Wire_Radius^2 + (Z(R_Y) - Z(R_X))^2);
    end
end

for R_X = 1:1:Section_Count+1
    
    for R_Y = 1:1:Section_Count
        A(R_X, R_Y) = exp(-1i * B_Zero * R_Distance_Matrix(R_X, R_Y)) / R_Distance_Matrix(R_X, R_Y);
    end

    A(R_X, R_Y + 1) = cos(B_Zero * Z(R_X));

end

B = transpose(-1i*Vg/(2*N_Zero)*sin(B_Zero*abs(Z)));
A = A/(4*pi);

X = A\B;

% V = IR -> R = V / I
% Change current to only the center current
Impedance = Vg/X(floor(length(A)/2))

plot(abs(X(1:Section_Count)),Z(1:Section_Count))

