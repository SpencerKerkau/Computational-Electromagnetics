clear; close all; clc;

Epsilon_Zero = 8.85*10^(-12);

Rod_Length = 1;
Rod_Radius = 1e-3;
Rod_Sections = 20;

Rod_Potential = 1;

Delta_L = Rod_Length/Rod_Sections;

x_i_values = zeros(1,Rod_Sections);
Phi_Column_Vector = ones(Rod_Sections,1) * Rod_Potential * (4*pi*Epsilon_Zero);
A = zeros(Rod_Sections, Rod_Sections);

for x_i = 1:1:Rod_Sections
    x_i_values(1,x_i) = (Delta_L / 2) + ((x_i - 1) * Delta_L);
end

for row = 1:1:Rod_Sections
    for col = 1:1:Rod_Sections
        
        if row == col
            A(row,col) = 4*pi*Rod_Length*log(Delta_L/Rod_Length);
        else
            A(row,col) = 2*pi*Rod_Length*Delta_L / abs(x_i_values(1,row) - x_i_values(1,col));
        end

    end
end

Rho_s_Vector = A \ Phi_Column_Vector;

stairs(x_i_values, Rho_s_Vector)


