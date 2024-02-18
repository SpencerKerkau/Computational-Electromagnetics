%% Reset Code each time you run
clear; close all; clc;

% Set structure with according to problem
Structure_Width = 20; 
Structure_Height = 10;
Mesh_H = 1;

% Create mesh dimensions
Mesh_Row_Length = Structure_Height / Mesh_H - 1;
Mesh_Col_Length = Structure_Width / Mesh_H - 1;

% Create -4 diagonal identity matrix based on mesh dimensions
Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;

% Call arm star program to find distribution across structure
[Mesh, Column_Vector] = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length);

% Call the eig function to get phi distribution and eigen values
[v, c] = eig(Mesh);

% Find the eigenvalues for modes 11 and 21 from the c vector
C11 = c(Mesh_Row_Length * Mesh_Col_Length, Mesh_Row_Length * Mesh_Col_Length);
C21 = c((Mesh_Row_Length * Mesh_Col_Length)-1, (Mesh_Row_Length * Mesh_Col_Length)-1);

% Speed of light constant
Vo = 299792458;

% Find the cutoff frequency using the equation and waveguide structure.
fc11 = Vo * sqrt(-C11) / (2*pi*Mesh_H)
fc_calculation_11 = (Vo / (2*pi)) * sqrt(((1*pi)/Structure_Width)^2 + ((1*pi)/Structure_Height)^2)

% Find the cutoff frequency using the equation and waveguide structure.
fc21 = Vo * sqrt(-C21) / (2*pi*Mesh_H)
fc_calculation_21 = (Vo / (2*pi)) * sqrt(((2*pi)/Structure_Width)^2 + ((1*pi)/Structure_Height)^2)

% Createa  vector to store the distributions for each mode. You can also
% just use one vector but im lazy
v_column_11 = zeros(Mesh_Row_Length*Mesh_Col_Length, 1);
v_column_21 = zeros(Mesh_Row_Length*Mesh_Col_Length, 1);

% Loop and store all the phi values for each mode
Idx = 1;
for Row = 1:1:Mesh_Row_Length*Mesh_Col_Length
    
    v_column_11(Idx,1) = v(Row, Mesh_Row_Length*Mesh_Col_Length);
    v_column_21(Idx,1) = v(Row, Mesh_Row_Length*Mesh_Col_Length-1);
    
    Idx = Idx + 1;

end

% Create new mesh for formatting purposes. If you dont do this, the color
% plot will look wack
Formatted_Output_Phi_11 = ones(Mesh_Row_Length, Mesh_Col_Length);
Formatted_Output_Phi_21 = ones(Mesh_Row_Length, Mesh_Col_Length);

% Fix the formatting of the phi distributions for each mode
Idx = 1;
for Row = 1:1:Mesh_Row_Length
    for Col = 1:1:Mesh_Col_Length

        Formatted_Output_Phi_11(Row,Col) = v_column_11(Idx,1);
        Formatted_Output_Phi_21(Row,Col) = v_column_21(Idx,1);

        Idx = Idx + 1;

    end
end

% Plot the distributions
figure(1)
imagesc(-Formatted_Output_Phi_11);
colorbar;
title('TM11 Distribution')

figure(2)
imagesc(-Formatted_Output_Phi_21);
colorbar;
title('TM21 Distribution')

%% Finding the electric field vectors

% Create mesh grid to store values for when you call the quiver function
[X, Y] = meshgrid(1:1:Mesh_Col_Length, 1:1:Mesh_Row_Length);

% Call the function to do find the partial derivatives at each point across
% your formatted Phi distribution matrix
[X_Vector, Y_Vector] = Vector_Electric_Field(Formatted_Output_Phi_11, Mesh_Row_Length, Mesh_Col_Length, Mesh_H);

% Create new column vectors to store the answer
Formatted_X_Vector = zeros(Mesh_Row_Length, Mesh_Col_Length); 
Formatted_Y_Vector = zeros(Mesh_Row_Length, Mesh_Col_Length);

% Loop and store the the values of each partial derivative across the
% matrix. We need to do it this way cause the way quiver works is strange.  
Idx = 1;
for Row = 1:1:Mesh_Row_Length
    for Col = 1:1:Mesh_Col_Length

        Formatted_X_Vector(Row,Col) = X_Vector(Idx,1);
        Formatted_Y_Vector(Row,Col) = Y_Vector(Idx,1);
        Idx = Idx + 1;
    end
end

% Plot
figure(3)
quiver(X,Y,Formatted_X_Vector,Formatted_Y_Vector)
set(gca, 'YDir','reverse')
title('TM11 Electric Field Vectors')

%% Repeat the above process but for the 21 mode

[X_Vector, Y_Vector] = Vector_Electric_Field(Formatted_Output_Phi_21, Mesh_Row_Length, Mesh_Col_Length, Mesh_H);


Idx = 1;
for Row = 1:1:Mesh_Row_Length
    for Col = 1:1:Mesh_Col_Length

        Formatted_X_Vector(Row,Col) = X_Vector(Idx,1);
        Formatted_Y_Vector(Row,Col) = Y_Vector(Idx,1);
        Idx = Idx + 1;
    end
end

figure(4)
quiver(X,Y,Formatted_X_Vector,Formatted_Y_Vector)
set(gca, 'YDir','reverse')
title('TM21 Electric Field Vectors')

%% Five Point Arm Star Algorithm

function [Mesh, Column_Vector] = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length)

    Column_Vector = zeros(Mesh_Row_Length * Mesh_Col_Length, 1);
    
    Left_Voltage = 0; Right_Voltage = 0; Top_Voltage = 0; Bottom_Voltage = 0;
    PHI_Matrix_Index = 1;
    
    for Row = 1:1:Mesh_Row_Length
        
        for Col = 1:1:Mesh_Col_Length
            
            Electric_Potential_Sum = 0;
            
            % Right Check
            if Col + 1 <= Mesh_Col_Length
                Mesh(PHI_Matrix_Index,PHI_Matrix_Index+1) = 1;
            else
                Electric_Potential_Sum = Electric_Potential_Sum - Right_Voltage;  
            end
            
            % Left Check. Same process as right check but for left neighbor.
            if Col - 1 == 0
                Electric_Potential_Sum = Electric_Potential_Sum - Left_Voltage;
            else
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1;
            end
    
            % Bottom Check
            if Row + 1 <= Mesh_Row_Length
                Mesh(PHI_Matrix_Index + Mesh_Col_Length, PHI_Matrix_Index) = 1;
            else
                Electric_Potential_Sum = Electric_Potential_Sum - Bottom_Voltage;
            end
    
            % Top Check. Same process as bottom check but for top neighbor.
            if Row - 1 == 0
                Electric_Potential_Sum = Electric_Potential_Sum - Top_Voltage;
            else
                Mesh(PHI_Matrix_Index - Mesh_Col_Length, PHI_Matrix_Index) = 1;
            end
    
            Column_Vector(PHI_Matrix_Index) = Electric_Potential_Sum;
    
            % Move the Phi index to the next phi in the mesh.
            PHI_Matrix_Index = PHI_Matrix_Index + 1;
        end
        
    end
    
end

%% 

function [X_Vector, Y_Vector] = Vector_Electric_Field(Formatted_Output_Phi, Mesh_Row_Length, Mesh_Col_Length, Mesh_H)
    
    Idx = 1;
    X_Vector = [];
    Y_Vector = [];
    
    for Row = 1:1:Mesh_Row_Length
        for Col = 1:1:Mesh_Col_Length

            % X Vector
            if Col == 1 
                X_Vector(Idx,1) = Formatted_Output_Phi(Row, Col + 1) / (2*Mesh_H); %#ok<*AGROW>
            elseif Col == Mesh_Col_Length
                X_Vector(Idx,1) = -Formatted_Output_Phi(Row, Col - 1) / (2*Mesh_H);
            else
                X_Vector(Idx,1) =  (Formatted_Output_Phi(Row,Col + 1) - Formatted_Output_Phi(Row, Col - 1)) / (2*Mesh_H);
            end

            % Y Vector
            if Row == 1 
                Y_Vector(Idx,1) = Formatted_Output_Phi(Row + 1, Col) / (2*Mesh_H);
            elseif Row == Mesh_Row_Length
                Y_Vector(Idx,1) = -Formatted_Output_Phi(Row - 1, Col) / (2*Mesh_H);
            else
                Y_Vector(Idx,1) =  (Formatted_Output_Phi(Row + 1 ,Col) - Formatted_Output_Phi(Row - 1, Col)) / (2*Mesh_H);
            end

            Idx = Idx + 1;

        end
    end

end
