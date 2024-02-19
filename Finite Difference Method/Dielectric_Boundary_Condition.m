%% Reset Code each time you run
clear; close all; clc;

%% Specify problem constraints according to IEEE Table I, third example from left
% a = 2.02, b = 7.0, h = 1.00, w = 1.00, t = 0.01, e1 = e0, e2 = 9.6*e0

% Set structure with according to problem
Structure_Width = 7; 
Structure_Height = 2;
Mesh_H = 0.1;

% Conductor Values
Conductor_Width = 1;
Conductor_Height = 1;
Conductor_Voltage = 10;

% Set permittivity values
Epsilon_Zero = 8.85*10^(-12);
Epsilon_One = 9.6 * Epsilon_Zero;

% Create mesh dimensions
Mesh_Row_Length = Structure_Height / Mesh_H - 1;
Mesh_Col_Length = Structure_Width / Mesh_H - 1;

% Create -4 diagonal identity matrix based on mesh dimensions
Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;

%% Specify Dielectric Location & Boundary Condition

% Calculate row that the boundary and conductor sit on within the mesh
Dielectric_Boundary_Row = (Structure_Height / Mesh_H) - (Conductor_Height/Mesh_H);
Conductor_Row = Dielectric_Boundary_Row;

% Find the left and right PHI on the boundary that the conductor takes up.
% Width of the conductor was set to 7 column spaces.

Conductor_Width_Left_Node = ceil(Mesh_Col_Length/2) - (Conductor_Width/Mesh_H/2);
Conductor_Width_Right_Node = Conductor_Width_Left_Node + (Conductor_Width/Mesh_H) - 2;

%% Develop the contour around the center conductor

% Left, Right, Top, Bottom wall of the contour calculated based on taking
% half of the distance from the center of the conductor to the wall of the
% mesh. 


% DEFAULT CONTOUR ADJUSTMENT = 4, YOU CAN ONLY ADD TO THE ADJUSTMENT. DO
% NOT GO UNDER 2.5 OR IT WILL BREAK. WHEN YOU ADD TO THE
% ADJUSTMENT, THE IMPEDANCE SHOWS AS THE SAME BUT ITS ACTUALLY VERY
% SLIGHTLY DIFFERENT.
Contour_Adjustment = 4;
Contour_Right_Wall = floor(Mesh_Col_Length/2) + ceil(Mesh_Col_Length/Contour_Adjustment);
Contour_Left_Wall = floor(Mesh_Col_Length/2) - ceil(Mesh_Col_Length/Contour_Adjustment);
Contour_Top_Wall = Dielectric_Boundary_Row - ceil(Mesh_Row_Length/Contour_Adjustment);
Contour_Bottom_Wall = Dielectric_Boundary_Row + ceil(Mesh_Row_Length/Contour_Adjustment);

%% Create column vector to store potentials

% Create the column vector and set all values to zero
Column_Vector = zeros(Mesh_Row_Length * Mesh_Col_Length, 1);

% Left and right node store where the conductor voltage should be within
% the column vector since the conductor is centered within the mesh.
Left_Node = (Dielectric_Boundary_Row - 1) * Mesh_Col_Length + (Conductor_Width_Left_Node);
Right_Node = (Dielectric_Boundary_Row - 1) * Mesh_Col_Length + (Conductor_Width_Right_Node);
Column_Vector(Left_Node:Right_Node, 1) = Conductor_Voltage;

%% Create the Phi mesh distribution given the values of epsilon for epsilon zero and epsilon one

% Call arm star function and pass through mesh characteristics, conductor
% location, and permittivity values.
Mesh_With_Epsilon = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, Dielectric_Boundary_Row, ...
    Conductor_Width_Left_Node, Conductor_Width_Right_Node, Epsilon_Zero, Epsilon_One);   

% Note that the "Mesh_Without_Epsilon" is the same calculation but the
% entire mesh has epsilon zero permittivity. 
Mesh_Without_Epsilon = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, Dielectric_Boundary_Row, ...
    Conductor_Width_Left_Node, Conductor_Width_Right_Node, Epsilon_Zero, Epsilon_Zero);  

%% Matrix Multiplication to Solve for Phi

% Pass both mesh's to the "Solve_For_Phi" function to find the phi
% distribution values and reformat to row x column format

Formatted_Output_Phi_With_Epsilon = Solve_For_Phi(Mesh_With_Epsilon, Column_Vector, Mesh_Row_Length, Mesh_Col_Length, ...
    Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall);

% Plot the distribution using the color format
figure(1)
imagesc(Formatted_Output_Phi_With_Epsilon);
colorbar;
title('Phi Distribution With Dielectric')

% Add the contour to the color plot
rectangle('Position',[Contour_Left_Wall, Contour_Top_Wall, ...
    Contour_Right_Wall - Contour_Left_Wall, Contour_Bottom_Wall - Contour_Top_Wall],'EdgeColor','r')

Formatted_Output_Phi_Without_Epsilon = Solve_For_Phi(Mesh_Without_Epsilon, Column_Vector, Mesh_Row_Length, Mesh_Col_Length, ...
    Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall);

figure(2)
imagesc(Formatted_Output_Phi_Without_Epsilon);
colorbar;
title('Phi Distribution Without Dielectric')

% Add the contour to the color plot
rectangle('Position',[Contour_Left_Wall, Contour_Top_Wall, ...
    Contour_Right_Wall - Contour_Left_Wall, Contour_Bottom_Wall - Contour_Top_Wall],'EdgeColor','r')

%% Calculate capacitance from contour

% Pass both formatted phi distributions to the "Calculate_Contour" function
% to find the capacitance of each mesh

q = Calculate_Contour(Formatted_Output_Phi_With_Epsilon, Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall, ...
    Epsilon_Zero, Epsilon_One, Dielectric_Boundary_Row);

q_zero = Calculate_Contour(Formatted_Output_Phi_Without_Epsilon, Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall, ...
    Epsilon_Zero, Epsilon_Zero, Dielectric_Boundary_Row);

%% Find the capacitance and characteristic impedance values

% Speed of light = 299792458

q
c = q / Conductor_Voltage

q_zero
c_zero = q_zero / Conductor_Voltage

Z_0 = 1 / (299792458 * sqrt(c * c_zero))

%% Arm Star Program

function Mesh = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, Dielectric_Boundary_Row, ...
    Conductor_Width_Left_Node, Conductor_Width_Right_Node, Epsilon_Zero, Epsilon_One)

    % PHI_Matrix_Index used to track the phi position within the mesh
    PHI_Matrix_Index = 1;

    for Row = 1:1:Mesh_Row_Length
    
        for Col = 1:1:Mesh_Col_Length
    
            % Right Check
            if Col + 1 <= Mesh_Col_Length
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index + 1) = 1;
            end
            
            % Left Check. Same process as right check but for left neighbor.
            if Col - 1 > 0
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1;
            end
    
            % Bottom Check
            if Row + 1 <= Mesh_Row_Length
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index + Mesh_Col_Length) = 1;
            end
    
            % Top Check. Same process as bottom check but for top neighbor.
            if Row - 1 > 0 
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index - Mesh_Col_Length) = 1;
            end
    
            % Boundary condition
            if Row == Dielectric_Boundary_Row
                
                % Phi position is on the conductor
                if Col >= Conductor_Width_Left_Node && Col <= Conductor_Width_Right_Node
    
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index) = 1;
                    
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index + 1) = 1 * (Epsilon_Zero + Epsilon_One);
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1 * (Epsilon_Zero + Epsilon_One);
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index + Mesh_Col_Length) = 1 * (2 * Epsilon_One);
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index - Mesh_Col_Length) = 1 * (2 * Epsilon_Zero);
    
                else
                    
                    % Phi position on the boundary but not on the conductor
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index) = -4 * (Epsilon_Zero + Epsilon_One);
    
                    if Col + 1 <= Mesh_Col_Length
                        Mesh(PHI_Matrix_Index, PHI_Matrix_Index + 1) = 1 * (Epsilon_Zero + Epsilon_One);
                    end
    
                    if Col - 1 > 0
                        Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1 * (Epsilon_Zero + Epsilon_One);
                    end
                    
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index + Mesh_Col_Length) = 1 * (2 * Epsilon_One);
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index - Mesh_Col_Length) = 1 * (2 * Epsilon_Zero);
    
                end
                    
            end
            
            % Move the Phi index to the next phi in the mesh.
            PHI_Matrix_Index = PHI_Matrix_Index + 1;
    
        end
        
    end
    
end

%% Solve for Phi

function Formatted_Output_Phi = Solve_For_Phi(Mesh, Column_Vector, Mesh_Row_Length, Mesh_Col_Length, ...
    Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall)

    % Solve for the potential at each phi within the mesh.
    Output_Phi = inv(Mesh) * Column_Vector; %#ok<MINV>
    
    % Create new mesh for formatting purposes
    Formatted_Output_Phi = ones(Mesh_Row_Length, Mesh_Col_Length);
    Idx = 1;
    
    % Use to for loops to fill out the formatted mesh. I'm not sure why but
    % using the other method messes with the output. Maybe my fault but idk
    % why it doesnt work. This is fail safe.
    for Row = 1:1:Mesh_Row_Length
        for Col = 1:1:Mesh_Col_Length
            Formatted_Output_Phi(Row,Col) = Output_Phi(Idx,1);
            Idx = Idx + 1;
        end
    end
    
end

%% Contour Potential Summation

function q = Calculate_Contour(Output_Phi_Matrix, Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall, ...
    Epsilon_Zero, Epsilon_One, Dielectric_Boundary_Row)

    q = 0;
    
    for Row = Contour_Top_Wall:1:Contour_Bottom_Wall
        for Col = Contour_Left_Wall:1:Contour_Right_Wall
    
            % Corners of contour
            if Row == Contour_Top_Wall && Col == Contour_Left_Wall % Top Left
                q = q + abs((Epsilon_Zero * (((Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2) + ...
                    (Output_Phi_Matrix(Row-1,Col) - Output_Phi_Matrix(Row+1,Col))/2))/2);
            elseif Row == Contour_Top_Wall && Col == Contour_Right_Wall % Top Right
                q = q + abs((Epsilon_Zero * (((Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2) + ...
                    (Output_Phi_Matrix(Row-1,Col) - Output_Phi_Matrix(Row+1,Col))/2))/2);
            elseif Row == Contour_Bottom_Wall && Col == Contour_Left_Wall % Bottom Left
                q = q + abs((Epsilon_One * (((Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2) + ...
                    (Output_Phi_Matrix(Row+1,Col) - Output_Phi_Matrix(Row-1,Col))/2))/2);
            elseif Row == Contour_Bottom_Wall && Col == Contour_Right_Wall % Bottom Right
                q = q + abs((Epsilon_One * (((Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2) + ...
                    (Output_Phi_Matrix(Row+1,Col) - Output_Phi_Matrix(Row-1,Col))/2))/2);
            end
    
            % Top wall without corners
            if Row == Contour_Top_Wall && Col > Contour_Left_Wall && Col < Contour_Right_Wall
                q = q + abs((Epsilon_Zero * (Output_Phi_Matrix(Row-1,Col) - Output_Phi_Matrix(Row+1,Col))/2));
            end
            
            % Bottom wall without corners
            if Row == Contour_Top_Wall && Col > Contour_Left_Wall && Col < Contour_Right_Wall
                q = q + abs((Epsilon_One * (Output_Phi_Matrix(Row+1,Col) - Output_Phi_Matrix(Row-1,Col))/2));
            end
    
            % Left wall without corners
            if Col == Contour_Left_Wall && Row > Contour_Top_Wall && Row < Contour_Bottom_Wall
    
                if Row < Dielectric_Boundary_Row
                    q = q + abs((Epsilon_Zero * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2));
                elseif Row > Dielectric_Boundary_Row
                    q = q + abs((Epsilon_One * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2));
                else
                    q = q + abs(((Epsilon_Zero * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/4) + ...
                        (Epsilon_One * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/4)));
                end
    
            end
    
            % Right wall without corners
            if Col == Contour_Right_Wall && Row > Contour_Top_Wall && Row < Contour_Bottom_Wall
    
                if Row < Dielectric_Boundary_Row
                    q = q + abs((Epsilon_Zero * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2));
                elseif Row > Dielectric_Boundary_Row
                    q = q + abs((Epsilon_One * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2));
                else
                    q = q + abs(((Epsilon_Zero * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/4) + ...
                        (Epsilon_One * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/4)));
                end
    
            end
           
        end
    end

end
