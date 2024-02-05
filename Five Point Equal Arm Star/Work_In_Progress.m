%% Reset Code

clear; close all; clc;

% Set the input provided by the user to individual varibles
Structure_Width = 7; 
Structure_Height = 2;
Mesh_H = 0.1;

Mesh_Row_Length = Structure_Height / Mesh_H - 1;
Mesh_Col_Length = Structure_Width / Mesh_H - 1;

Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;
Column_Vector = zeros(Mesh_Row_Length * Mesh_Col_Length, 1);

%% Specify Dielectric Location & Boundary Condition

Dielectric_Boundary_Row = ceil(Mesh_Row_Length / 2);

Epsilon_Zero = 8.85*10^(-12);
Epsilon_One = 9.6 * Epsilon_Zero;

Conductor_Row = Dielectric_Boundary_Row;
Conductor_Width_Left_Node = floor(Mesh_Col_Length/2) - 3;
Conductor_Width_Right_Node = ceil(Mesh_Col_Length/2) + 3;

Conductor_Voltage = 10;

%% Develop the contour

Contour_Right_Wall = floor(Mesh_Col_Length/2) + ceil(Mesh_Col_Length/4);
Contour_Left_Wall = floor(Mesh_Col_Length/2) - ceil(Mesh_Col_Length/4);
Contour_Top_Wall = Dielectric_Boundary_Row - ceil(Mesh_Row_Length/4);
Contour_Bottom_Wall = Dielectric_Boundary_Row + ceil(Mesh_Row_Length/4);

Mesh_With_Epsilon = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, Dielectric_Boundary_Row, ...
    Conductor_Width_Left_Node, Conductor_Width_Right_Node, Epsilon_Zero, Epsilon_One);   

Mesh_Without_Epsilon = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, Dielectric_Boundary_Row, ...
    Conductor_Width_Left_Node, Conductor_Width_Right_Node, Epsilon_Zero, Epsilon_Zero);  

%% Column Vector Voltage

Left_Node = (Dielectric_Boundary_Row - 1) * Mesh_Col_Length + (Conductor_Width_Left_Node);
Right_Node = (Dielectric_Boundary_Row - 1) * Mesh_Col_Length + (Conductor_Width_Right_Node);
Column_Vector(Left_Node:Right_Node, 1) = Conductor_Voltage;

%% Matrix Multiplication to Solve for Phi

Formatted_Output_Phi_With_Epsilon = Solve_For_Phi(Mesh_With_Epsilon, Column_Vector, Mesh_Row_Length, Mesh_Col_Length, ...
    Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall);

Formatted_Output_Phi_Without_Epsilon = Solve_For_Phi(Mesh_Without_Epsilon, Column_Vector, Mesh_Row_Length, Mesh_Col_Length, ...
    Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall);

%% Calculate capacitance from contour

q = Calculate_Contour(Formatted_Output_Phi_With_Epsilon, Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall, ...
    Epsilon_Zero, Epsilon_One, Dielectric_Boundary_Row);

q_zero = Calculate_Contour(Formatted_Output_Phi_Without_Epsilon, Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall, ...
    Epsilon_Zero, Epsilon_One, Dielectric_Boundary_Row);

%% Charge and Capacitance value with dielectric

q
c = q / Conductor_Voltage

q_zero
c_zero = q_zero / Conductor_Voltage

Z_0 = 1 / (299792458 * sqrt(c * c_zero))

%% Arm Star Program

function Mesh = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, Dielectric_Boundary_Row, ...
    Conductor_Width_Left_Node, Conductor_Width_Right_Node, Epsilon_Zero, Epsilon_One)

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
    
            if Row == Dielectric_Boundary_Row
    
                if Col >= Conductor_Width_Left_Node && Col <= Conductor_Width_Right_Node
    
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index) = 1;
                    
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index + 1) = 1 * (Epsilon_Zero + Epsilon_One);
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1 * (Epsilon_Zero + Epsilon_One);
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index + Mesh_Col_Length) = 1 * (2 * Epsilon_One);
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index - Mesh_Col_Length) = 1 * (2 * Epsilon_Zero);
    
                else
    
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
    
    Formatted_Output_Phi = ones(Mesh_Row_Length, Mesh_Col_Length);
    Idx = 1;
    
    for Row = 1:1:Mesh_Row_Length
        for Col = 1:1:Mesh_Col_Length
            Formatted_Output_Phi(Row,Col) = Output_Phi(Idx,1);
            Idx = Idx + 1;
        end
    end
    
    imagesc(Formatted_Output_Phi);
    colorbar;

    rectangle('Position',[Contour_Left_Wall, Contour_Top_Wall, ...
        Contour_Right_Wall - Contour_Left_Wall, Contour_Bottom_Wall - Contour_Top_Wall],'EdgeColor','r')

end

%% Contour Potential Summation

function q = Calculate_Contour(Output_Phi_Matrix, Contour_Top_Wall, Contour_Bottom_Wall, Contour_Left_Wall, Contour_Right_Wall, ...
    Epsilon_Zero, Epsilon_One, Dielectric_Boundary_Row)

    q = 0;
    
    for Row = Contour_Top_Wall:1:Contour_Bottom_Wall
        for Col = Contour_Left_Wall:1:Contour_Right_Wall
    
            % Corners of contour
            if Row == Contour_Top_Wall && Col == Contour_Left_Wall % Top Left
                q = q + (Epsilon_Zero * (((Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2) + ...
                    (Output_Phi_Matrix(Row-1,Col) - Output_Phi_Matrix(Row+1,Col))/2))/2;
            elseif Row == Contour_Top_Wall && Col == Contour_Right_Wall % Top Right
                q = q + (Epsilon_Zero * (((Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2) + ...
                    (Output_Phi_Matrix(Row-1,Col) - Output_Phi_Matrix(Row+1,Col))/2))/2;
            elseif Row == Contour_Bottom_Wall && Col == Contour_Left_Wall % Bottom Left
                q = q + (Epsilon_One * (((Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2) + ...
                    (Output_Phi_Matrix(Row+1,Col) - Output_Phi_Matrix(Row-1,Col))/2))/2;
            elseif Row == Contour_Bottom_Wall && Col == Contour_Right_Wall % Bottom Right
                q = q + (Epsilon_One * (((Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2) + ...
                    (Output_Phi_Matrix(Row+1,Col) - Output_Phi_Matrix(Row-1,Col))/2))/2;
            end
    
            % Top wall without corners
            if Row == Contour_Top_Wall && Col > Contour_Left_Wall && Col < Contour_Right_Wall
                q = q + (Epsilon_Zero * (Output_Phi_Matrix(Row-1,Col) - Output_Phi_Matrix(Row+1,Col))/2);
            end
            
            % Bottom wall without corners
            if Row == Contour_Top_Wall && Col > Contour_Left_Wall && Col < Contour_Right_Wall
                q = q + (Epsilon_One * (Output_Phi_Matrix(Row+1,Col) - Output_Phi_Matrix(Row-1,Col))/2);
            end
    
            % Left wall without corners
            if Col == Contour_Left_Wall && Row > Contour_Top_Wall && Row < Contour_Bottom_Wall
    
                if Row < Dielectric_Boundary_Row
                    q = q + (Epsilon_Zero * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2);
                elseif Row > Dielectric_Boundary_Row
                    q = q + (Epsilon_One * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/2);
                else
                    q = q + ((Epsilon_Zero * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/4) + ...
                        (Epsilon_One * (Output_Phi_Matrix(Row,Col-1) - Output_Phi_Matrix(Row,Col+1))/4));
                end
    
            end
    
            % Right wall without corners
            if Col == Contour_Right_Wall && Row > Contour_Top_Wall && Row < Contour_Bottom_Wall
    
                if Row < Dielectric_Boundary_Row
                    q = q + (Epsilon_Zero * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2);
                elseif Row > Dielectric_Boundary_Row
                    q = q + (Epsilon_One * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/2);
                else
                    q = q + ((Epsilon_Zero * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/4) + ...
                        (Epsilon_One * (Output_Phi_Matrix(Row,Col+1) - Output_Phi_Matrix(Row,Col-1))/4));
                end
    
            end
           
        end
    end

end
