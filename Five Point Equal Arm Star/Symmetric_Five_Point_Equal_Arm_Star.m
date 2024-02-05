%% Reset Code

clear; close all; clc;

Input_Text = {'Structure Width','Structure Height', 'Mesh Size H','Top Electric Potential (V)','Right Electric Potential (V)',...
              'Left Electric Potential (V)', 'Bottom Electric Potential (V)'};

Box_Text_Dimensions = [1 45];

Input = inputdlg(Input_Text, 'Input', [Box_Text_Dimensions; Box_Text_Dimensions; Box_Text_Dimensions; ...
    Box_Text_Dimensions; Box_Text_Dimensions; Box_Text_Dimensions; Box_Text_Dimensions]);

% Set the input provided by the user to individual varibles
Structure_Width = str2num(Input{1}); %#ok<*ST2NM>
Structure_Height = str2num(Input{2});
Mesh_H = str2num(Input{3});

Top_Voltage = str2num(Input{4});
Right_Voltage = str2num(Input{5});
Left_Voltage = str2num(Input{6});
Bottom_Voltage = str2num(Input{7});

Mesh_Row_Length = Structure_Height / Mesh_H - 1;
Mesh_Col_Length = Structure_Width / Mesh_H - 1;

%% Check for Poisson or Laplace application

if Top_Voltage == 0 & Bottom_Voltage == 0 & Right_Voltage == 0 & Left_Voltage == 0
    Column_Index_Multiplier = -2; % -2 based on sample problem
else
    Column_Index_Multiplier = 0;
end

%% Check for symmetry

if mod(Mesh_Row_Length, 2) == 1 || mod(Mesh_Col_Length, 2) == 1
    
    % Symmetry

    if mod(Mesh_Row_Length, 2) == 1 && mod(Mesh_Col_Length, 2) == 1
        
        Mesh_Row_Length = ceil(Mesh_Row_Length / 2);
        Mesh_Col_Length = ceil(Mesh_Col_Length / 2);

        Right_Symmetry = true;
        Bottom_Symmetry = true;

    elseif mod(Mesh_Row_Length, 2) == 1 && mod(Mesh_Col_Length, 2) == 0
        
        Mesh_Row_Length = ceil(Mesh_Row_Length / 2);

        Right_Symmetry = false;
        Bottom_Symmetry = true;

    else

        Mesh_Col_Length = ceil(Mesh_Col_Length / 2);
        
        Right_Symmetry = true;
        Bottom_Symmetry = false;

    end

    Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;
    Column_Vector = ones(Mesh_Row_Length * Mesh_Col_Length, 1) * Column_Index_Multiplier * Mesh_H.^2;

    Mesh = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, 2, Column_Vector, ...
        Top_Voltage, Bottom_Voltage, Right_Voltage, Left_Voltage, Right_Symmetry, Bottom_Symmetry)
    
else
    
    % No Symmetry
    Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;
    Column_Vector = ones(Mesh_Row_Length * Mesh_Col_Length, 1) * Column_Index_Multiplier * Mesh_H.^2;

    Mesh = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, 1, Column_Vector, ...
        Top_Voltage, Bottom_Voltage, Right_Voltage, Left_Voltage, Right_Symmetry, Bottom_Symmetry)

end

%%

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

Formatted_Output_Phi

% Plot the distribution using the color format
imagesc(Formatted_Output_Phi);
colorbar;

%%

function Mesh = Arm_Star(Mesh, Mesh_Row_Length, Mesh_Col_Length, PHI_Multiplier, Column_Vector, ...
    Top_Voltage, Bottom_Voltage, Right_Voltage, Left_Voltage, Right_Symmetry, Bottom_Symmetry) %#ok<*INUSD>

    % PHI_Matrix_Index used to track the phi position within the mesh
    PHI_Matrix_Index = 1;

    for Row = 1:1:Mesh_Row_Length
   
        for Col = 1:1:Mesh_Col_Length
    
            Electric_Potential_Sum = 0;
            
            % Left Check. Same process as right check but for left neighbor.
            if Col - 1 == 0
                Electric_Potential_Sum = Electric_Potential_Sum - Left_Voltage;
            else
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1;
            end
    
            % Right Check
            if Col + 1 <= Mesh_Col_Length
                Mesh(PHI_Matrix_Index,PHI_Matrix_Index + 1) = 1;
            else

                if Right_Symmetry
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1 * PHI_Multiplier;
                else
                    Electric_Potential_Sum = Electric_Potential_Sum - Right_Voltage;
                end

            end
        
            % Top check
            if Row - 1 == 0
                Electric_Potential_Sum = Electric_Potential_Sum - Top_Voltage;
            else
                Mesh(PHI_Matrix_Index - Mesh_Col_Length, PHI_Matrix_Index) = 1;
            end
            
            % Bottom Check
            if Row + 1 <= Mesh_Row_Length
                Mesh(PHI_Matrix_Index + Mesh_Col_Length, PHI_Matrix_Index) = 1;
            else

                if Bottom_Symmetry
                    Mesh(PHI_Matrix_Index, PHI_Matrix_Index - Mesh_Col_Length) = 1 * PHI_Multiplier;
                else
                    Electric_Potential_Sum = Electric_Potential_Sum - Bottom_Voltage;
                end

            end
    
            Column_Vector(PHI_Matrix_Index) = Column_Vector(PHI_Matrix_Index) - Electric_Potential_Sum;
            Electric_Potential_Sum = 0;
    
            PHI_Matrix_Index = PHI_Matrix_Index + 1;
    
        end

    end

end

