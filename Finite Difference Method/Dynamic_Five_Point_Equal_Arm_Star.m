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

% Create identity matrix to represent mesh. Each phi location within matrix
% set to -4 per equal arm algorithm.
Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;

% Create column vector to store sum of electric potential at each phi.
Column_Vector = zeros(Mesh_Row_Length * Mesh_Col_Length, 1);

%% Five Point Arm Star Algorithm

% Track which Phi to encircle throughout the mesh.
PHI_Matrix_Index = 1;

%                                   Top Voltage
%              | Phi[1] = (1,1) | Phi[2] = (1,2) | Phi[3] = (1,3) |
% Left Voltage | Phi[4] = (2,1) | Phi[5] = (2,2) | Phi[6] = (2,3) |   Right Voltage
%              | Phi[7] = (3,1) | Phi[8] = (3,2) | Phi[9] = (3,3) |
%                                 Bottom Voltage

% Iterate through each row in the mesh.
for Row = 1:1:Mesh_Row_Length
    
    % Iterate through each column within each row. This allows access to
    % each individual phi within the mesh.
    for Col = 1:1:Mesh_Col_Length
        
        % Column vector sum required when Phi is adjacent to the side of
        % the structure that has an electric potential.
        Electric_Potential_Sum = 0;
        
        % Right Check
        if Col + 1 <= Mesh_Col_Length

            % If right neighbor is another Phi point, set the right
            % neighbor within the mesh matrix to 1 to indicate existing
            % neighbor.

            Mesh(PHI_Matrix_Index,PHI_Matrix_Index+1) = 1;
        else

            % If right neighbor is the side of the structure, update the
            % column vector sum. Subtrack out electric potential as the
            % potential will be summed and subtracted out during math
            % calculation anyway.

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

            % If bottom neighbor is the bottom of the structure, update the
            % column vector sum. Subtrack out electric potential as the
            % potential will be summed and subtracted out during math
            % calculation anyway.
            
            % Instead of adding or subtracting 1 when checking left and
            % right, add and subtract the length of the column to move to
            % the next row. For example, looking at the visual on line 30,
            % there are 3 columns in each row. To go from Phi[2] to Phi[5],
            % you must do Phi[2 + # of columns].

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
        
        % Set the sum of the potential voltages to the Phi index in the column
        % vector. Then reset the potential sum.
        Column_Vector(PHI_Matrix_Index) = Electric_Potential_Sum;
        Electric_Potential_Sum = 0;

        % Move the Phi index to the next phi in the mesh.
        PHI_Matrix_Index = PHI_Matrix_Index + 1;
    end
    
end


%% Matrix Multiplication to Solve for Phi

% Solve for the potential at each phi within the mesh.

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

% Plot the distribution using the color format
imagesc(Formatted_Output_Phi);
colorbar;

