%% Reset Code

clear; close all; clc;

%% Receive Input for Mesh Size and Voltage

% Take input on the dimensions and electric potential of the structure.

Input_Text = {'Structure Width','Structure Height', 'Mesh Size H','Top Electric Potential (V)','Right Electric Potential (V)',...
              'Left Electric Potential (V)', 'Bottom Electric Potential (V)', 'Center Conductor Potential (V)'};

Box_Text_Dimensions = [1 45];

Input = inputdlg(Input_Text, 'Input', [Box_Text_Dimensions; Box_Text_Dimensions; Box_Text_Dimensions; ...
    Box_Text_Dimensions; Box_Text_Dimensions; Box_Text_Dimensions; Box_Text_Dimensions; Box_Text_Dimensions]);

% Set the input provided by the user to individual varibles
Structure_Width = str2num(Input{1}); %#ok<*ST2NM>
Structure_Height = str2num(Input{2});
Mesh_H = str2num(Input{3});

Top_Voltage = str2num(Input{4});
Right_Voltage = str2num(Input{5});
Left_Voltage = str2num(Input{6});
Bottom_Voltage = str2num(Input{7});
Conductor_Voltage = str2num(Input{8});

Mesh_Row_Length = Structure_Height / Mesh_H - 1;
Mesh_Col_Length = Structure_Width / Mesh_H - 1;


%% Check for Poisson or Laplace application

if Top_Voltage == 0 & Bottom_Voltage == 0 & Right_Voltage == 0 & Left_Voltage == 0 %#ok<*AND2>
    Column_Index_Multiplier = -2;
else
    Column_Index_Multiplier = 0;
end

%% Generate Mesh with centered metal conductor

Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;
Column_Vector = ones(Mesh_Row_Length * Mesh_Col_Length, 1) * Column_Index_Multiplier * Mesh_H.^2;

Mesh_Dimensions = size(Mesh);

fprintf('Mesh Row Length: %d Mesh Col Length: %d \n', Mesh_Col_Length, Mesh_Col_Length);
fprintf('Mesh matrix dimensions: %d %d' , size(Mesh));

Metal_Conductor_Row = floor(Mesh_Row_Length / 2) -1;
Metal_Conductor_Width = (ceil(Mesh_Col_Length/2) - 2: ceil(Mesh_Col_Length/2) + 2);

%% Five Point Arm Star Algorithm

% Track which Phi to encircle throughout the mesh.
PHI_Matrix_Index = 1;

% Iterate through each row in the mesh.
for Row = 1:1:Mesh_Row_Length
    
    for Col = 1:1:Mesh_Col_Length
        
        Electric_Potential_Sum = 0;
        
        % Right Check
        if Col + 1 <= Mesh_Col_Length
            
            if ismember(Row, Metal_Conductor_Row) & ismember(Col+1, Metal_Conductor_Width)
                Electric_Potential_Sum = Electric_Potential_Sum - Conductor_Voltage;
            else
                Mesh(PHI_Matrix_Index,PHI_Matrix_Index+1) = 1;
            end

        else
            Electric_Potential_Sum = Electric_Potential_Sum - Right_Voltage;
        end
        
        % Left Check. Same process as right check but for left neighbor.
        if Col - 1 == 0
            Electric_Potential_Sum = Electric_Potential_Sum - Left_Voltage;
        else

            if ismember(Row, Metal_Conductor_Row) & ismember(Col-1, Metal_Conductor_Width)
                Electric_Potential_Sum = Electric_Potential_Sum - Conductor_Voltage;
            else
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 1;
            end
        end

        % Bottom Check
        if Row + 1 <= Mesh_Row_Length

            if ismember(Row+1, Metal_Conductor_Row) & ismember(Col+1, Metal_Conductor_Width)
                Electric_Potential_Sum = Electric_Potential_Sum - Conductor_Voltage;
            else
                Mesh(PHI_Matrix_Index + Mesh_Col_Length, PHI_Matrix_Index) = 1;
            end
        else
            Electric_Potential_Sum = Electric_Potential_Sum - Bottom_Voltage;
        end

        % Top Check. Same process as bottom check but for top neighbor.
        if Row - 1 == 0
            Electric_Potential_Sum = Electric_Potential_Sum - Top_Voltage;
        else

            if ismember(Row+1, Metal_Conductor_Row) & ismember(Col+1, Metal_Conductor_Width)
                Electric_Potential_Sum = Electric_Potential_Sum - Conductor_Voltage;
            else
                Mesh(PHI_Matrix_Index - Mesh_Col_Length, PHI_Matrix_Index) = 1;
            end
        end
        
        Column_Vector(PHI_Matrix_Index) = Electric_Potential_Sum;
        Electric_Potential_Sum = 0;

        % Move the Phi index to the next phi in the mesh.
        PHI_Matrix_Index = PHI_Matrix_Index + 1;

    end
    
end

%% Matrix Multiplication to Solve for Phi

% Solve for the potential at each phi within the mesh.

Output_Phi = inv(Mesh) * Column_Vector; %#ok<MINV>
Output_Phi = (reshape(Output_Phi,Mesh_Col_Length,Mesh_Row_Length));
imagesc(transpose(Output_Phi));
