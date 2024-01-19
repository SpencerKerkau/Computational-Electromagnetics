%% Reset Code

clear; close all; clc;

%% Receive Input for Mesh Size and Voltage

% Take input on the dimensions and electric potential of the structure.
Input = inputdlg({'Structure Width','Structure Height', 'Mesh Size H','Top Electric Potential (V)','Right Electric Potential (V)',...
              'Left Electric Potential (V)', 'Bottom Electric Potential (V)'}, 'Input',[1 45; 1 45; 1 45; 1 45; 1 45; 1 45; 1 45]);

% Set the input provided by the user to individual varibles
Structure_Width = str2num(Input{1}); %#ok<*ST2NM>
Structure_Height = str2num(Input{2});
Mesh_H = str2num(Input{3});

Top_Voltage = str2num(Input{4});
Right_Voltage = str2num(Input{5});
Left_Voltage = str2num(Input{6});
Bottom_Voltage = str2num(Input{7});

Total_Number_of_Square_Mesh = (Structure_Width * Structure_Height) / (Mesh_H.^2);
Dimension_Factors = factor(Total_Number_of_Square_Mesh);

%% Create Proper Mesh Sizes

if size(Dimension_Factors,2) == 2

    Lower_Factors = Dimension_Factors(1);
    Higher_Factors = Dimension_Factors(2);

else if size(Dimension_Factors,2) == 3

    Lower_Factors = Dimension_Factors(1) * Dimension_Factors(2);
    Higher_Factors = Dimension_Factors(3);

else

    Lower_Factors = 1;

    for i = 1:1:size(Dimension_Factors,2) - 2  
        Lower_Factors = Lower_Factors * Dimension_Factors(i);
    end
    
    Higher_Factors = Dimension_Factors(end) * Dimension_Factors(end-1);

end

end

Lower_Factors = Lower_Factors - 1;
Higher_Factors = Higher_Factors - 1;

if Structure_Height >= Structure_Width

    Mesh_Row_Length = Lower_Factors;
    Mesh_Col_Length = Higher_Factors;

else

    Mesh_Row_Length = Higher_Factors;
    Mesh_Col_Length = Lower_Factors;

end

%% Check for Poisson or Laplace application

if Top_Voltage == 0 & Bottom_Voltage == 0 & Right_Voltage == 0 & Left_Voltage == 0
    Column_Index_Multiplier = -2;
else
    Column_Index_Multiplier = 0;
end

%% Check for symmetry

PHI_Matrix_Index = 1;

if mod(Mesh_Row_Length - 1, 2) == 0 & mod(Mesh_Col_Length - 1, 2) == 0 & Top_Voltage == Bottom_Voltage & Right_Voltage == Left_Voltage%#ok<*AND2>
        
    Symmetry = true;

    Mesh_Row_Length = ceil(Mesh_Row_Length / 2);
    Mesh_Col_Length = ceil(Mesh_Col_Length / 2);
    
    Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;
    Column_Vector = ones(Mesh_Row_Length * Mesh_Col_Length, 1) * Column_Index_Multiplier * Mesh_H.^2;

else

    Symmetry = false;
    Mesh = eye(Mesh_Row_Length * Mesh_Col_Length) * -4;
    Column_Vector = ones(Mesh_Row_Length * Mesh_Col_Length, 1) * Column_Index_Multiplier * Mesh_H.^2;
end

%% Non Symmetric Arm Star

if Symmetry == false
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
            
            % Set the sum of the potential voltages to the Phi index in the column
            % vector. Then reset the potential sum.
            Column_Vector(PHI_Matrix_Index) = Column_Vector(PHI_Matrix_Index) + Electric_Potential_Sum;
            Electric_Potential_Sum = 0;
    
            % Move the Phi index to the next phi in the mesh.
            PHI_Matrix_Index = PHI_Matrix_Index + 1;
        end
        
    end
end
%% Symmetric Arm Star

if Symmetry == true

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
            if Col + 1 > Mesh_Col_Length
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index - 1) = 2;
            else
                Mesh(PHI_Matrix_Index,PHI_Matrix_Index + 1) = 1;
            end
        
            % Top check
            if Row - 1 == 0
                Electric_Potential_Sum = Electric_Potential_Sum - Top_Voltage;
            else
                Mesh(PHI_Matrix_Index - Mesh_Col_Length, PHI_Matrix_Index) = 1;
            end
            
            % Bottom Check
            if Row + 1 > Mesh_Row_Length
                Mesh(PHI_Matrix_Index, PHI_Matrix_Index - Mesh_Col_Length) = 2;
            else
                Mesh(PHI_Matrix_Index + Mesh_Col_Length, PHI_Matrix_Index) = 1;
            end

            Column_Vector(PHI_Matrix_Index) = Column_Vector(PHI_Matrix_Index) - Electric_Potential_Sum;
            Electric_Potential_Sum = 0;
    
            PHI_Matrix_Index = PHI_Matrix_Index + 1;

        end
    
    end
end

%% Matrix Multiplication to Solve for Phi

% Solve for the potential at each phi within the mesh.

Output_Phi = inv(Mesh) * Column_Vector; %#ok<MINV>

figure
mesh(Mesh);

disp(transpose(reshape(Output_Phi,Mesh_Col_Length,Mesh_Row_Length)));


