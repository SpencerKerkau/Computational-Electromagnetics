clear; close all; clc;

Number_Of_Sections = 10;
Rod_Length = 1;

X_Distance = 0:(Rod_Length/(Number_Of_Sections-1)):Rod_Length;

B_Matrix = zeros(Number_Of_Sections, 1);
B_Matrix(Number_Of_Sections,1) = 1;

Mesh_Matrix = eye(Number_Of_Sections, Number_Of_Sections);

for Row = 1:1:Number_Of_Sections
    for Col = 1:1:Number_Of_Sections

        if Row == Col
            Mesh_Matrix(Row,Col) = -2 + 1 + X_Distance(Row);
        end
        
        if Row - 1 > 0
            Mesh_Matrix(Row, Row - 1) = 1;
        end

        if Row + 1 <= Number_Of_Sections
            Mesh_Matrix(Row, Row + 1) = 1;
        end

    end
end

A_Matrix = Mesh_Matrix\B_Matrix


imagesc(A_Matrix);
colorbar;
