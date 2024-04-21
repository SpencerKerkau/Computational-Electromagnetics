clear; close all; clc;

warning('off','all')

Epsilon_Zero = 8.85*10^(-12);
Voltage = 1;

Square_Plate_Dimensions = 1;

Plate_Rows = 10;
Plate_Columns = 10;

Section_Count = Plate_Rows * Plate_Columns;

Width_Of_Sections = Square_Plate_Dimensions / Plate_Columns;
Height_Of_Sections = Square_Plate_Dimensions / Plate_Rows;

C1_X = zeros(Plate_Rows, Plate_Columns);
C1_Y = zeros(Plate_Rows, Plate_Columns);

for Row = 1:1:Plate_Rows
    for Col = 1:1:Plate_Columns
        C1_X(Row,Col) = (Width_Of_Sections / 2) + ((Col - 1) * Width_Of_Sections);
        C1_Y(Row,Col) = (Height_Of_Sections / 2) + ((Row - 1) * Height_Of_Sections);       
    end
end

C1X_Surf = C1_X;
C1Y_Surf = C1_Y;

C1_X = reshape(C1_X.', 1, []);
C1_X = [C1_X, C1_X];

C1_Y = reshape(C1_Y', 1, []);
C1_Y = [C1_Y, C1_Y];

d_over_w = zeros(1,11);
dIdx = 1;

for Distance = 0:0.1:1
    
    B = ones(Section_Count*2, 1) * Voltage;

    for B_Idx = Plate_Rows*Plate_Columns + 1:1: Plate_Rows*Plate_Columns*2
        B(B_Idx, 1) = -Voltage;
    end

    A = zeros(Section_Count*2, Section_Count*2);
    R_Distance_Matrix = zeros(Plate_Rows*Plate_Columns*2, Plate_Rows*Plate_Columns*2);

    for X_Cell = 1:1:Plate_Rows*Plate_Columns*2
        for Y_Cell = 1:1:Plate_Rows*Plate_Columns*2
    
            for Idx = 1:1:Plate_Rows*Plate_Columns*2
                
                if Idx <= Plate_Rows*Plate_Columns
                    R_Distance_Matrix(X_Cell, Idx) = abs(sqrt((C1_X(X_Cell) - C1_X(Idx))^2 + (C1_Y(X_Cell) - C1_Y(Idx))^2));
                else
                    R_Distance_Matrix(X_Cell, Idx) = abs(sqrt((C1_X(X_Cell) - C1_X(Idx))^2 + (C1_Y(X_Cell) - C1_Y(Idx))^2) + Distance^2);
                end
            end
    
        end
    end


    for Row = 1:1:Plate_Rows*Plate_Columns*2
        for Col = 1:1:Plate_Rows*Plate_Columns*2
            
            if Row == Col
                A(Row,Col) = (0.282 * sqrt(pi*(Width_Of_Sections)^2)) / Epsilon_Zero;
            else
                A(Row, Col) = (pi*(Width_Of_Sections)^2) / (4*pi*Epsilon_Zero*R_Distance_Matrix(Row,Col));
            end
    
        end
    end
    
    Rho_s_Vector = A \ B;

    C = sum(Rho_s_Vector) / 2 * Voltage;
    C_Theo = (Epsilon_Zero * (Square_Plate_Dimensions^2)) / Distance;
    
    d_over_w(1, dIdx) = C / C_Theo;

    figure(dIdx)
    Reshaped_Rho = reshape(Rho_s_Vector(101:200,1), Plate_Rows, []);
    surf(C1X_Surf, C1Y_Surf, Reshaped_Rho)

    dIdx = dIdx + 1;

end

figure(dIdx + 1)
plot(0:0.1:1, d_over_w)


