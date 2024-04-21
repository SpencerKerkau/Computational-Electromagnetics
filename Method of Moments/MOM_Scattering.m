clear; close all; clc;

Eo = 1;
Dielectric_Constant = 4;

Inner_Radius = 0.25;
Outer_Radius = 0.3;
Radius = (Inner_Radius + Outer_Radius) / 2;

Section_Num = 100;

Phi_Angles = 360 / Section_Num;
Phi = 0:Phi_Angles:360-Phi_Angles;
x = Radius * cosd(Phi);
y = Radius * sind(Phi);

k = 2*pi;

am = sqrt((Outer_Radius^2 - Inner_Radius^2) / Section_Num);
B = Eo * exp(1j * 2 * pi * x)';
A = zeros(Section_Num, Section_Num);

for n = 1:Section_Num
    for m = 1:Section_Num
        if n == m
            A(n, m) = 1 + (Dielectric_Constant- 1) * (1j / 2) * (pi* k * am * besselh(1, 2, k * am) - 2 * 1j);
        else
            A(n, m) = (1j * pi * k * am / 2) * (Dielectric_Constant - 1) * besselj(1, k * am) * besselh(0, 2, k * sqrt((x(n) - x(m))^2 + (y(n) - y(m))^2));
        end
    end
end

X = A \ B;

plot(Phi(1:50), abs(X(1:50)))

hold on
grid on
title('Scattered Electric Field Distribution')
xlabel('Phi Angle (degrees)')
ylabel('Electric Field')
