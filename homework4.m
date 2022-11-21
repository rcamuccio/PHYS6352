% PHYS 6352: Computational Physics
% Homework 4
% Author: Richard Camuccio

clc

% [Problem 1]
%
% In class we calculated numerical derivative of a function y(x) using the
% straightforward formula
%
% z(n) = (y(n+1) - y(n)) / dx
%
% Write an M-function that calculates the numerical derivative of a given
% function y(x) using this formula and assign either to x(n) or to x(n+1)
% and thus obtain left and right derivatives.

clear

N = 100;
a = 0;
b = 10;

m = 1;
N = m*N;
b = m*b;

x = linspace(a, b, N);
y = x.^3 - 15*x - 5;
u = 3*(x.^2) - 15;

[z, xl, xr] = lrderiv(x, y, a, b, N);

figure(1)
plot(x, y, 'k', x, u, 'k:', xl, z, 'b-*', xr, z, 'r-*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, left)', 'dy/dx (numerical, right)', 'Location', 'NorthWest')
hold on

% [Problem 2]
%
% In class we calculated numerical derivative using a more accurate,
% mid-point rule:
%
% z(n+1) = (y(n+2) - y(n)) / (2 * dx)
%
% also known as central derivative. Write an M-function that calculates the
% derivative using this formula and assign z(n) to x(n) directly, but
% remember the issue with the first point. Note: if you don't take care of
% the first point (or z(1)) MATLAB will give it zero value.

clear

N = 100;
a = 0;
b = 10;

m = 1;
N = m*N;
b = m*b;

x = linspace(a, b, N);
y = x.^3 -15*x - 5;
u = 3*(x.^2) - 15;

uz = -15;

[zc, xc] = cderiv(x, y, uz, a, b, N);

figure(2)
plot(x, y, 'k', x, u, 'k:', xc, zc, 'r-*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, central)', 'Location', 'NorthWest')
hold on

% [Problem 3]
%
% Use the two M-functions from Problem 1 and Problem 2 to calculate
% numerically the derivatives of the following functions:
%
% y = (3/2)*x^(2/3) + (18/7)*x^(1/6) + (9/5)*x^(2/3) + (6/13)*x^(13/6)
%
% y = (x + sqrt(x))^(1/3)
%
% y = (9/5)*(x+2)^(-5) - 3*(x+2)^(-4) + 2*(x+2)^(-3) - (1/2)*(x+2)^(-2)
%
% y = sin(x) - sin(x)^3
%
% For each function choose a suitable range for its argument (x). Compare
% the numerical derivatives calculated with the above methods (left, right,
% and central) with the exact analytical derivatives. (Use the formula for
% the error we introduced in class.)
%
% Note: In this problem you need to write an M-script that defines the
% x-variable, then specifies the function (one of the four above) and then
% calculates the numerical derivatives by calling the M-functions you wrote
% in Problem 1 and 2.

clear

N = 100;
a = 0.1;
b = 10;

m = 1;
N = m*N;
b = m*b;

x = linspace(a, b, N);

% Function 1
y1 = (3/2)*x.^(2/3) + (18/7)*x.^(1/6) + (9/5)*x.^(2/3) + (6/13)*x.^(13/6);
u1 = x.^(-1/3) + (3/7)*x.^(-5/6) + (6/5)*x.^(-1/3) + x.^(7/6);
uz1 = 7.72771;

[z1, xl1, xr1] = lrderiv(x, y1, a, b, N);
[zc1, xc1] = cderiv(x, y1, uz1, a, b, N);

ul1 = u1(1:(N-1));
ur1 = u1(2:N);

pl1 = z1 - ul1;
pr1 = z1 - ur1;
pc1 = zc1 - ul1;

error_l1 = sqrt(sum(pl1.^2) / N);
error_r1 = sqrt(sum(pr1.^2) / N);
error_c1 = sqrt(sum(pc1.^2) / (N-1));

% Function 2
y2 = (x + x.^(1/2)).^(1/3);
u2 = (1/3) * ((x + x.^(1/2)).^(-2/3)) .* (1 + ((1/2) * x.^(-1/2)));
uz2 = 1.54337;

[z2, xl2, xr2] = lrderiv(x, y2, a, b, N);
[zc2, xc2] = cderiv(x, y2, uz2, a, b, N);

ul2 = u2(1:(N-1));
ur2 = u2(2:N);

pl2 = z2 - ul2;
pr2 = z2 - ur2;
pc2 = zc2 - ul2;

error_l2 = sqrt(sum(pl2.^2) / N);
error_r2 = sqrt(sum(pr2.^2) / N);
error_c2 = sqrt(sum(pc2.^2) / (N-1));

% Function 3
y3 = (9/5)*((x+2).^(-5)) - 3*((x+2).^(-4)) + 2*((x+2).^(-3)) - (1/2)*((x+2).^(-2));
u3 = -9*((x+2).^(-6)) + 12*((x+2).^(-5)) - 6*((x+2).^(-4)) + ((x+2).^(-3));
uz3 = -0.01164;

[z3, xl3, xr3] = lrderiv(x, y3, a, b, N);
[zc3, xc3] = cderiv(x, y3, uz3, a, b, N);

ul3 = u3(1:(N-1));
ur3 = u3(2:N);

pl3 = z3 - ul3;
pr3 = z3 - ur3;
pc3 = zc3 - ul3;

error_l3 = sqrt(sum(pl3.^2) / N);
error_r3 = sqrt(sum(pr3.^2) / N);
error_c3 = sqrt(sum(pc3.^2) / (N-1));

% Function 4
y4 = sin(x) - (sin(x)).^3;
u4 = cos(x) - 3 * ((sin(x)).^2) .* cos(x);
uz4 = 1;

[z4, xl4, xr4] = lrderiv(x, y4, a, b, N);
[zc4, xc4] = cderiv(x, y4, uz4, a, b, N);

ul4 = u4(1:(N-1));
ur4 = u4(2:N);

pl4 = z4 - ul4;
pr4 = z4 - ur4;
pc4 = zc4 - ul4;

error_l4 = sqrt(sum(pl4.^2) / N);
error_r4 = sqrt(sum(pr4.^2) / N);
error_c4 = sqrt(sum(pc4.^2) / (N-1));

% Display
disp('Function 1:')
disp(['  Left derivative error: ', num2str(error_l1)])
disp(['  Right derivative error: ', num2str(error_r1)])
disp(['  Central derivative error: ', num2str(error_c1)])
disp(' ')
disp('Function 2:')
disp(['  Left derivative error: ', num2str(error_l2)])
disp(['  Right derivative error: ', num2str(error_r2)])
disp(['  Central derivative error: ', num2str(error_c2)])
disp(' ')
disp('Function 3:')
disp(['  Left derivative error: ', num2str(error_l3)])
disp(['  Right derivative error: ', num2str(error_r3)])
disp(['  Central derivative error: ', num2str(error_c3)])
disp(' ')
disp('Function 4:')
disp(['  Left derivative error: ', num2str(error_l4)])
disp(['  Right derivative error: ', num2str(error_r4)])
disp(['  Central derivative error: ', num2str(error_c4)])
disp(' ')

% Plots

figure(3)
plot(x, y1, 'k', x, u1, 'k:', xl1, z1, 'b*', xr1, z1, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, left)', 'dy/dx (numerical, right)', 'Location', 'NorthWest')
hold on

figure(4)
plot(x, y1, 'k', x, u1, 'k:', xc1, zc1, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, central)', 'Location', 'NorthWest')
hold on

figure(5)
plot(x, y2, 'k', x, u2, 'k:', xl2, z2, 'b*', xr2, z2, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, left)', 'dy/dx (numerical, right)', 'Location', 'NorthWest')
hold on

figure(6)
plot(x, y2, 'k', x, u2, 'k:', xc2, zc2, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, central)', 'Location', 'NorthWest')
hold on

figure(7)
plot(x, y3, 'k', x, u3, 'k:', xl3, z3, 'b*', xr3, z3, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, left)', 'dy/dx (numerical, right)', 'Location', 'SouthEast')
hold on

figure(8)
plot(x, y3, 'k', x, u3, 'k:', xc3, zc3, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, central)', 'Location', 'SouthEast')
hold on

figure(9)
plot(x, y4, 'k', x, u4, 'k:', xl4, z4, 'b*', xr4, z4, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, left)', 'dy/dx (numerical, right)', 'Location', 'SouthEast')
hold on

figure(10)
plot(x, y4, 'k', x, u4, 'k:', xc4, zc4, 'r*')
grid
legend('y(x)', 'dy/dx (exact)', 'dy/dx (numerical, central)', 'Location', 'SouthEast')
hold on

% [Functions]

function [z, xl, xr] = lrderiv(x, y, a, b, N)
    
    dx = (b - a) / (N - 1);

    for n = 1:(N-1)

        dy = y(n+1) - y(n);
        z(n) = dy / dx;

    end

    xl = x(1:(N-1));
    xr = x(2:N);

end

function [zc, xc] = cderiv(x, y, uz, a, b, N)

    dx = (b - a) / (N - 1);

    for n = 1:(N-2)

        dy = y(n+2) - y(n);
        zc(n+1) = dy / (2 * dx);

    end

    zc(1) = uz;
    xc = x(1:(N-1));

end