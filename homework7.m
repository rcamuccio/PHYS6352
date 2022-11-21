% PHYS 6352: Computational Physics
% Homework 7
% Author: Richard Camuccio

% This assignment contains the same integrals as in HW 6 except that now
% each integral is a function of x. (The variable x takes values in some
% interval [a, b].) So, you can simply take your M-scripts from HW 6 and
% modify them accordingly.

% Calculate numerically the following integrals. For each integral use all
% three methods: rule of rectangles, rule of trapezoids, and Simpson rule.
% Choose a suitable interval [a, b] for the calculations and also choose a
% suitable number of points N within each interval. (Since N does not
% matter much, you can simply set N = 100 or N = 101 for Simpson to work.)
%
% Compare the approximate numerical integral I(x) with the exact analytical
% answer J(x) by plotting them on the same figure and also by using the
% formula for the error we introduced in class.

clc
clear

N = 101;

% --- Integral 1 ---
a = -10;
b = 10;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = (x + 3) ./ (x.^2 - 2*x + 5);
wrapper(N, a, h, x, y, 1);

% --- Integral 2 ---
a = 0;
b = 30;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = (5*x + 3) ./ sqrt(x.^2 + 4*x + 10);
wrapper(N, a, h, x, y, 2);

% --- Integral 3 ---
a = -10;
b = 10;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = atan(x);
wrapper(N, a, h, x, y, 3);

% --- Integral 4 ---
a = 7;
b = 10;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = (x.^2) .* exp(x);
wrapper(N, a, h, x, y, 4);

% --- Integral 5 ---
a = -10;
b = 10;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = (x.^2 + 7*x - 5) .* cos(2*x);
wrapper(N, a, h, x, y, 5);

% --- Integral 6 ---
a = -2;
b = 2;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = sqrt(4 - x.^2);
wrapper(N, a, h, x, y, 6);

% --- Integral 7 ---
a = -10;
b = 10;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = (x.^2) ./ ((x.^2 + 2).^2);
wrapper(N, a, h, x, y, 7);

% --- Integral 8 ---
a = -5;
b = 0;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = x ./ ((x.^2 + 1) .* (x - 1));
wrapper(N, a, h, x, y, 8);

% --- Integral 9 ---
a = -2*pi;
b = 2*pi;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = 1 ./ (2 - sin(x).^2);
wrapper(N, a, h, x, y, 9);

% --- Integral 10 ---
a = pi/3;
b = pi/2;
h = (b - a) / (N - 1);
x = linspace(a, b, N);
y = (cos(x).^3) ./ (sin(x).^4);
wrapper(N, a, h, x, y, 10);

% --- Functions ---

function wrapper(N, a, h, x, y, n)

    Ir = int_rect(y, h, N);
    It = int_trap(y, h, N);
    [Js, u] = int_simp(x, y, h, N);
    Ia = int_exact(x, a, n);

    figure(n)
    subplot(1, 2, 1)
    plot(x, y, 'k')
    grid
    axis padded

    subplot(1, 2, 2)
    plot(x, Ir, 'r', x, It, 'g', u, Js, 'b', x, Ia, 'k')
    grid
    axis padded
    legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

    Rr = Ir - Ia;
    Rt = It - Ia;
    error_Rr = sqrt(sum(Rr.^2) / N);
    error_Rt = sqrt(sum(Rt.^2) / N);

    Ja = Ia(1:2:N);
    Rs = Js - Ja;
    error_Rs = sqrt(sum(Rs.^2) / N);

    disp(['[Integral ', num2str(n), ']'])
    disp(['  Error (rect): ', num2str(error_Rr)])
    disp(['  Error (trap): ', num2str(error_Rt)])
    disp(['  Error (Simp): ', num2str(error_Rs)])

end

function I = int_rect(y, h, N)

    I(1) = 0;
    
    for n = 1:1:(N-1)
        I(n+1) = I(n) + y(n) * h;
    end

end

function I = int_trap(y, h, N)

    I(1) = 0;

    for n = 1:1:(N-1)
        I(n+1) = I(n) + ((y(n) + y(n+1)) / 2) * h;
    end

end

function [J, u] = int_simp(x, y, h, N)

    J(1) = 0;

    for n = 2:2:(N-1)
        J(n/2+1) = J(n/2) + (y(n-1) + 4*y(n) + y(n+1)) * h/3;
   
    end

    u = x(1:2:N);

end

function I = int_exact(x, a, int_num)

    if (int_num == 1)
        
        alpha = 1 + 2i;
        beta = 1 - 2i;
        chi = (3 + alpha) / (alpha - beta);
        xi = (3 + beta) / (alpha - beta);

        I = chi * log((x-alpha) / (a-alpha)) - xi * log((x-beta) / (a-beta));

    elseif (int_num == 2)

        I = (-5*atanh((x + 2) / sqrt(x.^2 + 4*x + 10))) + (5*sqrt(x.^2 + 4*x + 10)) + (5*atanh((a + 2) / sqrt(a^2 + 4*a + 10))) - (5*sqrt(a^2 + 4*a + 10));

    elseif (int_num == 3)

        I = x.*atan(x) - a*atan(a) - (1/2)*log((1+x.^2) / (1+a^2));

    elseif (int_num == 4)

        I = exp(x).*(x.^2 - 2*x + 2) - exp(a)*(a^2 - 2*a + 2);

    elseif (int_num == 5)

        I = (x/2).*cos(2*x) + ((x.^2)/2).*sin(2*x) - (1/4)*sin(2*x) - (a/2)*cos(2*a) - ((a^2)/2)*sin(2*a) + (1/4)*sin(2*a) + (7/4)*cos(2*x) + (7*x/2).*sin(2*x) - (7/4)*cos(2*a) - (7*a/2)*sin(2*a) - (5/2)*sin(2*x) + (5/2)*sin(2*a);

    elseif (int_num == 6)

        I = (x/2).*sqrt(4-x.^2) + 2*asin(x/2) - (a/2)*sqrt(4-a^2) - 2*asin(a/2);

    elseif (int_num == 7)

        I = atan(x/sqrt(2))/(2*sqrt(2)) - x./(4 + 2*x.^2) - atan(a/sqrt(2))/(2*sqrt(2)) + a/(4 + 2*a^2);

    elseif (int_num == 8)

        I = (1/4) * (-log(x.^2 + 1) + 2*log(1 - x) + 2*atan(x) + log(a^2 + 1) - 2*log(1 - a) - 2*atan(a));

    elseif (int_num == 9)

        I = (x/sqrt(2)) - (a/sqrt(2));

    elseif (int_num == 10)

        I = (1 ./ sin(x)) - (1 / sin(a)) + (1 / (3*(sin(a)^3))) - (1 ./ (3*(sin(x).^3)));

    end

end