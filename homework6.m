% PHYS 6352: Computational Physics
% Homework 6
% Author: Richard Camuccio

% In this assignment you need to submit your M-scripts and analytical
% calculations for each integral. Analytical calculations need to be done
% with pen and paper.

% Calculate numerically the following integrals

clc
clear

N = 100;

% --- Integral 1 ---------------------------------------------------------
a1 = -10;
b1 = 10;
h1 = (b1 - a1) / (N - 1);
x1 = linspace(a1, b1, N);
y1 = (x1 + 3) ./ (x1.^2 - 2*x1 + 5);

[I_nr1, I_r1] = int_rect(x1, y1, h1, N);
[I_nt1, I_t1] = int_trap(x1, y1, h1, N);
[J_ns1, I_s1, u1] = int_simp(x1, y1, h1, N);
[I_na1, I_a1] = int_exact(x1, a1, b1, 1);

figure(1)
subplot(1, 2, 1)
plot(x1, y1, 'k')
grid
axis padded

subplot(1, 2, 2)
plot(x1, I_nr1, 'r', x1, I_nt1, 'g', u1, J_ns1, 'b', x1, I_na1, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r1 = I_r1 - I_a1;
diff_t1 = I_t1 - I_a1;
diff_s1 = I_s1 - I_a1;

R_r1 = I_nr1 - I_na1;
R_t1 = I_nt1 - I_na1;

error_R_r1 = sqrt(sum(R_r1.^2)/N);
error_R_t1 = sqrt(sum(R_t1.^2)/N);

J_na1 = I_na1(1:2:N);
R_s1 = J_ns1 - J_na1;
error_R_s1 = sqrt(sum(R_s1.^2)/N);

disp('[Integral 1]')
disp(['  Analytical: ', num2str(I_a1)])
disp(['  Rectangles: ', num2str(I_r1), ', Diff: ', num2str(diff_r1), ', Error: ', num2str(error_R_r1)])
disp(['  Trapezoids: ', num2str(I_t1), ', Diff: ', num2str(diff_t1), ', Error: ', num2str(error_R_t1)])
disp(['  Simpsons:   ', num2str(I_s1), ', Diff: ', num2str(diff_s1), ', Error: ', num2str(error_R_s1)])

% --- Integral 2 ---------------------------------------------------------
a2 = 0;
b2 = 30;
h2 = (b2 - a2) / (N - 1);
x2 = linspace(a2, b2, N);
y2 = (5*x2 + 3) ./ sqrt(x2.^2 + 4*x2 + 10);

[I_nr2, I_r2] = int_rect(x2, y2, h2, N);
[I_nt2, I_t2] = int_trap(x2, y2, h2, N);
[J_ns2, I_s2, u2] = int_simp(x2, y2, h2, N);
[I_na2, I_a2] = int_exact(x2, a2, b2, 2);

figure(2)
subplot(1, 2, 1)
plot(x2, y2, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x2, I_nr2, 'r', x2, I_nt2, 'g', u2, J_ns2, 'b', x2, I_na2, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r2 = I_r2 - I_a2;
diff_t2 = I_t2 - I_a2;
diff_s2 = I_s2 - I_a2;

R_r2 = I_nr2 - I_na2;
R_t2 = I_nt2 - I_na2;

error_R_r2 = sqrt(sum(R_r2.^2)/N);
error_R_t2 = sqrt(sum(R_t2.^2)/N);

J_na2 = I_na2(1:2:N);
R_s2 = J_ns2 - J_na2;
error_R_s2 = sqrt(sum(R_s2.^2)/N);

disp('[Integral 2]')
disp(['  Analytical: ', num2str(I_a2)])
disp(['  Rectangles: ', num2str(I_r2), ', Diff: ', num2str(diff_r2), ', Error: ', num2str(error_R_r2)])
disp(['  Trapezoids: ', num2str(I_t2), ', Diff: ', num2str(diff_t2), ', Error: ', num2str(error_R_t2)])
disp(['  Simpsons:   ', num2str(I_s2), ', Diff: ', num2str(diff_s2), ', Error: ', num2str(error_R_s2)])

% --- Integral 3 ---------------------------------------------------------
a3 = -10;
b3 = 10;
h3 = (b3 - a3) / (N - 1);
x3 = linspace(a3, b3, N);
y3 = atan(x3);

[I_nr3, I_r3] = int_rect(x3, y3, h3, N);
[I_nt3, I_t3] = int_trap(x3, y3, h3, N);
[J_ns3, I_s3, u3] = int_simp(x3, y3, h3, N);
[I_na3, I_a3] = int_exact(x3, a3, b3, 3);

figure(3)
subplot(1, 2, 1)
plot(x3, y3, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x3, I_nr3, 'r', x3, I_nt3, 'g', u3, J_ns3, 'b', x3, I_na3, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r3 = I_r3 - I_a3;
diff_t3 = I_t3 - I_a3;
diff_s3 = I_s3 - I_a3;

R_r3 = I_nr3 - I_na3;
R_t3 = I_nt3 - I_na3;

error_R_r3 = sqrt(sum(R_r3.^2)/N);
error_R_t3 = sqrt(sum(R_t3.^2)/N);

J_na3 = I_na3(1:2:N);
R_s3 = J_ns3 - J_na3;
error_R_s3 = sqrt(sum(R_s3.^2)/N);

disp('[Integral 3]')
disp(['  Analytical: ', num2str(I_a3)])
disp(['  Rectangles: ', num2str(I_r3), ', Diff: ', num2str(diff_r3), ', Error: ', num2str(error_R_r3)])
disp(['  Trapezoids: ', num2str(I_t3), ', Diff: ', num2str(diff_t3), ', Error: ', num2str(error_R_t3)])
disp(['  Simpsons:   ', num2str(I_s3), ', Diff: ', num2str(diff_s3), ', Error: ', num2str(error_R_s3)])

% --- Integral 4 ---------------------------------------------------------
a4 = 7;
b4 = 10;
h4 = (b4 - a4) / (N - 1);
x4 = linspace(a4, b4, N);
y4 = (x4.^2) .* exp(x4);

[I_nr4, I_r4] = int_rect(x4, y4, h4, N);
[I_nt4, I_t4] = int_trap(x4, y4, h4, N);
[J_ns4, I_s4, u4] = int_simp(x4, y4, h4, N);
[I_na4, I_a4] = int_exact(x4, a4, b4, 4);

figure(4)
subplot(1, 2, 1)
plot(x4, y4, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x4, I_nr4, 'r', x4, I_nt4, 'g', u4, J_ns4, 'b', x4, I_na4, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r4 = I_r4 - I_a4;
diff_t4 = I_t4 - I_a4;
diff_s4 = I_s4 - I_a4;

R_r4 = I_nr4 - I_na4;
R_t4 = I_nt4 - I_na4;

error_R_r4 = sqrt(sum(R_r4.^2)/N);
error_R_t4 = sqrt(sum(R_t4.^2)/N);

J_na4 = I_na4(1:2:N);
R_s4 = J_ns4 - J_na4;
error_R_s4 = sqrt(sum(R_s4.^2)/N);

disp('[Integral 4]')
disp(['  Analytical: ', num2str(I_a4)])
disp(['  Rectangles: ', num2str(I_r4), ', Diff: ', num2str(diff_r4), ', Error: ', num2str(error_R_r4)])
disp(['  Trapezoids: ', num2str(I_t4), ', Diff: ', num2str(diff_t4), ', Error: ', num2str(error_R_t4)])
disp(['  Simpsons:   ', num2str(I_s4), ', Diff: ', num2str(diff_s4), ', Error: ', num2str(error_R_s4)])

% --- Integral 5 ---------------------------------------------------------
a5 = -10;
b5 = 10;
h5 = (b5 - a5) / (N - 1);
x5 = linspace(a5, b5, N);
y5 = (x5.^2 + 7*x5 - 5) .* cos(2*x5);

[I_nr5, I_r5] = int_rect(x5, y5, h5, N);
[I_nt5, I_t5] = int_trap(x5, y5, h5, N);
[J_ns5, I_s5, u5] = int_simp(x5, y5, h5, N);
[I_na5, I_a5] = int_exact(x5, a5, b5, 5);

figure(5)
subplot(1, 2, 1)
plot(x5, y5, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x5, I_nr5, 'r', x5, I_nt5, 'g', u5, J_ns5, 'b', x5, I_na5, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r5 = I_r5 - I_a5;
diff_t5 = I_t5 - I_a5;
diff_s5 = I_s5 - I_a5;

R_r5 = I_nr5 - I_na5;
R_t5 = I_nt5 - I_na5;

error_R_r5 = sqrt(sum(R_r5.^2)/N);
error_R_t5 = sqrt(sum(R_t5.^2)/N);

J_na5 = I_na5(1:2:N);
R_s5 = J_ns5 - J_na5;
error_R_s5 = sqrt(sum(R_s5.^2)/N);

disp('[Integral 5]')
disp(['  Analytical: ', num2str(I_a5)])
disp(['  Rectangles: ', num2str(I_r5), ', Diff: ', num2str(diff_r5), ', Error: ', num2str(error_R_r5)])
disp(['  Trapezoids: ', num2str(I_t5), ', Diff: ', num2str(diff_t5), ', Error: ', num2str(error_R_t5)])
disp(['  Simpsons:   ', num2str(I_s5), ', Diff: ', num2str(diff_s5), ', Error: ', num2str(error_R_s5)])

% --- Integral 6 ---------------------------------------------------------
a6 = -2;
b6 = 2;
h6 = (b6 - a6) / (N - 1);
x6 = linspace(a6, b6, N);
y6 = sqrt(4 - x6.^2);

[I_nr6, I_r6] = int_rect(x6, y6, h6, N);
[I_nt6, I_t6] = int_trap(x6, y6, h6, N);
[J_ns6, I_s6, u6] = int_simp(x6, y6, h6, N);
[I_na6, I_a6] = int_exact(x6, a6, b6, 6);

figure(6)
subplot(1, 2, 1)
plot(x6, y6, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x6, I_nr6, 'r', x6, I_nt6, 'g', u6, J_ns6, 'b', x6, I_na6, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r6 = I_r6 - I_a6;
diff_t6 = I_t6 - I_a6;
diff_s6 = I_s6 - I_a6;

R_r6 = I_nr6 - I_na6;
R_t6 = I_nt6 - I_na6;

error_R_r6 = sqrt(sum(R_r6.^2)/N);
error_R_t6 = sqrt(sum(R_t6.^2)/N);

J_na6 = I_na6(1:2:N);
R_s6 = J_ns6 - J_na6;
error_R_s6 = sqrt(sum(R_s6.^2)/N);

disp('[Integral 6]')
disp(['  Analytical: ', num2str(I_a6)])
disp(['  Rectangles: ', num2str(I_r6), ', Diff: ', num2str(diff_r6), ', Error: ', num2str(error_R_r6)])
disp(['  Trapezoids: ', num2str(I_t6), ', Diff: ', num2str(diff_t6), ', Error: ', num2str(error_R_t6)])
disp(['  Simpsons:   ', num2str(I_s6), ', Diff: ', num2str(diff_s6), ', Error: ', num2str(error_R_s6)])

% --- Integral 7 ---------------------------------------------------------
a7 = -10;
b7 = 10;
h7 = (b7 - a7) / (N - 1);
x7 = linspace(a7, b7, N);
y7 = (x7.^2) ./ ((x7.^2 + 2).^2);

[I_nr7, I_r7] = int_rect(x7, y7, h7, N);
[I_nt7, I_t7] = int_trap(x7, y7, h7, N);
[J_ns7, I_s7, u7] = int_simp(x7, y7, h7, N);
[I_na7, I_a7] = int_exact(x7, a7, b7, 7);

figure(7)
subplot(1, 2, 1)
plot(x7, y7, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x7, I_nr7, 'r', x7, I_nt7, 'g', u7, J_ns7, 'b', x7, I_na7, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r7 = I_r7 - I_a7;
diff_t7 = I_t7 - I_a7;
diff_s7 = I_s7 - I_a7;

R_r7 = I_nr7 - I_na7;
R_t7 = I_nt7 - I_na7;

error_R_r7 = sqrt(sum(R_r7.^2)/N);
error_R_t7 = sqrt(sum(R_t7.^2)/N);

J_na7 = I_na7(1:2:N);
R_s7 = J_ns7 - J_na7;
error_R_s7 = sqrt(sum(R_s7.^2)/N);

disp('[Integral 7]')
disp(['  Analytical: ', num2str(I_a7)])
disp(['  Rectangles: ', num2str(I_r7), ', Diff: ', num2str(diff_r7), ', Error: ', num2str(error_R_r7)])
disp(['  Trapezoids: ', num2str(I_t7), ', Diff: ', num2str(diff_t7), ', Error: ', num2str(error_R_t7)])
disp(['  Simpsons:   ', num2str(I_s7), ', Diff: ', num2str(diff_s7), ', Error: ', num2str(error_R_s7)])

% --- Integral 8 ---------------------------------------------------------
a8 = -5;
b8 = 0;
h8 = (b8 - a8) / (N - 1);
x8 = linspace(a8, b8, N);
y8 = x8 ./ ((x8.^2 + 1) .* (x8 - 1));

[I_nr8, I_r8] = int_rect(x8, y8, h8, N);
[I_nt8, I_t8] = int_trap(x8, y8, h8, N);
[J_ns8, I_s8, u8] = int_simp(x8, y8, h8, N);
[I_na8, I_a8] = int_exact(x8, a8, b8, 8);

figure(8)
subplot(1, 2, 1)
plot(x8, y8, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x8, I_nr8, 'r', x8, I_nt8, 'g', u8, J_ns8, 'b', x8, I_na8, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r8 = I_r8 - I_a8;
diff_t8 = I_t8 - I_a8;
diff_s8 = I_s8 - I_a8;

R_r8 = I_nr8 - I_na8;
R_t8 = I_nt8 - I_na8;

error_R_r8 = sqrt(sum(R_r8.^2)/N);
error_R_t8 = sqrt(sum(R_t8.^2)/N);

J_na8 = I_na8(1:2:N);
R_s8 = J_ns8 - J_na8;
error_R_s8 = sqrt(sum(R_s8.^2)/N);

disp('[Integral 8]')
disp(['  Analytical: ', num2str(I_a8)])
disp(['  Rectangles: ', num2str(I_r8), ', Diff: ', num2str(diff_r8), ', Error: ', num2str(error_R_r8)])
disp(['  Trapezoids: ', num2str(I_t8), ', Diff: ', num2str(diff_t8), ', Error: ', num2str(error_R_t8)])
disp(['  Simpsons:   ', num2str(I_s8), ', Diff: ', num2str(diff_s8), ', Error: ', num2str(error_R_s8)])

% --- Integral 9 ---------------------------------------------------------
a9 = -2*pi;
b9 = 2*pi;
h9 = (b9 - a9) / (N - 1);
x9 = linspace(a9, b9, N);
y9 = 1 ./ (2 - sin(x9).^2);

[I_nr9, I_r9] = int_rect(x9, y9, h9, N);
[I_nt9, I_t9] = int_trap(x9, y9, h9, N);
[J_ns9, I_s9, u9] = int_simp(x9, y9, h9, N);
[I_na9, I_a9] = int_exact(x9, a9, b9, 9);

figure(9)
subplot(1, 2, 1)
plot(x9, y9, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x9, I_nr9, 'r', x9, I_nt9, 'g', u9, J_ns9, 'b', x9, I_na9, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r9 = I_r9 - I_a9;
diff_t9 = I_t9 - I_a9;
diff_s9 = I_s9 - I_a9;

R_r9 = I_nr9 - I_na9;
R_t9 = I_nt9 - I_na9;

error_R_r9 = sqrt(sum(R_r9.^2)/N);
error_R_t9 = sqrt(sum(R_t9.^2)/N);

J_na9 = I_na9(1:2:N);
R_s9 = J_ns9 - J_na9;
error_R_s9 = sqrt(sum(R_s9.^2)/N);

disp('[Integral 9]')
disp(['  Analytical: ', num2str(I_a9)])
disp(['  Rectangles: ', num2str(I_r9), ', Diff: ', num2str(diff_r9), ', Error: ', num2str(error_R_r9)])
disp(['  Trapezoids: ', num2str(I_t9), ', Diff: ', num2str(diff_t9), ', Error: ', num2str(error_R_t9)])
disp(['  Simpsons:   ', num2str(I_s9), ', Diff: ', num2str(diff_s9), ', Error: ', num2str(error_R_s9)])

% --- Integral 10 --------------------------------------------------------
a10 = pi/3;
b10 = pi/2;
h10 = (b10 - a10) / (N - 1);
x10 = linspace(a10, b10, N);
y10 = (cos(x10).^3) ./ (sin(x10).^4);

[I_nr10, I_r10] = int_rect(x10, y10, h10, N);
[I_nt10, I_t10] = int_trap(x10, y10, h10, N);
[J_ns10, I_s10, u10] = int_simp(x10, y10, h10, N);
[I_na10, I_a10] = int_exact(x10, a10, b10, 10);

figure(10)
subplot(1, 2, 1)
plot(x10, y10, 'k')
grid
axis padded
subplot(1, 2, 2)
plot(x10, I_nr10, 'r', x10, I_nt10, 'g', u10, J_ns10, 'b', x10, I_na10, 'k')
grid
axis padded
legend('num. (rectangle)', 'num. (trapezoid)', 'num. (Simpson)', 'exact')

diff_r10 = I_r10 - I_a10;
diff_t10 = I_t10 - I_a10;
diff_s10 = I_s10 - I_a10;

R_r10 = I_nr10 - I_na10;
R_t10 = I_nt10 - I_na10;

error_R_r10 = sqrt(sum(R_r10.^2)/N);
error_R_t10 = sqrt(sum(R_t10.^2)/N);

J_na10 = I_na10(1:2:N);
R_s10 = J_ns10 - J_na10;
error_R_s10 = sqrt(sum(R_s10.^2)/N);

disp('[Integral 10]')
disp(['  Analytical: ', num2str(I_a10)])
disp(['  Rectangles: ', num2str(I_r10), ', Diff: ', num2str(diff_r10), ', Error: ', num2str(error_R_r10)])
disp(['  Trapezoids: ', num2str(I_t10), ', Diff: ', num2str(diff_t10), ', Error: ', num2str(error_R_t10)])
disp(['  Simpsons:   ', num2str(I_s10), ', Diff: ', num2str(diff_s10), ', Error: ', num2str(error_R_s10)])

%
% --- Functions ----------------------------------------------------------
%

function [I_n, I] = int_rect(x, y, h, N)
    
    I_n(1) = 0;
    I = 0;

    for n = 1:1:(N-1)
        I_n(n+1) = I_n(n) + y(n) * h;
        I = I + y(n) * h;
    end

end

function [I_n, I] = int_trap(x, y, h, N)
   
    I_n(1) = 0;
    I = 0;
    
    for n = 1:1:(N-1)
        I_n(n+1) = I_n(n) + ((y(n) + y(n+1)) / 2) * h;
        I = I + ((y(n) + y(n+1)) / 2) * h;
    end

end

function [J_n, I, u] = int_simp(x, y, h, N)

    J_n(1) = 0;
    I = 0;

    for n = 2:2:(N-1)
        J_n(n/2+1) = J_n(n/2) + (y(n-1) + 4*y(n) + y(n+1)) * h/3;
        I = I + (h/3)*(y(n-1) + 4*y(n) + y(n+1));
    end

    u = x(1:2:N);

end

function [I_n, I] = int_exact(x, a, b, int_num)

    if (int_num == 1)
        
        alpha = 1 + 2i;
        beta = 1 - 2i;
        chi = (3 + alpha) / (alpha - beta);
        xi = (3 + beta) / (alpha - beta);

        I = chi * log((b-alpha) / (a-alpha)) - xi * log((b-beta) / (a-beta));
        I_n = chi * log((x-alpha) / (a-alpha)) - xi * log((x-beta) / (a-beta));

    elseif (int_num == 2)

        I = (-5*atanh((b + 2) / sqrt(b^2 + 4*b + 10))) + (5*sqrt(b^2 + 4*b + 10)) + (5*atanh((a + 2) / sqrt(a^2 + 4*a + 10))) - (5*sqrt(a^2 + 4*a + 10));
        I_n = (-5*atanh((x + 2) / sqrt(x.^2 + 4*x + 10))) + (5*sqrt(x.^2 + 4*x + 10)) + (5*atanh((a + 2) / sqrt(a^2 + 4*a + 10))) - (5*sqrt(a^2 + 4*a + 10));

    elseif (int_num == 3)

        I = b*atan(b) - a*atan(a) - (1/2)*log((1+b^2) / (1+a^2));
        I_n = x.*atan(x) - a*atan(a) - (1/2)*log((1+x.^2) / (1+a^2));

    elseif (int_num == 4)

        I = exp(b)*(b^2 - 2*b + 2) - exp(a)*(a^2 - 2*a + 2);
        I_n = exp(x).*(x.^2 - 2*x + 2) - exp(a)*(a^2 - 2*a + 2);

    elseif (int_num == 5)

        I = (b/2)*cos(2*b) + ((b^2)/2)*sin(2*b) - (1/4)*sin(2*b) - (a/2)*cos(2*a) - ((a^2)/2)*sin(2*a) + (1/4)*sin(2*a) + (7/4)*cos(2*b) + (7*b/2)*sin(2*b) - (7/4)*cos(2*a) - (7*a/2)*sin(2*a) - (5/2)*sin(2*b) + (5/2)*sin(2*a);
        I_n = (x/2).*cos(2*x) + ((x.^2)/2).*sin(2*x) - (1/4)*sin(2*x) - (a/2)*cos(2*a) - ((a^2)/2)*sin(2*a) + (1/4)*sin(2*a) + (7/4)*cos(2*x) + (7*x/2).*sin(2*x) - (7/4)*cos(2*a) - (7*a/2)*sin(2*a) - (5/2)*sin(2*x) + (5/2)*sin(2*a);

    elseif (int_num == 6)

        I = (b/2)*sqrt(4-b^2) + 2*asin(b/2) - (a/2)*sqrt(4-a^2) - 2*asin(a/2);
        I_n = (x/2).*sqrt(4-x.^2) + 2*asin(x/2) - (a/2)*sqrt(4-a^2) - 2*asin(a/2);

    elseif (int_num == 7)

        I = atan(b/sqrt(2))/(2*sqrt(2)) - b/(4 + 2*b^2) - atan(a/sqrt(2))/(2*sqrt(2)) + a/(4 + 2*a^2);
        I_n = atan(x/sqrt(2))/(2*sqrt(2)) - x./(4 + 2*x.^2) - atan(a/sqrt(2))/(2*sqrt(2)) + a/(4 + 2*a^2);

    elseif (int_num == 8)

        I = (1/4) * (-log(b^2 + 1) + 2*log(1 - b) + 2*atan(b) + log(a^2 + 1) - 2*log(1 - a) - 2*atan(a));
        I_n = (1/4) * (-log(x.^2 + 1) + 2*log(1 - x) + 2*atan(x) + log(a^2 + 1) - 2*log(1 - a) - 2*atan(a));

    elseif (int_num == 9)

        I = (b/sqrt(2)) - (a/sqrt(2));
        I_n = (x/sqrt(2)) - (a/sqrt(2));

    elseif (int_num == 10)

        I = (1 / sin(b)) - (1 / sin(a)) + (1 / (3*(sin(a)^3))) - (1 / (3*(sin(b)^3)));
        I_n = (1 ./ sin(x)) - (1 / sin(a)) + (1 / (3*(sin(a)^3))) - (1 ./ (3*(sin(x).^3)));

    end

end