% PHYS 6352: Computational Physics
% Homework 8
% Author: Richard Camuccio

clear
clc

% [Problem 1]
%
% The period of oscillations of a simple pendulum (mass on a string) is
% given by the formula
%
% T = 2*pi*sqrt(L/g)
%
% where L is the length of the string and g is the acceleration of free
% fall. Thus, by measuring the period and the length one can find the value
% for g.
%
% A somewhat different approach would be to measure the period for
% different values of the length to see how well T = 2*pi*sqrt(L/g) is
% satisfied. In this case, it is more convenient to use the following form
% of this equation:
%
% L = (g/(4*pi^2))*T^2
%
% Assuming that the x-data is period (T_i) and the y-data is length (L_i)
% we have the following record:

x = [0.91, 1.11, 1.23, 1.42, 1.55, 1.70, 1.84];
y = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];

% (a) Fit the following curve to the data by hand as best you can
%
% y = c*x^2
%
% by changing the value of parameter c and looking at R. In other words,
% try to obtain the value for R as low as you can.

c = 0.2487;
yfit = c*x.^2;

dy = yfit - y;
F = sum(dy.^2);
R = sqrt(F);

g = 4*pi*pi*c;

disp('[Problem 1]')
disp(['  Estimated  R: ', num2str(R)])
disp(['             c: ', num2str(c), ' m/s^2'])
disp('             --------------------')
disp(['             g: ', num2str(g), ' m/s^2'])
disp(' ')

figure(1)
subplot(1,2,1)
plot(x, y, 'ro', x, yfit, 'b')
grid
xlabel('Period [s]')
ylabel('Length [m]')
title('Estimated fit to L = cT^2')
legend('data', 'fit', 'Location', 'NorthWest')

% (b) Calculate the optimal value for the parameter c analytically from the
% condition
%
% dR/dc = 0
%
% and obtain the least squares fit to the data. Convince yourself that the
% value of R for the optimal parameter is smaller than the one you obtained
% by tuning the parameter by hand in part (a).

c = sum(y / x.^2);
yfit = c*x.^2;

dy = yfit - y;
F = sum(dy.^2);
R = sqrt(F);

g = 4 * pi * pi * c;

disp(['  Analytical R: ', num2str(R)])
disp(['             c: ', num2str(c), ' m/s^2'])
disp('             --------------------')
disp(['             g: ', num2str(g), ' m/s^2'])

subplot(1,2,2)
plot(x, y, 'ro', x, yfit, 'b')
grid
xlabel('Period [s]')
ylabel('Length [m]')
title('Analytical fit to L = cT^2')
legend('data', 'fit', 'Location', 'NorthWest')

% [Problem 2]
%
% The force of viscous friction F acting on an object moving thruogh medium
% is proportional to its velocity v for slow motion and proportional to v^2
% for fast motion. Mathematically, this can be expressed by the formula
%
% F = k1*v + k2*v^2
%
% where k1 and k2 are some constants. The characteristic velocity at which
% the quadratic part dominates is given by v0 = k1/k2. Assuming that the
% x-data is velocity (v_i) and the y-data is force (F_i) we have the
% following record:

x = [8, 12, 16, 20, 24, 28];
y = [13, 17, 25, 42, 47, 62];

% (a) Fit the following curve to the data by hand as best as you can
%
% y = b*x + c*x^2
%
% by changing the value of parameters b and c and looking at R. In other
% words, try to obtain the value for R as low as you can.

b = 1.3893;
c = 0.0295;
yfit = b*x + c*x.^2;

dy = yfit - y;
F = sum(dy.^2);
R = sqrt(F);

v0 = b/c;

disp('[Problem 2]')
disp(['  Estimated  R: ', num2str(R)])
disp(['             b: ', num2str(b), ' kg/s'])
disp(['             c: ', num2str(c), ' kg/m'])
disp('             --------------------')
disp(['             v0: ', num2str(v0), ' m/s'])
disp(' ')

figure(2)
subplot(1,2,1)
plot(x, y, 'ro', x, yfit, 'b')
grid
xlabel('Velocity [m/s]')
ylabel('Force [N]')
title('Estimated fit to F = bv + cv^2')
legend('data', 'fit', 'Location', 'NorthWest')

% (b) Calculate the optimal values for parameters b and c analytically from
% the conditions
%
% dR/db = 0 and dR/dc = 0
%
% and thus obtain the least squares fit to the data. Convince yourself that
% the value of R for the optimal parameters is smaller than the one you
% obtained by tuning the parameters by hand.

v = [sum(x.*y); sum(y.*x.^2)];
M = [sum(x.^2), sum(x.^3); sum(x.^3), sum(x.^4)];
u = inv(M) * v;

b = u(1);
c = u(2);

yfit = b*x + c*x.^2;

dy = yfit - y;
F = sum(dy.^2);
R = sqrt(F);

v0 = b/c;

disp(['  Analytical R: ', num2str(R)])
disp(['             b: ', num2str(b), ' kg/s'])
disp(['             c: ', num2str(c), ' kg/m'])
disp('             --------------------')
disp(['             v0: ', num2str(v0), ' m/s'])

subplot(1,2,2)
plot(x, y, 'ro', x, yfit, 'b')
grid
xlabel('Velocity [m/s]')
ylabel('Force [N]')
title('Analytical fit to F = bv + cv^2')
legend('data', 'fit', 'Location', 'NorthWest')

% [Problem 3]
%
% A cat fell from the balcony of tall building and its fall was captured by
% digital camera with a clock. Neglecting the air friction, we can describe
% the cat's height with respect to the ground by the formula:
%
% h(t) = h0 + v0*t + (1/2)*beta*t^2
%
% where h0, v0, and beta are constants. Assuming that the x-data is time
% (t_i) and the y-data is height (h_i) we have the following record:

x = [0.0, 0.4, 0.8, 1.2, 1.6];
y = [27.4, 22.1, 17.3, 9.9, 1.3];

% (a) Fit the following curve to the data by hand as best as you can
%
% y = a + bx + cx^2
%
% by changing the value of parameters a, b, and c and looking at R. In
% other words, try to obtain the value for R as low as you can.

a = 27.4;
b = -8.9375;
c = -4.61081;
yfit = a + b*x + c*x.^2;

dy = yfit - y;
F = sum(dy.^2);
R = sqrt(F);

beta = 2*c;

disp('[Problem 3]')
disp(['  Estimated  R: ', num2str(R)])
disp(['             a: ', num2str(a), ' m'])
disp(['             b: ', num2str(b), ' m/s'])
disp(['             c: ', num2str(c), ' m/s^2'])
disp('             --------------------')
disp(['             h0: ', num2str(a), ' m'])
disp(['             v0: ', num2str(b), ' m/s'])
disp(['             beta: ', num2str(beta), ' m/s^2'])
disp(' ')

figure(3)
subplot(1,2,1)
plot(x, y, 'ro', x, yfit, 'b')
grid
xlabel('Time [s]')
ylabel('Height [m]')
title('Estimated fit to h = a + bt + ct^2')
legend('data', 'fit', 'Location', 'NorthEast')

% (b) Calculate the optimal values for parameters a, b, and c analytically
% from the conditions
%
% dR/da = 0, dR/db = 0, dR/dc = 0
%
% and thus obtain the least squares fit to the data. Convince yourself that
% the value for R for the optimal parameters is smaller than the one you
% obtained by tuning the parameters by hand.

N = length(x);
v = [sum(y); sum(x.*y); sum(y.*x.^2)];
M = [N, sum(x), sum(x.^2); sum(x), sum(x.^2), sum(x.^3); sum(x.^2), sum(x.^3), sum(x.^4)];
u = inv(M) * v;

a = u(1);
b = u(2);
c = u(3);
yfit = a + b*x + c*x.^2;

dy = yfit - y;
F = sum(dy.^2);
R = sqrt(F);

beta = 2*c;

disp(['  Analytical R: ', num2str(R)])
disp(['             a: ', num2str(a), ' m'])
disp(['             b: ', num2str(b), ' m/s'])
disp(['             c: ', num2str(c), ' m/s^2'])
disp('             --------------------')
disp(['             h0: ', num2str(a), ' m'])
disp(['             v0: ', num2str(b), ' m/s'])
disp(['             beta: ', num2str(beta), ' m/s^2'])

subplot(1,2,2)
plot(x, y, 'ro', x, yfit, 'b')
grid
xlabel('Time [s]')
ylabel('Height [m]')
title('Analytical fit to h = a + bt + ct^2')
legend('data', 'fit', 'Location', 'NorthEast')

% (d) Assuming that the cat fell without initial velocity find the height
% of the balcony. Note that the beginning of the recording (t = 0) is
% arbitrary and may not coincide with the beginning of the fall. Also, find
% the total time the cat was in free fall.

x = linspace(-5, 5, 100);
y = a + b*x + c*x.^2;

xp = -b/(2*c);
yp = a + b*xp + c*xp^2;

t = sqrt(abs(a/c));

disp(' ')
disp(['  Initial height: ', num2str(yp), ' meters'])
disp(['  Time of fall: ', num2str(t), ' seconds'])