% PHYS 6352: Computational Physics
% Homework 3
% Author: Richard Camuccio

% [Introduction]
%
% Consider an arbitrary continuous function y = f(x) where x is on the
% interval [a, b]. In MATLAB we have to discretize the interval of x and
% obtain x_n for n = 1,2,3,...,N. We then find the correpsonding values of
% y and obtain y_n = f(x_n). At this point we can plot it and check
% visually if it has zero crossings. If we see the zero crossing then we
% can find the small segments [x_n, x_(n+1)] where these zero crossings
% take place. That was the first and simplest approach we did (Method 1).
%
% Later we came up with a more sophisticated technique that was based on
% splitting the interval into halves. This method was giving us a much more
% accurate result (Method 2). Method 2 does not work if the function has
% the same sign at both ends of the interval. However, it still can cross
% zero inside the interval. For example, y = x^2 - 1 defined in interval
% [-4, 4] has two zeros: x = -1 and x = 1 inside the interval and has the
% same sign at the boundaries.

clc

% [Problem 1]
%
% Combine Method 1 and Method 2 to obtain an accurate method of finding the
% zeros of an arbitrary function. Practice your method on a simple example:
%
% y = sin(x) + (1/2), defined in x : [0, 6*pi],
%
% and show that it works. Check the accuracy of your algorithm by comparing
% it with the analytical result.

clear

N = 100;

a = 0;
b = 6*pi;

x = linspace(a, b, N);
y = sin(x) + (1/2);

figure(1)
plot(x, y, 'k')
grid
axis padded
hold on

Xs = [];
Ys = [];

for n = 1:(N-1)

    prod = y(n) * y(n+1);

    if prod < 0

        Xs = x(n);
        Ys = y(n);

        plot(Xs, Ys, 'r*')
        
    end

end

% [Problem 2]
%
% Complete the M-script that animates the expansion of the cycloid as the
% wheel is rolling from left to right. Add the wheel to the figure.
%
% The parameters for the cycloid are: radius of the wheel (r) and angle of
% rotation (phi) defined in the interval [a, b]. We can always set a = 0,
% but you have to choose b sufficiently large.

a = 0;
b = 4*pi;
h = (b - a) / (N - 1);
x = linspace(a, b, N);

radius = 1;
phi(1) = 0;

% center of wheel (x, y) = (ac, bc)
ac = 0;
bc = radius;
phi_c = linspace(0, 2*pi, N);

for n = 1:1:N

    % cycloid position
    x = radius * phi - radius * sin(phi);
    y = radius - radius * cos(phi);

    % end point position
    xe = x(end);
    ye = y(end);

    % calculate circle position
    xc = radius^2 - (bc + radius*sin(phi_c));
    yc = radius^2 - (ac + radius*cos(phi_c));

    % plot with equal axis
    figure(2)
    plot(x, y, 'k', xe, ye, 'r*')
    hold on
    plot(xc, yc, 'b')
    hold off
    axis equal
    axis([a*radius b*radius 0 2.5])

    % jog animation
    phi(n+1) = phi(n) + h;
    bc = bc - h;
    pause(0.01)

end

% [Problem 3]
%
% Modify your previous M-script to animate expansion of the epicycloid as
% the wheel is rolling counterclockwise on the surface of the disk.
%
% The parameters for the epicycloid are: radius of the wheel (r), radius of
% the disk (R), and angle of rotation (phi) defined in the interval [a, b].
%
% Can you figure out under which conditions the trajectory closes onto
% itself?

% [Problem 4]
%
% Modify your previous M-script to animate expansion of the hypocycloid as
% the wheel is rolling counterclockwise inside the larger circle. The
% parameters for the hypocycloid are: radius of the wheel (r), radius of
% the larger circle (R), and phi.
%
% Can you figure out under which conditions the trajectory closes onto
% itself?

% [Problem 5]
%
% Using your previous M-script observe Copernicus Theorem: if the radius of
% the wheel is half the radius of the larger circle the hypocycloid becomes
% a straight line.
%
% Can you figure out why this is happening? Or, can you prove this fact
% mathematically?