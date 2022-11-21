% PHYS 6352: Computational Physics
% Homework 5
% Author: Richard Camuccio

% [Introduction]
%
% In calculus we have an arbitrary function y = f(x) defined in the
% interval [a, b]. The corresponding thing on the computer will be a set of
% points (x_n, y_n) where n = 1,2,3,...,N. The x_n array can equally spaced
% or just general. If it is equally spaced then we will have a fixed step
% size: h = x_(n+1) - x_n. But the computers will work even if the space
% between the points x_n is not equal.
%
% If we know the function, or if we have a formula for y = f(x), we can use
% discrete sampling to put x_n in it and get y_n out. Often in practice we
% do not know this formula or it does not exist. In this case we cannot use
% the techniques from calculus but the computers don't care and they
% proceed just as well. Therefore, we will write the algorithms using only
% the discrete values x_n and y_n. However, to test our algorithms we can
% use known functions y = f(x) and obtain the exact values for y_n using
% y_n = f(x_n).

clc

% [Problem 1]
%
% One of the simplest methods of resampling is based on adding mid-points
% to the originally sampled data. The algorithm works as follows. Take the
% first pair of points, x_1 and x_2, and define the mid-point,
%
% x* = (x_1 + x_2) / 2
%
% Now you have to assign a value for y corresponding to this new point.
% Denote it by y*. You can approximate this value using linear
% interpolation between y_1 and y_2. Obviously, the answer will be
%
% y* = (y_1 + y_2) / 2
%
% Now we can move on to the next pair: x_2 and x_3 and repeat the process.
% Keep doing it until you reached the end of the x-sequence.
%
% At the end we will have 2N-1 points of which N were original samples and
% N-1 were calculated by linear interpolation. For large N this process
% approximately doubles the number of points.
%
% Write an M-script that implements this algorithm. Test it on some known
% function, say sin(x) or cos(x). Plot the exact and interpolated functions
% on the same plot. Also calculate the error that comes from linear
% interpolation using the formula we used in class.

clear

N = 10;
a = 0;
b = 2*pi;

x = linspace(a, b, N);
y = sin(x);

figure(1)
plot(x, y, 'k')
grid
axis padded
hold on

for n = 1:1:(N-1)

    xa(n) = (x(n) + x(n+1)) / 2;

    ya(n) = (y(n) + y(n+1)) / 2;

end

xn = [x, xa]
yn = [y, ya]

plot(xn, yn, 'r*')

% [Problem 2]
%
% A natural extension of the previous algorithm would be to add 2 points to
% each segment [x_n, x_(n+1)] instead of one. These new points can be
% equally spaced. In this case you will add two points to the small segment
% [x_1, x_2]; then move to the next segment [x_2, x_3], and so on.
%
% However from the computer's point of view this would not be the best
% approach. Instead, we should add m points to each segment [x_n, x_(n+1)].
% For each new point, we can approximate the corresponding value of y using
% linear interpolation.
%
% Here, the parameter m can take values 2,3,4, etc. Note that for m = 1 you
% should have the same answer as in Problem 1.
%
% Write an M-script that implements that algorithm. Test it on some known
% function, say sin(x) or cos(x). Plot the exact and interpolated functions
% on the same plot. Also calculate the error that comes from linear
% interpolation using the formula we used in class.