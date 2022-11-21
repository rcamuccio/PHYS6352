% PHYS 6352: Computational Physics
% Homework 2
% Author: Richard Camuccio

clc

% [Problem 1]
%
% Viete's theorem states that a cubic polynomial can be uniquely defined by
% its roots x_1, x_2, x_3 and the coefficient a:
%
% ax^3 + bx^2 + cx + d = a(x - x_1)(x - x_2)(x - x_3)
%
% Viete gives explicit formulas for coefficients b, c, and d in terms of a
% and x1, x2, x3. Viete's formulas can be very useful if we want to know
% the answer to the question: "What kind of cubic polynomial has the roots
% that we specified?"
%
% Write an M-script (or M-function) that calculates the coefficients b, c,
% d from the known roots and the coefficient a that you must specify by
% hand. You can use this M-script to test your calculations in Problem 3.

clear

a = 1;
x1 = 2;
x2 = 4;
x3 = 6;

b = -a * (x1 + x2 + x3);
c = a * ((x1 * x2) + (x1 * x3) + (x2 * x3));
d = - a * x1 * x2 * x3;

disp('[Problem 1]')
disp(['  Given a = ', num2str(a), ' and roots ', num2str(x1), ', ', num2str(x2), ', ', num2str(x3)])
disp(['   b = ', num2str(b)])
disp(['   c = ', num2str(c)])
disp(['   d = ', num2str(d)])
disp(' ')

% [Problem 2]
%
% Beginning with a general expression for cubic equation:
%
% ax^3 + bx^2 + cx + d = 0
%
% perform the operations to bring it to the reduced form:
%
% y^3 + py + q = 0
%
% Give explicit formulas for the coefficients p and q in terms of a, b, c,
% d.

clear

a = 1;
b = 3;
c = 5;
d = 7;

p = (c / a) - (b^2 / (3 * a^2));
q = ((2/27) * (b^3 / a^3)) - ((b * c) / (3 * a^2)) + (d / a);

disp('[Problem 2]')
disp('  Given')
disp(['   a = ', num2str(a)])
disp(['   b = ', num2str(b)])
disp(['   c = ', num2str(c)])
disp(['   d = ', num2str(d)])
disp('  Reduced')
disp(['   p = ', num2str(p)])
disp(['   q = ', num2str(q)])
disp(' ')

% [Problem 3]
%
% Write an M-script (or M-function) that calculates the roots of the
% reduced cubic equation:
%
% y^3 + py + q = 0
%
% Combine it with the solution of Problem 2 to have an M-script that solves
% a general cubic equation:
%
% ax^3 + bx^2 + cx + d = 0
%
% Include in your M-script checks of the results. Also include in your
% M-script a plotting routine that produces a graph of the function
%
% f(x) = ax^3 + bx^2 + cx + d
%
% to convince yourself that the number of real roots you found is correct.
% You can overlay your roots with the graph of this function.

clear

N = 1e6;
xl = -10;
xr = 10;
x = linspace(xl, xr, N);

a = 1;
b = 0;
c = -25;
d = -10;

p = (c / a) - (b^2 / (3 * a^2));
q = ((2/27) * (b^3 / a^3)) - ((b * c) / (3 * a^2)) + (d / a);

alpha = (-1 + 1i*sqrt(3)) / 2;
beta = (-1 - 1i*sqrt(3)) / 2;

u = ((-q/2) + sqrt(((q^2)/4 + (p^3)/27)))^(1/3);
v = ((-q/2) - sqrt(((q^2)/4 + (p^3)/27)))^(1/3);

x1 = u + v;
x2 = (u * alpha) + (v * beta);
x3 = (u * beta) + (v * alpha);

y1 = (a * x1.^3) + (b * x1.^2) + (c * x1) + d;
y2 = (a * x2.^3) + (b * x2.^2) + (c * x2) + d;
y3 = (a * x3.^3) + (b * x3.^2) + (c * x3) + d;

y_cubic = (a * x.^3) + (b * x.^2) + (c * x) + d;

xr1 = real(x1);
xr2 = real(x2);
xr3 = real(x3);

yr1 = real(y1);
yr2 = real(y2);
yr3 = real(y3);

disp('[Problem 3]')
disp('  Given starting values for a cubic polynomial')
disp(['   a = ', num2str(a)])
disp(['   b = ', num2str(b)])
disp(['   c = ', num2str(c)])
disp(['   d = ', num2str(d)])
disp('  We have the roots')
disp(['   x1 = ', num2str(x1), ', y1 = ', num2str(y1)])
disp(['   x2 = ', num2str(x2), ', y2 = ', num2str(y2)])
disp(['   x3 = ', num2str(x3), ', y3 = ', num2str(y3)])
disp(' ')

figure(1)
grid
hold on
plot(x, y_cubic, 'k')
plot(xr1, yr1, 'b*', xr2, yr2, 'r*', xr3, yr3, 'g*')

% [Problem 4]
%
% Viete's theorem states that a quartic polynomial can be uniquely defined
% by its roots x1, x2, x3, x4 and the coefficient a:
%
% ax^4 + bx^3 + cx^2 + dx + e = a(x - x1)(x - x2)(x - x3)(x - x4)
%
% Viete gives explicit formulas for coefficients b, c, d, and e in terms of
% a and x1, x2, x3, x4. Viete's formulas can be very useful if we want
% to know the answer to the question: "What kind of quartic polynomial has
% the roots that we specified?"
%
% Write an M-script (or M-function) that calculates the coefficients b, c,
% d, e from the known roots and the coefficient a that you must specify by
% hand. You can use this M-script to test your calculations in Problem 6.

clear

a = 1;
x1 = 3;
x2 = 5;
x3 = 7;
x4 = 9;

b = -a * (x1 + x2 + x3 + x4);
c = a * ((x1 * x2) - (x1 * x3) - (x2 * x3) + (x1 * x4) + (x2 * x4) + (x3 * x4));
d = a * ((-x1 * x2 * x3) - (x1 * x2 * x4) + (x1 * x3 * x4) + (x2 * x3 * x4));
e = a * x1 * x2 * x3 * x4;

disp('[Problem 4]')
disp(['  Given a = ', num2str(a), ' and roots ', num2str(x1), ', ', num2str(x2), ', ', num2str(x3), ', ', num2str(x4)])
disp(['   b = ', num2str(b)])
disp(['   c = ', num2str(c)])
disp(['   d = ', num2str(d)])
disp(['   e = ', num2str(e)])
disp(' ')

% [Problem 5]
%
% Beginning with a general expression for quartic equation:
%
% ax^4 + bx^3 + cx^2 + dx + e = 0
%
% perform the operations to bring it to the reduced form:
%
% y^4 + py^2 + qy + r = 0
%
% Give explicit formulas for the coefficients p, q, and r in terms of a, b,
% c, d, e

clear

a = 1;
b = 1;
c = 2;
d = 3;
e = 5;

bp = b/a;
cp = c/a;
dp = d/a;
ep = e/a;

p = cp - (6/16)*(bp^2);
q = (((8 * (bp^3)) - (32 * bp * cp)) / 64) + dp;
r = ((bp^4) / 256) - ((bp^3) / 64) + ((bp^2) * cp / 16) - (bp * dp / 4) + cp;

disp('[Problem 5]')
disp('  Given')
disp(['   a = ', num2str(a)])
disp(['   b = ', num2str(b)])
disp(['   c = ', num2str(c)])
disp(['   d = ', num2str(d)])
disp(['   e = ', num2str(e)])
disp('  Reduced')
disp(['   p = ', num2str(p)])
disp(['   q = ', num2str(q)])
disp(['   r = ', num2str(r)])
disp(' ')

% [Problem 6] (Bonus)
%
% Write an M-script (or M-function) that calculates the roots of the
% reduced quartic equation:
%
% y^4 + py^2 + qx + r = 0
%
% Combine it with the solution of Problem 5 to have an M-script that solve
% a general quartic equation:
%
% ax^4 + bx^3 + cx^2 + dx + e = 0
%
% Include in your M-script checks of the results. Also include in your
% M-script a plotting routine that produces a graph of the function
%
% f(x) = ax^4 + bx^3 + cx^2 + dx + e
%
% to convince yourself that the number of real roots you found is correct.
% You can overlay your roots with the graph of this function.

% [Problem 7] Parabola
%
% In Cartesian coordinate system parabola is given by
%
% x^2 = 2py
%
% Plot this curve in MATLAB for different values of p. Find out the meaning
% of parameter p. Identify the focal point F (focus) and prove numerically
% that for every point P on the parabola |PF| = |PD|. Here D is the point
% on a straight line: y = -p/2 such that |PD| is the shortest distance.

N = 1e6;
a = -4;
b = 4;
x = linspace(a, b, N);

% semi-latus rectum - half of a chord going through focus of an ellipse;
% increasing p widens the parabola opening; p < 0 flips the parabola down
p1 = 1;
p2 = 2;
p3 = 3;

% focal length
fl1 = p1 / 2;
fl2 = p2 / 2;
fl3 = p3 / 2;

% directrix
d1 = -fl1;
d2 = -fl2;
d3 = -fl3;

% parabola
y1 = x.^2 / (2 * p1);
y2 = x.^2 / (2 * p2);
y3 = x.^2 / (2 * p3);

figure(2)
hold on

% plot parabola
plot(x, y1, 'b', x, y2, 'g', x, y3, 'r')

% plot focus
plot(0, fl1, 'b*', 0, fl2, 'g*', 0, fl3, 'r*')

% plot directrix
yline(d1, 'b')
yline(d2, 'g')
yline(d3, 'r')

% plot points P (on parabola) and D (on directrix)
plot(p1, fl1, 'bx', 0, -fl1, 'bx')
plot(p2, fl2, 'gx', 0, -fl2, 'gx')
plot(p3, fl3, 'rx', 0, -fl3, 'rx')

% plot line segment |PF1| and |PD1| respectively
plot([0, p1], [fl1, fl1], 'k')
plot([0, 0], [d1, fl1], 'k:')

legend('p=1', 'p=2', 'p=3')
grid

% check if semi-latus rectum equals distance between focus and directrix
delta_1 = fl1 - d1;
check_1 = delta_1 - p1;

delta_2 = fl2 - d2;
check_2 = delta_2 - p2;

delta_3 = fl3 - d3;
check_3 = delta_3 - p3;

disp('[Problem 7]')
disp(['  Check for |PF1| = |PD1| (should be zero): ', num2str(check_1)])
disp(['  Check for |PF2| = |PD2| (should be zero): ', num2str(check_2)])
disp(['  Check for |PF3| = |PD3| (should be zero): ', num2str(check_3)])
disp(' ')

% [Problem 8] Ellipse
%
% In Cartesian coordinate system ellipse is given by
%
% (x^2 / a^2) + (y^2 / b^2) = 1
%
% where a and b are called semi-axes. For a = b ellipse becomes circle.
%
% Plot ellipse in MATLAB. Identify its two focal points F1 and F2 (foci)
% and prove numerically that the sum |PF1| + |PF2| is the same for every
% point P on the ellipse.

clear

N = 1e6;
xl = -10;
xr = 10;
x = linspace(xl, xr, N);

a = 3;
b = 2;
c = sqrt((a^2) - (b^2));

y1 = real(sqrt(b * (1 - (x.^2/a^2))));
y2 = real(-sqrt(b * (1 - (x.^2/a^2))));

y1(y1 == 0) = nan;
y2(y2 == 0) = nan;

figure(3)
hold on

plot(x, y1, 'k:', x, y2, 'k:')
plot(c, 0, 'kx', -c, 0, 'kx')

grid

% [Problem 9] Hyperbola
%
% In Cartesian coordinate system hyperbola is given by
%
% (x^2 / a^2) - (y^2 / b^2) = 1
%
% or
%
% -(x^2 / a^2) + (y^2 / b^2) = 1
%
% where a and b are called semi-axes.
%
% Plot hyperbola in MATLAB. Identify its two focal points F1 and F2
% (foci) and prove numerically that the difference ||PF1| - |PF2|| is the
% same for every point on the hyperbola.

clear

N = 1e6;
xl = -10;
xr = 10;
x = linspace(xl, xr, N);

a = 3;
b = 2;
c = sqrt((a^2) + (b^2));

y1 = real(b * sqrt(((x.^2)/(a^2)) - 1));
y2 = real(-b * sqrt(((x.^2)/(a^2)) - 1));

y1(y1 == 0) = nan;
y2(y2 == 0) = nan;

figure(4)
hold on

plot(x, y1, 'k:', x, y2, 'k:')
plot(c, 0, 'kx', -c, 0, 'kx')

grid