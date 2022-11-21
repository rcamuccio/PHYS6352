% PHYS 6352: Computational Physics
% Homework 1
% Author: Richard Camuccio

% [Problem 1]
% In many grading schemes the final grade is calculated from the overall
% performance in percents using the following correspondence: grade "A" is
% given for 90-100%, grade "B" for 80-90%, and so on.
%
% Write an M-function that takes the percentile performance as input
% variable and produces the output in the form of a statement: "The final
% grade for this class is X", where "X" stands for one of the letters: "A",
% "B", ....

clear

test_percentage = 83;

final_grade = percentage(test_percentage);

disp('[Problem 1]')
disp(['  The final grade for this class is ', final_grade])
disp(' ')

% [Problem 2]
% A function is defined as follows
%
% f(x) = 0 (if x < 0) and f(x) = x (if x >= 0)
% 
% (a) Using MATLAB commands "if" and "else" write an M-function for it. Use
% a for-loop to calculate N = 100 values of this function in the interval
% x : [-4, 4] and plot it.
%
% (b) Derive a single analytical expression for this function using the
% absolute value sign. Compare the M-function in part (a) with your
% analytical expression by plotting both within an M-script. (Note: in
% MATLAB the absolute value is realized by "abs".)

clf
figure(1)
grid
axis equal
hold on

for x = -4:0.08:4

    y = strange(x);
    plot(x, y, 'r.')

end

xl = linspace(-4, 4, 100);
yl = (xl + abs(xl)) / 2;
plot(xl, yl, 'k')
title('Figure 1')

hold off

disp('[Problem 2]')
disp('  See Figure 1')
disp(' ')

% [Problem 3]
% A function is defined as follows:
%
% f(x) = sqrt(x) for x > 0 and f(x) = -sqrt(-x) for x <= 0
%
% (a) Using MATLAB commands "if" and "else" write an M-function for it. Use
% a for-loop to calculate N = 100 values of this function in the interval
% x : [-10, 10] and plot it.
%
% [b] Derive a single analytical expression for this function using the
% absolute value sign. Compare the M-function in part [a] with your
% analytical expression by plotting both within an M-script. (Note: you may
% need to treat the point x = 0 separately.)

figure(2)
grid
axis equal
hold on

for x = -10:0.1:10

    y = charm(x);
    plot(x, y, 'r.')

end

xl = linspace(-10, 10, 100);
yl = abs(sqrt(xl)) .* sign(xl);
plot(xl, yl, 'k')
title('Figure 2')

hold off

disp('[Problem 3]')
disp('  See Figure 2')
disp(' ')

% [Problem 4]
% A function is defined as follows:
%
% f(x) = sqrt(x - 5), x >= 5
% f(x) = 25 - x^2, -5 < x < 5
% f(x) = sqrt(-x - 5), x <= -5
%
% Write an M-function which describes it. (Use "if", "elseif", and "else"
% commands.) Using a for-loop calculate N = 1000 values of this function in
% the interval x : [-20, 20] and plot it.

figure(3)
grid
axis equal
hold on

for x = -20:0.04:20

    y = beauty(x);
    plot(x, y, 'r.')

end

title('Figure 3')
hold off

disp('[Problem 4]')
disp('  See Figure 3')
disp(' ')

% [Problem 5]
% Use a while-loop to determine how many terms (N) in the series
% a_n = 2^n + 3^n, n = 1, 2, 3, ... are required for the sum of the terms
%
% S = sum(a_n) from n = 1 to N
%
% to exceed A = 10^6. Think of a method to check your calculations and
% implement the check procedure in your code to show that your answer is
% correct.

A = 10^6;
S = 0;
N = 0;

a = 2;
b = 3;

while (S < 10^6)

    N = N + 1;
    S = S + ((a^N) + (b^N));

end

T_1 = a * (1 - a^N) / (1 - a);
T_2 = b * (1 - b^N) / (1 - b);
T = T_1 + T_2;

check = S - T;

disp('[Problem 5]')
disp(['  ', num2str(S), ' is the sum'])
disp(['  ', num2str(N), ' terms'])
disp(['  Check: ', num2str(check)])
disp(' ')

% [Problem 6]
% Write an M-script to calculate how many years it will take to accumulate
% more than $100000 in a savings account in a bank if the initial deposit
% is $1000 and the annual interest rate is 5.547 percent. (Use a
% while-loop.)

p = 1000;
A = 100000;
r = 0.05547;

yr = 0;

while (p < A)

    yr = yr + 1;

    p = p + (p * 0.05547);

end

disp('[Figure 6]')
disp(['  You will have $', num2str(A), ' after ', num2str(yr), ' years'])
disp(' ')

% [Problem 7]
% Write an M-script for calculating interest on a loan from a bank. Let the
% initial sum borrowed be S = $18000. Assume that the loan has a fixed
% annual interest rate of 6% or r = 0.06 (the rate per month is f = r /12)
% and also assume that the person pays a fixed amount every month p = $200.
% How many months will it take for this person to repay the loan? What is
% the total amount he will pay to the bank, and what is the total interest
% paid? (Use either for-loop or while-loop.)

S = 18000;
r = 0.06;
n = 12;
f = r / n;
p = 200;

A = 0; % goal
m = 0; % months

p_tot = 0; % total amount paid
p_int = 0; % total interest paid

while (S > A)

    m = m + 1;

    cur_int = S * f;

    S = S + cur_int;
    S = S - p;

    p_int = p_int + cur_int;
    p_tot = p_tot + p;

end

disp('[Problem 7]')
disp(['  Months until paid: ', num2str(m)])
disp(['  Total amount paid: $', num2str(p_tot)])
disp(['  Total interest paid: $', num2str(p_int)])
disp(' ')

% --- Functions

function y = percentage(x)

    if (x < 60)

        y = 'F';

    elseif (x < 70)

        y = 'D';

    elseif (x < 80)

        y = 'C';

    elseif (x < 90)

        y = 'B';

    else

        y = 'A';

    end

end

function y = strange(x)

    if (x < 0)

        y = 0;

    else

        y = x;

    end

end

function y = charm(x)

    if (x <= 0)

        y = -sqrt(-x);

    else

        y = sqrt(x);
    
    end

end

function y = beauty(x)

    if (x <= -5)

        y = sqrt(-x - 5);

    elseif (x < 5)

        y = 25 - x^2;

    else

        y = sqrt(x - 5);

    end

end