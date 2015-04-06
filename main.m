%% Deel 1: Drie veeltermbasissen

clc
format short

% equidistant points
nbpoints = 10000;
x = linspace(-1,1,nbpoints);
c = [1:10];

% Chebychev evaluation
[T, dT] = chebyshev(x,5);

% Check if Chebyshev evaluation == cos(k bgcos(x))
[T, dT] = chebyshev_check(x,5);

% Plot Chebychev polynomials
figure(1),clf
plot(x,T)
hold on
title('Chebyshev veeltermen t.e.m. graad 5 in [-1,1]')
legend('graad 0','graad 1', 'graad 2', 'graad 3', 'graad 4', 'graad 5')
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('Chebyshev veelterm evaluaties') % y-axis label

[y, dy] = chebpolyval(c,x)

%% Deel 2: Veelterminterpolatie

%% Deel 3: Methode van Newton-Raphson