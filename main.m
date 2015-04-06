%% Deel 1: Drie veeltermbasissen

clc
format short

% equidistant points
nbpoints = 10;
x = linspace(-1,1,nbpoints);
c = [1:10]

% Chebychev evaluation
[T, dT] = chebyshev(x,5);

% Check if Chebyshev evaluation == cos(k bgcos(x))
[T, dT] = chebyshev_check(x,5);

%[y, dy] = chebpolyval(c,x)

%% Deel 2: Veelterminterpolatie

%% Deel 3: Methode van Newton-Raphson