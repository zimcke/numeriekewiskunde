%% 3.1 Equidistante punten en het Runge fenomeen

clc
f = @fx;

nbpoints = 10000;
z = linspace(-1,1,nbpoints);
y = fx(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lagrange interpolations of the Runge function through equidistant nodes,  % 
% degrees 8, 10, 12 and 14                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lagrange interpolation degree 8
x = linspace(-1,1,9);
lagrange8 = evalueer_lagrange(x, f, z);

% Lagrange interpolation degree 10
x = linspace(-1,1,11);
lagrange10 = evalueer_lagrange(x, f, z);

% Lagrange interpolation degree 12
x = linspace(-1,1,13);
lagrange12 = evalueer_lagrange(x, f, z);

% Lagrange interpolation degree 14
x = linspace(-1,1,15);
lagrange14 = evalueer_lagrange(x, f, z);

% Figure 1: Lagrange interpolation
figure(1),clf
plot(z,y)
hold on
plot(z, lagrange8)
plot(z, lagrange10)
plot(z, lagrange12)
plot(z, lagrange14)
title('Veelterminterpolatie in de Lagrange basis');
legend('f(x)','Graad 8', 'Graad 10', 'Graad 12', 'Graad 14');
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('Functiewaarden interpolatie') % y-axis label
print('report\lagrange_interpolatie','-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relative error for Lagrange interpolation, when using equidistant nodes   %
%(n = 20,50) and Chebyshev nodes (n = 20,50)                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lagrange interpolation degree 20, equidistant nodes 
x = linspace(-1,1,21);
lagrange20 = evalueer_lagrange(x, f, z);
abs_error20 = abs(lagrange20 - y);
rel_error20 = abs_error20 ./ y;

% Lagrange interpolation degree 50, equidistant nodes 
x = linspace(-1,1,51);
lagrange50 = evalueer_lagrange(x, f, z);
abs_error50 = abs(lagrange50 - y);
rel_error50 = abs_error50 ./ y;

% Lagrange interpolation degree 20, Chebyshev nodes 
xc = chebyshev_nodes(21);
lagrange20c = evalueer_lagrange(xc, f, z);
abs_error20c = abs(lagrange20c - y);
rel_error20c = abs_error20c ./ y;

% Lagrange interpolation degree 50, Chebyshev nodes 
xc = chebyshev_nodes(51);
lagrange50c = evalueer_lagrange(xc, f, z);
abs_error50c = abs(lagrange50c - y);
rel_error50c = abs_error50c ./ y;

% Figure 2: relative error with Lagrange interpolation (equidistant and
% chebychev nodes)
figure(2),clf
semilogy(z, rel_error20)
hold on
semilogy(z, rel_error50)
semilogy(z, rel_error20c)
semilogy(z, rel_error50c)
title('Relatieve fout bij Langrange interpolatie');
legend('Equidistante punten n=20', 'Equidistante punten n=50', 'Chebyshev punten n=20', 'Chebyshev punten n=50', 'Location','southwest')
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('Relative fout') % y-axis label
print('report\relatieve_fout','-depsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation error for Lagrange interpolation, when using equidistant    %
% nodes and Chebyshev nodes, for degrees 0 through 50                       %             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolation error for equidistant nodes
E_equidistant = zeros(1,51);
z = linspace(-1,1,10000);
for i = 0:50
    x = linspace(-1,1,i+1);
    E_equidistant(i+1) = max(abs(fx(z) - evalueer_lagrange(x,f,z)));
end

% Interpolation error for Chebyshev nodes
E_chebyshev = zeros(1,51);
for i = 0:50
    Echebnodes = chebyshev_nodes(i+1);
    E_chebyshev(i+1) = max(abs(fx(z) - evalueer_lagrange(Echebnodes,f,z)));
    clear xc;
end

% Figure 3: Interpolation error
figure(3),clf
semilogy(0:50, E_equidistant);
hold on
semilogy(0:50, E_chebyshev);
title('Benaderingsfout E voor toenemende graden 0,1,...,50');
legend('Door equidistante punten', 'Door Chebyshev punten', 'Location', 'southwest')
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('E') % y-axis label
print('report\benaderingsfout_interpolatiepunten','-depsc')


%% 3.2 Verschillende basissen

clc

f = @fx;
nbpoints = 100;
z = linspace(-1,1,nbpoints);
E_lagrange = zeros(1,80);
E_monomial = zeros(1,80);
E_chebyshev = zeros(1,80);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation error for Lagrange, Monomial and Chebyshev interpolation,   % 
% for degrees 1 through 80                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:80
    warning('Off', 'MATLAB:nearlySingularMatrix');
    chebnodes = chebyshev_nodes(i+1);
    
    % calculate monomial coefficients
    M = monomiaal(chebnodes, i);
    b = f(chebnodes)';
    coeff_monomial = M\b;
    
    % calculate Chebyshev coefficients
    [T, dT] = chebyshev(chebnodes, i);
    coeff_chebyshev = T\b;
    
    % interpolation error
    E_lagrange(i) = max(abs(fx(z) - evalueer_lagrange(chebnodes,f,z)));
    E_monomial(i) = max(abs(fx(z) - polyval(coeff_monomial,z)));
    E_chebyshev(i) = max(abs(fx(z) - chebpolyval(coeff_chebyshev, z)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The matrix M above is known as the Vandermonde matrix and is notoriously  % 
% ill-conditioned, if the nodes are close to each other or the matrix is    %
% of high-order. In such a case the system may not be accurately solved.    % 
% => check with cond(M)                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 5: interpolation error for Lagrange, Monomial and Chebyshev
% interpolation
figure(5), clf
semilogy(1:80, E_lagrange);
hold on
semilogy(1:80, E_monomial);
semilogy(1:80, E_chebyshev);
title('Interpolatiefout bij verschillende basissen, graad 1 t.e.m 80')
legend('E Lagrange','E monomiaal', 'E Chebyshev', 'Location','southwest');
print('report\verschillende_basissen','-depsc');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison of computation speed between Chebyshev and Lagrange            % 
% interpolation for degrees 50 and 100                                      % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = linspace(-1,1,10000);
i = 50;
chebyshev_nodes50 = chebyshev_nodes(i+1);

% Chebyshev n = 50
gemiddelde_chebyshev = 0;
tic
for j = 1:10
    b = f(chebyshev_nodes50)';
    [T, dT] = chebyshev(chebyshev_nodes50, i);
    coeff_chebyshev = T\b;
    chebpolyval(coeff_chebyshev, z);
end
gemiddelde_chebyshev = toc/10


% Lagrange n = 50
gemiddelde_lagrange = 0;
tic
for j = 1:10
    evalueer_lagrange(chebyshev_nodes50,f,z);
end
gemiddelde_lagrange = toc/10


i = 100;
chebyshev_nodes100 = chebyshev_nodes(i+1);

% Chebyshev n = 100
tic
for j = 1:10
    b = f(chebyshev_nodes100)';
    [T, dT] = chebyshev(chebyshev_nodes100, i);
    coeff_chebyshev = T\b;
    chebpolyval(coeff_chebyshev, z);
end
gemiddelde_chebyshev = toc/10

% Lagrange n = 100
tic
for j = 1:10
    evalueer_lagrange(chebyshev_nodes100,f,z);
end
gemiddelde_lagrange = toc/10





