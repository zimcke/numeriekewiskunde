%% 3.1 Equidistante punten en het Runge fenomeen

clc
f = @fx;

nbpoints = 10000;
z = linspace(-1,1,nbpoints);
y = fx(z);

% Lagrange interpolation degree 8
x = linspace(-1,1,9);
polynom8 = evalueer_lagrange(x, f, z);

% Lagrange interpolation degree 10
x = linspace(-1,1,11);
polynom10 = evalueer_lagrange(x, f, z);

% Lagrange interpolation degree 12
x = linspace(-1,1,13);
polynom12 = evalueer_lagrange(x, f, z);

% Lagrange interpolation degree 14
x = linspace(-1,1,15);
polynom14 = evalueer_lagrange(x, f, z);

% Figure 1: Lagrange interpolation
figure(1),clf
plot(z,y)
hold on
plot(z, polynom8)
plot(z, polynom10)
plot(z, polynom12)
plot(z, polynom14)
legend('f(x)','L8', 'L10', 'L12', 'L14')
xlabel('x-values in [-1,1]') % x-axis label
ylabel('Lagrange interpolation') % y-axis label
print('report\3_1_1','-deps');

% Lagrange interpolation degree 20, equidistantial nodes 
x = linspace(-1,1,21);
polynom20 = evalueer_lagrange(x, f, z);
abs_error20 = abs(polynom20 - y);
rel_error20 = abs_error20 ./ y;

% Lagrange interpolation degree 50, equidistantial nodes 
x = linspace(-1,1,51);
polynom50 = evalueer_lagrange(x, f, z);
abs_error50 = abs(polynom50 - y);
rel_error50 = abs_error50 ./ y;

% Lagrange interpolation degree 20, chebychev nodes 
xc = 0;
for i = 1:21
   xc(i) = cos(((2*i - 1)*pi)/ (2*21));
end
polynom20c = evalueer_lagrange(xc, f, z);
abs_error20c = abs(polynom20c - y);
rel_error20c = abs_error20c ./ y;

% Lagrange interpolation degree 50, chebychev nodes 
for i = 1:51
   xc(i) = cos(((2*i - 1)*pi)/ (2*51));
end
polynom50c = evalueer_lagrange(xc, f, z);
abs_error50c = abs(polynom50c - y);
rel_error50c = abs_error50c ./ y;

% Figure 2: relative error with Lagrange interpolation (equidistant and
% chebychev nodes)
figure(2),clf
semilogy(z, rel_error20)
hold on
semilogy(z, rel_error50)
semilogy(z, rel_error20c)
semilogy(z, rel_error50c)
legend('Equidistant nodes n=20', 'Equidistant nodes n=50', 'Chebychev nodes n=20', 'Chebychev nodes n=50', 'Location','southwest')
xlabel('x-values in [-1,1]') % x-axis label
ylabel('Relative error') % y-axis label
print('report\3_1_2','-deps');

benaderingsfout_equidistant = zeros(1,50);
for i = 1:50
    x = linspace(-1,1,i+1);
    z = linspace(-1,1,10000);
    benaderingsfout_equidistant(i) = max(norm(fx(z) - evalueer_lagrange(x,f,z)));
end

benaderingsfout_chebychev = zeros(1,50);
for i = 1:50
    xc = cos( (2*(1:(i+1))-1)*pi/2/(i+1));
    benaderingsfout_chebychev(i) = max(norm(fx(z) - evalueer_lagrange(xc,f,z)));
end

figure(3),clf
semilogy(1:50, benaderingsfout_equidistant);
hold on
semilogy(1:50, benaderingsfout_chebychev);
legend('Benaderingsfout equidistant nodes', 'Benaderingsfout chebychev nodes')
xlabel('x-values in [-1,1]') % x-axis label
ylabel('Benaderingsfout') % y-axis label
print('report\3_1_3','-deps')


%% 3.2 Verschillende basissen

clc

f = @fx;
nbpoints = 100;
z = linspace(-1,1,nbpoints);

E_lagrange = zeros(1,80);
E_monomial = zeros(1,80);
E_chebyshev = zeros(1,80);
for i = 1:80
    warning('Off');
    chebnodes = cos( (2*(1:(i+1))-1)*pi/2/(i+1));
    M = monomiaal(chebnodes, i);
    b = f(chebnodes)';
    coeff_monomial = M\b;
    [T, dT] = chebyshev(chebnodes, i);
    coeff_chebyshev = T\b;
    warning('On');
    E_lagrange(i) = max(abs(fx(z) - evalueer_lagrange(chebnodes,f,z)));
    E_monomial(i) = max(abs(fx(z) - polyval(coeff_monomial,z)));
    E_chebyshev(i) = max(abs(fx(z) - chebpolyval(coeff_chebyshev, z)));
end


figure(5), clf
semilogy(1:80, E_lagrange);
hold on
semilogy(1:80, E_monomial);
semilogy(1:80, E_chebyshev);
legend('E Lagrange','E monomiaal', 'E Chebyshev', 'Location','southwest');
print('report\verschillende_basissen','-deps');

z = linspace(-1,1,10000);
gemiddelde_chebyshev = 0;
tic
for j = 1:10
    i = 50;
    chebnodes = cos( (2*(1:(i+1))-1)*pi/2/(i+1));
    b = f(chebnodes)';
    [T, dT] = chebyshev(chebnodes, i);
    coeff_chebyshev = T\b;
    chebpolyval(coeff_chebyshev, z);
end
gemiddelde_chebyshev = toc/10


gemiddelde_lagrange = 0;
tic
for j = 1:10
    i = 50;
    chebnodes = cos( (2*(1:(i+1))-1)*pi/2/(i+1));
    evalueer_lagrange(chebnodes,f,z);
end
gemiddelde_lagrange = toc/10

tic
for j = 1:10
    i = 100;
    chebnodes = cos( (2*(1:(i+1))-1)*pi/2/(i+1));
    b = f(chebnodes)';
    [T, dT] = chebyshev(chebnodes, i);
    coeff_chebyshev = T\b;
    chebpolyval(coeff_chebyshev, z);
end
gemiddelde_chebyshev = toc/10

tic
for j = 1:10
    i = 100;
    chebnodes = cos( (2*(1:(i+1))-1)*pi/2/(i+1));
    evalueer_lagrange(chebnodes,f,z);
end
gemiddelde_lagrange = toc/10

%  chebnodes = chebnodes';
% [T,dT] = chebychev(chebnodes,80);
%  y = fx(chebnodes);
%  x = y\T;
%  
%  figure(6), clf
%  plot(chebnodes,x);
% T = T';
%benaderingsfout_chebychev = zeros(80);
%for i = 1:80
%    benaderingsfout_chebychev(i) = max(norm(y - T(i+1)));
%end

% lagrange80 = evalueer_lagrange(xc, f, z);
% p = [25 0 1];
% monomial80 = 1./polyval(p, xc);
% c = [27/2 0 25/2];
% cheb80 = 1./chebpolyval(c, xc);

% figure(4),clf
% plot(z,y);
% hold on
% plot(z, lagrange80);
% plot(xc, monomial80);
% plot(xc, cheb80);
% legend('f(x)','L80', 'M80', 'Chebychev 80');
% xlabel('x-values in [-1,1]') % x-axis label
% ylabel('Polynomial interpolation') % y-axis label



