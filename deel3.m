%% Deel 3: Methode van Newton-Raphson

%variabelen toekennen
format long
load('cheb_coeffs.mat')
nbpoints = 10000;
tol = 10^-50;
nmax = 100;
h = [10^-1 10^-1.5 10^-2 10^-3 10^-6];
x_start = linspace(-1,1,100);
x = linspace(-1,1,nbpoints);
[y,dy] = chebpolyval(c,x);

%nulpunten benaderen
clearvars nulp
N_R = newton_cheb(c,x_start(1),tol,nmax,'exact',h);
nulp = N_R(end);
for k = 2:100
    N_R = newton_cheb(c,x_start(k),tol,nmax,'exact',h);
    nulp(end+1)=N_R(end);
    
end

%nulpunten plotten op de veelterm
y1 = zeros(length(nulp));
figure(1)
plot(x,y,nulp,y1,'o')
title('Chebyshevveelterm met nulpunten')
legend('Chebyshevveelterm','nulpunten met Newton-Raphson benadering')
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('Chebyshevveelterm') % y-axis label

%plot: residu's i.f.v. de startwaarden
figure(2)
plot(x_start,nulp,'.')
title('residus i.f.v. de startwaarden')
xlabel('startwaarden') % x-axis label
ylabel('nulpunten met Newton-Raphson benadering') % y-axis label

%plot: raaklijn aan het punt waar de veelterm zijn absolute maximum bereikt
x_max = x(find(y==max(y)));
[y_max,dy_max] = chebpolyval(c,x_max);
raaklijn = y_max+dy_max.*(x-x_max);
figure(3)
plot(x,y,x,raaklijn)
title('raaklijn aan het absolute maximum')
legend('Chebyshevveelterm','raaklijn')
xlabel('x-waarden in [-1,1]') % x-axis label
ylabel('Chebyshevveelterm') % y-axis label


%plot de fout i.f.v. de iteratiestap voor startwaarde x0=-0.9, en x0=0.5
x_0 = [-0.9,0.5];
x_nul = zeros(length(x_0),nmax);
nulp_ben = zeros(length(x_0),nmax);
for k = 1:length(x_0)
    x_nul(k,:) = newton_cheb(c,x_0(k), tol, nmax, 'exact', h);
    nulp_ben(k,:) = chebpolyval(c,x_nul(k,:));
end
x_stap = (1:nmax);
figure(4)
semilogy(x_stap,abs(nulp_ben))
title('fout i.f.v de iteratiestap')
legend('x0=-0.9','x0=0.5')
xlabel('iteratiestap') % x-axis label
ylabel('fout') % y-axis label

%plot: fout i.f.v. de iteratiestap voor startwaarde x0=0.5 en h
nulp_exact = nulp_ben(2,:);
x_diff = zeros(length(h),nmax);
nulp_diff = zeros(length(h),nmax);
for k = 1:length(h)
    x_diff(k,:) = newton_cheb(c,0.5, tol, nmax, 'differenties', h(k));
    nulp_diff(k,:) = chebpolyval(c,x_diff(k,:));
end

figure(5)
semilogy(x_stap,abs(nulp_exact),x_stap,abs(nulp_diff))
title('fout i.f.v de iteratiestap')
legend('exact','h=10^-1','h=10^-1.5','h=10^-2','h=10^-3','h=10^-6')
xlabel('iteratiestap') % x-axis label
ylabel('fout') % y-axis label
    
