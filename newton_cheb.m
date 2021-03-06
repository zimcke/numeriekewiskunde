% Newton Raphson method for approximating the zero of a real-valued
% function
% input c = coefficients of Chebyshev polynomial
% input x0 = start value
% input tol = tolerance between two iteration nodes
% input nmax = maximum number of iteration steps
% input afgeleide = string with value 'exact' or 'differenties'
% input h
% output x = vector containing the successive approximations of the zero
function [ x ] = newton_cheb( c, x0, tol,nmax, afgeleide, h )
x = (1:nmax);
x(1)= x0;

if strcmp(afgeleide,'exact')
    for k = 2:nmax
        [f0,df0] = chebpolyval(c,x0);
        x1 = x0 -(f0/df0);
        x(k) = x1;
        if abs(x1-x0) < tol
            x(k:nmax) = x1;            
            break
        end
        x0 = x1;
    end
elseif strcmp(afgeleide,'differenties')
    for k = 2:nmax
        [f0,~] = chebpolyval(c,x0); % ~ means this output argument is ignored => we use voorwaartse diff instead
        df0 = afg_voorwaartse_diff(c,10,x0,h);
        x1 = x0 -(f0/df0);
        x(k) = x1;
        if abs(x1-x0) < tol
            x(k:nmax) = x1;
            break
        end
        x0 = x1;
    end
end

end

