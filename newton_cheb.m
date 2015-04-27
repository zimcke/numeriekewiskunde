function [ x ] = newton_cheb( c, x0, tol,nmax, afgeleide, h )
x=[1:nmax];
x(1)=x0;
if length(afgeleide)==5
    for k=2:nmax
        [f0,df0]=chebpolyval(c,x0);
        x1=x0-(f0/df0);
        x(k)=x1;
        if abs(x1-x0)<tol
            x(k:nmax)=x1;            
            break
        end
        x0=x1;
    end
elseif length(afgeleide)==12
    for k=2:nmax
        [f0,df0]=chebpolyval(c,x0);
        df0=afg_voorwaartse_diff(c,10,x0,h);
        x1=x0-(f0/df0);
        x(k)=x1;
        if abs(x1-x0)<tol
            x(k:nmax)=x1;
            break
        end
        x0=x1;
    end
end
end

