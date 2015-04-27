function [ delta_i ] = voorwaartse_diff( c, i, x0, h )

if i==0;
    delta_i=chebpolyval(c,x0);
else
    delta_i=voorwaartse_diff(c,i-1,x0+h,h)-voorwaartse_diff(c,i-1,x0,h);
end

