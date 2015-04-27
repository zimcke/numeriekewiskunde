function [ df0 ] = afg_voorwaartse_diff( c, k, x0 ,h )
df0=0;
for i=0:k
    if i==0;
        t=0;
    else
        t=voorwaartse_diff( c, i, x0, h);
    end
    df0=df0-((-1)^i)*t/h;
    
end

