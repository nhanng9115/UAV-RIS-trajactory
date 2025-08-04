function F_bil = get_FbilLB(x,y,c,x0,y0)

if c > 0
    F_bil = -0.25*(x-y)^2 - 0.25*(x0+y0)^2 + 0.5*(x0+y0)*(x+y);
else
    F_bil = -0.25*(x+y)^2 - 0.25*(x0-y0)^2 + 0.5*(x0-y0)*(x-y);
end

end