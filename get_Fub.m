function F_bil = get_Fub(x,y,z,c,x0,y0,z0)


F_bil = 0.25*abs(c)*(x-y)^2 + 0.25*abs(c)*(x0+y0)^2 - 0.5*(x0+y0)*(x+y);

end % EOF