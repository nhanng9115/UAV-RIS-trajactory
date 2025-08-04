function F_pow = get_F0pow_all(x,c,x0)

F_pow = (c-1).*x0.^c - c.*x0.^(c-1).*x;

end % EOF