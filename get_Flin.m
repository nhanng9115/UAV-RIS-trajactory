function F_lin = get_Flin(v,r,c,v0)
    F_qua = get_Fqua(v,r,v0);
    x0 = norm(v0-r)^2;
    F_lin = (c-1)*x0^c + c*x0^(c-1)*F_qua;
end % EOF