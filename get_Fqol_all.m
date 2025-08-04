function F_qol = get_Fqol_all(x,y,x0,y0)
% if length(size(x)) == 2
%     x = repmat(x,1,1,2);
%     x0 = repmat(x0,1,1,2);
% end
F_qol = (x0./y0).^2.*y - 2.*(x0./y0).*x;
end % EOF