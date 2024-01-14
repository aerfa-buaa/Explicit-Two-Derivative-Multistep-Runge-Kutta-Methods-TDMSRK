function [r,g]=Objective_fun(coeffs)
% Used by opt_TDMSRK for Objective Function
% Same as Ketchesons original code
r=coeffs(end);
g=zeros(size(coeffs));
g(end)=1;
end

 