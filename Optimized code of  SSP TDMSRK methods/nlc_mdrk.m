function [con,coneq]=nlc_mdrk(x,step,stage,order,K)
% The Nonlinear Constraints for our optimization routine are a combination
% of the equality constraints, coneq, which come from the order conditions,
% and Inequality Constraints, con, which are from our SSP conditions. 
% IE Modified Shu Osher Decomposition having positive Coeficients

%=====================================================
% Extract arrays A,d,b,theta from x
%Unpacking decision variables from x
% KK=0.6888921;
 r=-x(end);
%r=-x(end-2:end);
 [A,Ahat,v,vhat,d,b] =  unpackMSMDRK_all(x,step,stage,order);
% [A,Ahat,v,vhat,d,b] =  unpackMSMDRK(x,s);
[Re,P,Q] = matrix_vector_form(A,Ahat,v,vhat,d,b,r,K);
% xx=[Re P Q];  
xx=[Re P Q ];  
con=-xx(:); %require entries in xx and r to be positive
 coneq = Order_MSTDRK(A,Ahat,v,vhat,d,b,step,stage,order);
 % coneq=[ ];


%This follows our SSP conditions defined in our paper. 
%Fmincon Syntax requires -xx<=0 to guarantee positivity of our coeficients

end % of function


