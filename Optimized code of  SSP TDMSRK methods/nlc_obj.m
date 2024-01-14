function [con,coneq]=nlc_obj(x,step,stage,order,K)
%=====================================================
% Extract arrays A,Ahat,v,vhat,d,b,theta from x
%Unpacking decision variables from x
 r=-x(end);
 [A,Ahat,v,vhat,d,b] =  unpackMSMDRK_all(x,step,stage,order);
[Re,P,Q] = matrix_vector_form(A,Ahat,v,vhat,d,b,r,K);
xx=[Re P Q ];  
con=-xx(:); %require entries in xx and r to be positive
coneq = Order_MSTDRK(A,Ahat,v,vhat,d,b,step,stage,order);

%This follows our SSP conditions defined in our paper. 
%Fmincon Syntax requires -xx<=0 to guarantee positivity of our coeficients
end % of function


