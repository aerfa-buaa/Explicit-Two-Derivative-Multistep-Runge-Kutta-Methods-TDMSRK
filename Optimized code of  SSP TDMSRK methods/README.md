    Optimized code of SSP TDMSRK methods:

    Main.m is the main program code.
%___________________________________________
stage=2;          %Number of Stages
step=2;            %Number of step
order=3;           %Number of order
K=1/sqrt(2);   %Second Derivative Coefficient 

%The optimization call:
[X,FVAL,info]=fmincon(@Objective_fun,x0,[],[],[],[],lb,ub,@(x) opt_TDMSRK(x,step,stage,order,K),opts);

%_______________________________________________
    opt_TDMSRK.m is the main optimization process
    Used by opt_TDMSRK for Objective Function
    unpackMSMDRK_all.m is used to generate the coefficient variable of the TDMSRK methods.
    matrix_vector_form.m is to convert the coefficient variable into matrices and vectors of the TDMSRK methods.
    Order_MSTDRK.m is the order conditions for the TDMSRK methods.
    
    


