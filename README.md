# Explicit-Two-Derivative-Multistep-Runge-Kutta-Methods-TDMSRK

## Optimized code of SSP TDMSRK methods
 
This folder is a collection of codes that we used to find Optimal SSP two-derivative multistep Runge-Kutta methods 

    # Main.m is the main program code.
`stage=2;          %Number of Stages `  
`step=2;            %Number of step`  
`order=3;           %Number of order`  
`K=1/sqrt(2);   %Second Derivative Coefficient `  
`%The optimization call:`  
`[X,FVAL,info]=fmincon(@Objective_fun,x0,[],[],[],[],lb,ub,@(x) opt_TDMSRK(x,step,stage,order,K),opts);`  

    # opt_TDMSRK.m is the main optimization process
    # Used by opt_TDMSRK for Objective Function
    # unpackMSMDRK_all.m is used to generate the coefficient variable of the TDMSRK methods.
    # matrix_vector_form.m is to convert the coefficient variable into matrices and vectors of the TDMSRK methods.
    # Order_MSTDRK.m is the order condition for the TDMSRK methods.
 


These Matlab codes are based on David Ketchesons and Zachary Grant's routines.  
This code evolved from https://github.com/ketch/RK-opt and https://github.com/SSPmethods/SSPMultiStageTwoDerivativeMethods.

## Data of optimal SSP TDMSRK

This folder gives the optimal coefficient matrices of SSP TDMSRK methods.

The TDMSRK_sqp method denotes the p-order TDMSRK scheme with q-step and s-stage.

## Numerical experiments



