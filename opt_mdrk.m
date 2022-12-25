%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimization Driver File for Finding optimal SSP Two Derivative multistep Runge Kutta Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
restart=0;   % Start with randomly generated coeficients
stage=2;          %Number of Stages
step=2;            %Number of step
order=3;           %Number of order
K=1/sqrt(2);   %Second Derivative Coefficient (dtVV/dtFE)

minreff =1.e-4; %Keep looking until method with at least this value is found
if restart==0
    %Because its generated from random starting point
    %we are being less restrictive on satisfying constraints.
%     clear A Ahat b bhat x X r   %Be sure all variables are reset
    X1=[0];
    opts=optimset('MaxFunEvals',10000000,'TolCon',1.e-15,'TolFun',1.e-15,'TolX',1.e-15,...
        'GradObj','on','MaxIter',10000000,'Diagnostics','on','Display','on',...
      'UseParallel','never','Algorithm','sqp'); % 'UseParallel','never','Algorithm','active-set');  %    
  
else
    X1=X;% stores original coefficients from loaded method
    opts=optimset('MaxIter',1000000,'MaxFunEvals',1000000,...
        'TolCon',1.e-15,'TolFun',1.e-15,'TolX',1.e-15,...
        'GradObj','on','Diagnostics','off','Display','off',...
        'Algorithm','sqp','UseParallel','never');
end
%%
n=stage*step+(2*step+stage-2)*(stage-1)+2*(step+stage-1)+1;
    lb=0+zeros(1,n);    lb(end)=-1.8; 
    ub=1+zeros(1,n);     ub(end)=-0.0501;        %requires r>=0


%==============================================
count=0;                                     %Count tracks the number of times optimizer has failed to find a method
info=-2;

%This While loop requires optimizer to keep running until r>minreff is found while satisfying all constraints
while (info==-2 || (r)<minreff || info==0)
    if count==1 %If fails to find a method after 100 times, stop routine
        ('exceed count')
        x=X1;
        r=-x(end);
        r=101;
        break
    end
    %defining initial starting point for Fmincon
    if restart==1
        x0=X1;
    elseif restart==2
        x0=X1+.4*rand(1,length(X1));
    else
        x0=[(2*rand(1,n-1)),-.01];
    end

    %==============================================
    %The optimization call:
    %  [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],Aeq,beq,lb,ub,@(x) nlc_mdrk(x,s,p,CC),opts);
    [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],[],[],lb,ub,@(x) nlc_mdrk(x,step,stage,order,K),opts);
    r=-FVAL;
    count=count+1;
end %while loop
%==============================================
  [A,Ahat,v,vhat,d,b] =  unpackMSMDRK_all(X,step,stage,order);
   coneq = Order_MSTDRK(A,Ahat,v,vhat,d,b,step,stage,order);
    r0=-X(end);
  [Re,P,Q] = Butcher2ShuOsher(A,Ahat,v,vhat,d,b,r0,K);


rr=X;


 