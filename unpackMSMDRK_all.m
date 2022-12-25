%%%Unpacking Explicit 2 Derivative  multistep Runge Kutta Method
function [A,Ahat,v,vhat,d,b] =  unpackMSMDRK_all(x,step,stage,order)
count1= (stage-1)*step;  %d
count2= count1 + step;   %b
count3= count2 + (2*step+stage-2)*(stage-1)/2; %A
count4= count3 + (step+stage-1);    %v
count5 = count4 +(2*step+stage-2)*(stage-1)/2;  %Ahat 
count6= count5 + (step+stage-1);     %vhat 


di=eye(step);
d=[di;reshape(x(1: count1),stage-1,step)];
b= reshape(x(count1+1:count2),1,step) ;
A=zeros(step+stage-1, step+stage-1);
count=count2;
for i = step+1 : step+stage-1
A(i,1:i-1)=x(count+1:count+i-1);
count= count+i-1;
end

v=[ reshape(x(count3+1:count4),1,step+stage-1)];

Ahat=zeros(step+stage-1, step+stage-1);
count=count4;
for i = step+1 : step+stage-1
Ahat(i,1:i-1)=x(count+1:count+i-1);
count= count+i-1;
end

vhat=[ reshape(x(count5+1:count6),1,step+stage-1)];


end

