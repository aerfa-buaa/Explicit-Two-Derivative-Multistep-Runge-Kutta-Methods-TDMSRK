function [Re,P,Q] = Butcher2ShuOsher(A,Ahat,v,vhat,d,b,r,K)
%function [Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,K);
%Converting Butcher form to Modified Shu Osher
% w=S*x+dt*T*F(Un)+ dt^2*That*Fdot(Un)


%Y= Un+dt*A*F(Un)+dt^2*Ahat*Fdot(Un)
%Y= R*e*Un+ P*(Un+(dt/r)*F(Un))+Q*(Un+(dt^2/r2)*Fdot(Un))     

    r2=r^2/K^2; n=length(A); z=zeros(n+1,1); I=eye(n+1);e=ones(n+1,1);
    
	T=[d;b];   S=[[A;v],z];   Shat=[[Ahat;vhat],z];       
    %Converting butcher to Modified Shu Osher
    R=inv((I+r*S+r2*Shat));  
    Re=R*T;
    P=R*r*S;
    Q=R*r2*Shat;
    
end
