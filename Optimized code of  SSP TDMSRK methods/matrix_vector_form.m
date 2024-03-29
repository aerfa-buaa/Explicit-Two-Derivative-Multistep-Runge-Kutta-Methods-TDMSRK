function [Re,P,Q] = matrix_vector_form(A,Ahat,v,vhat,d,b,r,K)

%Y= Un+dt*A*F(Un)+dt^2*Ahat*Fdot(Un)
%Y= R*e*Un+ P*(Un+(dt/r)*F(Un))+Q*(Un+(dt^2/r2)*Fdot(Un))     

    r2=r^2/K^2; s=length(A); z=zeros(s+1,1); I=eye(s+1);e=ones(s+1,1);
    
	T=[[A;v],z];That=[[Ahat;vhat],z]; S=[d;b];         
    %Converting butcher to Modified Shu Osher
    R=inv((I+r*T+r2*That));  
    Re=R*S;
    P=R*r*T;
    Q=R*r2*That;
    
end
