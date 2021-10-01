clear
format long;
k=0.1%spatial
h=0.56%time
c=1;
a=[];
for i=0:1:(30/k)
    a=[a,0.4955*exp((-log(2))*(((-15+i*k)/3)^2))]
end
t=0;%set t_0=0

h=0.005;%set time step length h
%This function may take around 1 min to provide result. Be patient pls.
%Change h may results in the bad behavior of the solution.
temp=[];
u=[a;]
%4-Stage ERK, straightforward algorithom.
for n=1:1:200
   
    s=u(n,:);
    k_1=f(s);%calculate k1
    k_2=f(s+(h/2)*(k_1));
    k_3=f(s+(h/2)*k_2);
    k_4=f(s+h*k_3);
    temp=s+h*(k_1./6+k_2./3+k_3./3+k_4./6);
    u=[u;temp];
end
%{
for j=1:1:200
     plot(-15:k:15,u(j,:))
     pause(0.03)
end
%}
%Coefficients
a(1,1)=0.377847764031163;
a(2,1)=0.385232756462588;
a(2,2)=0.461548399939329;
a(3,1)=0.675724855841358;
a(3,2)=-0.061710969841169;
a(3,3)=0.241480233100410;
b(1)=0.750869573741408;
b(2)=-0.362218781852651;
b(3)=0.611349208111243;
c(1)=0.257820901066211;
c(2)=0.434296446908075;
c(3)=0.758519768667167;

%Define the function for the odes system.
function f_n=f(x)
 f_n=[];
 l=length(x);
 f_n=x; 
 for j=2:l-1
 if j==1
     f_n(j)=10*(-x(j+1)+x(j))/(0.1);

 elseif j==l
     f_n(j)=10*(-x(j)+x(j-1))/(0.1);
 else
        
        f_n(j)=10*(-x(j+1)+x(j-1))/(2*0.1);
 end
 end
end



