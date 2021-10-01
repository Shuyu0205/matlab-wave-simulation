clear
clear
format long;
V=[1;0;0];%set v_0
t=0;%set t_0=0
h=0.005;%set time step length h
N=100/h;%get the step number
m=1/4-sqrt(3)/6;%for shorter jacobian
p=1/4+sqrt(3)/6;%for shorter jacobian
q=1/4;%for shorter jacobian
temp=[];
b_1=0.5; %set b_1
b_4=0.5; %set b_1
T=25;
c=1;
F_k=[];
h=0.016;
k=0.008;
y_0=0.5;
q=(c*k)/(h);%here we let q=0.5
N=(4/h)+1;
M=(4/h)+1;
J=(T/k)+1;
delta=0.2;%wave source scale
d=0.5;%the width of slit
%build up the empty solution matrix;
U= zeros(M,N,J);

for m=1:1:M %y
    for n=1:1:N%x
        r=sqrt((-2+((n-1)*h))^2+(((-2+(m-1)*h))-y_0)^2);
        if (r>=-delta) && (r<=delta)%locate the wave source
            U(m,n,1)=cos((pi*r)/(2*delta));
        end
    end
end
%Second page of the solution
for m=2:1:M-1%y
    for n=2:1:N-1%x
        U(m,n,2)=0.5*((q^2)*((U(m,n+1,1)-2*U(m,n,1)+U(m,n-1,1))+(U(m+1,n,1)-2*U(m,n,1)+U(m-1,n,1)))+2*U(m,n,1));%first time layer
    end
end
%Bottom solid wall 
for j=1:1:J
    for n=1:1:N
        U(1,n,j)=0;
    end
end
%Left solid wall
for j=1:1:J
    for m=1:1:M
        if (-2+(m-1)*h)<=2
            U(m,1,j)=0;
        end
    end
end
%right solid wall boundary
for j=1:1:J
    for m=1:1:M
        if (-2+(m-1)*h)<=2
            U(m,N,j)=0;
            m_location=m;%save it for upper solid wall
        end
    end
end
%Upper solid wall boundary
for j=1:1:J
    for n=1:1:N
        if (-2+(n-1)*h)<=0 |  (-2+(n-1)*h)>=0
            U(m_location,n,j)=0;
        end
    end
end
%main calculation
for j=3:1:J
    for m=2:1:M-1
        for n=2:1:N-1
            if m==m_location
                if (-2+(n-1)*h)<=d && (-2+(n-1)*h)>=-d
                    U(m,n,j)=(q^2)*((U(m,n+1,j-1)-2*U(m,n,j-1)+U(m,n-1,j-1))+(U(m+1,n,j-1)-2*U(m,n,j-1)+U(m-1,n,j-1)))+2*U(m,n,j-1)-U(m,n,j-2);
                end
            else
                U(m,n,j)=(q^2)*((U(m,n+1,j-1)-2*U(m,n,j-1)+U(m,n-1,j-1))+(U(m+1,n,j-1)-2*U(m,n,j-1)+U(m-1,n,j-1)))+2*U(m,n,j-1)-U(m,n,j-2);
            end                    
        end
    end
    %non-reflection boundary
    %left and right
    for m=m_location+1:1:M-1
        U(m,1,j)=(1/(1+(q/2)))*((q/2)*(U(m,2,j)-U(m,2,j-2)+U(m,1,j-2))+((q^2)/2)*(U(m+1,1,j-1)-2*U(m,1,j-1)+U(m-1,1,j-1))+2*U(m,1,j-1)-U(m,1,j-2));
        U(m,N,j)=(1/(1+(q/2)))*(-(q/2)*(-U(m,N-1,j)-U(m,N,j-2)+U(m,N-1,j-2))+((q^2)/2)*(U(m+1,N,j-1)-2*U(m,N,j-1)+U(m-1,N,j-1))+2*U(m,N,j-1)-U(m,N,j-2));
    end
    %Top
    for n=2:1:N-1
        U(M,n,j)=(1/(1+(q/2)))*(-(q/2)*(-U(M-1,n,j)-U(M,n,j-2)+U(M-1,n,j-2))+((q^2)/2)*(U(M,n+1,j-1)-2*U(M,n,j-1)+U(M,n-1,j-1))+2*U(M,n,j-1)-U(M,n,j-2));
    end
    %Two corner
    U(M,1,j)=(1/(-3-2*q))*(q*((-U(M-1,1,j)-U(M,1,j-2)+U(M-1,1,j-2))-(U(M,2,j)-U(M,2,j-2)+U(M,1,j-2)))-6*U(M,1,j-1)+3*U(M,1,j-2));%(-2,2) point
    U(M,N,j)=(1/(-3-2*q))*(q*((-U(M-1,N,j)-U(M,N,j-2)+U(M-1,N,j-2))+(-U(M,N-1,j)-U(M,N,j-2)+U(M,N-1,j-2)))-6*U(M,N,j-1)+3*U(M,N,j-2));%(2,2) point
end
for n=2:1:15
  e=1;%set a initial error, can be arbitrary, better big.
  k_np=[1;0;0;0;0;0];%set an initial guess. It's a safe guess for good reason.
  %Newton's iteration
  while e>1 %test the tolerance
      preserve=k_np; %keep the old value so that can calculte new tolerance
      %Jacobian
      J=[1+0.04*q*h,-(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*q*h,-(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,0.04*m*h,-(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*m*h,-(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h;
         -0.04*q*h,1+(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*q*h+3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,-0.04*m*h,(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*m*h+3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h,(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h;
         0,-3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,1,0,-3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h,0;
         0.04*p*h,-(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*p*h,-(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h,1+0.04*q*h,-(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*q*h,-(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h;
         -0.04*p*h,(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*p*h+3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h  ,(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h,-0.04*q*h,1+(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*q*h+3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h,(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h;
         0,-3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h,0,0,-3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h,1];
     %Target function F
     F_k=[k_np(1,:)+0.04*(V(1,n-1)+q*h*k_np(1,:)+m*h*k_np(4,:))-(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:));
           k_np(2,:)-0.04*(V(1,n-1)+q*h*k_np(1,:)+m*h*k_np(4,:))+(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))+3*(10^7)*((V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))^2);
           k_np(3,:)-3*(10^7)*((V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))^2);
           k_np(4,:)+0.04*(V(1,n-1)+p*h*k_np(1,:)+q*h*k_np(4,:))-(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:));
           k_np(5,:)-0.04*(V(1,n-1)+p*h*k_np(1,:)+q*h*k_np(4,:))+(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))+3*(10^7)*((V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))^2);
           k_np(6,:)-3*(10^7)*((V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))^2)];
      %Newton's method
      k_np=k_np-inv(J)*F_k;
      %New Tolerance
      e=norm(k_np-preserve)/norm(k_np);
  end
  %Update result
  k_1=k_np(1:3);%get k_1
  k_4=k_np(4:6);%get k_4
  %Perform the final step of the RK method to get an solution at one point.
  temp=V(:,n-1)+b_1*h.*k_1+b_4*h.*k_4;
  V=[V,temp];% update and n=n+1 then.
end

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

subplot(2,2,1)
surf(-2:h:2,-2:h:2,U(:,:,35));
shading flat
subplot(2,2,2)
surf(-2:h:2,-2:h:2,U(:,:,150));
shading flat
subplot(2,2,3)
surf(-2:h:2,-2:h:2,U(:,:,210));
shading flat
subplot(2,2,4)
surf(-2:h:2,-2:h:2,U(:,:,300));
shading flat
