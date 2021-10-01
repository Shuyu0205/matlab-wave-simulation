iy=117;
ix=133;
 
 
 fid=fopen('fort.122')
 X2=fscanf(fid,'%e',[ix,iy]);
 
 fid=fopen('fort.123')
 X3=fscanf(fid,'%e',[ix,iy]);
 
 fid=fopen('fort.125')
 X4=fscanf(fid,'%e',[ix,iy]);
  fid=fopen('fort.124')
 X4_1=fscanf(fid,'%e',[ix,iy]);
 fid=fopen('fort.127')
 X5=fscanf(fid,'%e',[ix,iy]);
 
 fid=fopen('fort.128')
 X6=fscanf(fid,'%e',[ix,iy]);
  
 fid=fopen('fort.129')
 X6_1=fscanf(fid,'%e',[ix,iy]);
 P_3=X3(:,60);
 
 P_4=X4(:,60);
 P_4_1=X4_1(:,60);
 P_5=X5(:,60);
 P_6=X6(:,60);
   P_6_1=X6_1(:,60);
subplot(2,1,1)
plot(P_4)
hold on
plot(P_4_1)
hold off
subplot(2,1,2)
 plot(P_6)
 hold on
 plot(P_6_1)
 
