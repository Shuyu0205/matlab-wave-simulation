k_x=0:0.01:3.5;
a=[1/12,-8/12,0,8/12,-1/12];
k_1_x=[];
for w=k_x
    for j=-2:1:2
        temp_1=a(j+3)*exp(i*j*w);
    end
    temp=-i*temp_1;
    k_1_x=[k_1_x,temp];
end
plot(k_1_x/k_x)