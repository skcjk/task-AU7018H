clear all
clc
close all
simulationiter1 =100;
simulationiter2 =200;
simulationiter3 = 300; 
DMC_sp1 = 5; %设定值
DMC_sp2 = 10; %设定值
DMC_sp3 = 15; %设定值
for i=1:simulationiter1
    y_set(1,i)=DMC_sp1;
end%init setpoint of y
for i=(simulationiter1+1):simulationiter2
    y_set(1,i)=DMC_sp2;
end%init setpoint of y
for i=(simulationiter2+1):simulationiter3
    y_set(1,i)=DMC_sp3;
end%init setpoint of y
 
% 实际系统
a_real=[1 -1.5 0.5];%A(q-1) of process 0.2
b_real=[1 0.6];%B(q-1) of process 
k_real=1;%time delay of process
%实际系统包含噪声
for i=1:simulationiter3
    u(1,i)=0;%init U
    noise(1,i)=0.5*randn()/10;%init random noise
end

for i=2:simulationiter3
    y_real(1,i)=noise(1,i);%init real y
end%init sum of noise

[a_model,b_model]  = RLS(y_real,u);
for t=3:simulationiter3     
k_model = k_real;
na=length(a_model)-1;
nb=length(b_model)-1;

%GPC配置参数%(实控的时候 需要)
GPC_p=3;%predict horizon
GPC_m=3;%control horizon
GPC_lambda=0.1;%control weight
GPC_alfa=0;%soften parameter 柔化因子
GPC_beta=0;%step scale 阶梯因子
 
matrix_e=zeros(k_model+GPC_p-1,1); % E 丢番图方程   （6,1） 生成零矩阵
matrix_f=zeros(na+1,k_model+GPC_p-1); %F 丢番图方程   （2.6）
matrixg_whole=zeros(k_model+GPC_p-1,nb+k_model+GPC_p-1); % G= E*B 丢番图方程
%（6,6）
%计算F 丢番图方程
matrix_f(1,1)=1-a_model(1,2);
for i=1:1:na-1
    matrix_f(i+1,1)=a_model(1,i+1)-a_model(1,i+2);
end
matrix_f(na+1,1)=a_model(1,na+1);

%计算E 丢番图方程
matrix_e(1,1)=1;
for j=2:1:k_model+GPC_p-1
    matrix_e(j,1)=matrix_f(1,j-1);
    matrix_f(1,j)=matrix_f(2,j-1)-matrix_e(j,1)*(a_model(1,2)-1);
    for i=1:1:na-1
        matrix_f(i+1,j)=matrix_f(i+2,j-1)-matrix_e(j,1)*(a_model(1,i+2)-a_model(1,i+1));
    end
    matrix_f(na+1,j)=matrix_e(j,1)*a_model(1,na+1);
end%init e,f
 
for i=1:nb+1
    matrixg_whole(1,i)=b_model(1,i);
end
for j=2:k_model+GPC_p-1
    for i=1:nb+j-1+1
        if i<=j-1
            matrixg_whole(j,i)=matrixg_whole(j-1,i);
        elseif i<=nb+j-1
            matrixg_whole(j,i)=matrixg_whole(j-1,i)+matrix_e(j,1)*b_model(1,i-j+1);
        elseif i==nb+j
            matrixg_whole(j,i)=matrix_e(j,1)*b_model(1,nb+1);
        end
    end
end%inint g
 
for i=1:GPC_p
    g_single(1,i)=matrixg_whole(k_model+GPC_p-1,i);
end
for j=1:GPC_p
    for i=1:j
        g_part(j,i)=g_single(1,j-i+1);
    end
end   
for j=1:GPC_p
    for i=1:GPC_m
        g_part2(j,i)=g_part(j,i);
    end
end%init g_part for teh use of control of different forms

% 离线计算D矩阵
temp = inv((g_part2)'*(g_part2)+GPC_lambda*eye(GPC_m))*(g_part2)';%inv求逆矩阵 
matrix_d = temp(1,:);  %取temp第一行
 
    y_real(1,t)=noise(1,t);

    for i=1:nb+1
        y_real(1,t)=y_real(1,t)+b_real(1,i)*u(1,t-k_real-i+1);
    end
    for i=1:na
        y_real(1,t)=y_real(1,t)-a_real(1,i+1)*y_real(1,t-i);
    end%sampleing
    
    for i=0:GPC_p-1
        y1(i+1,1)=0;
        for l=1:na+1
            y1(i+1,1)=y1(i+1,1)+matrix_f(l,i+k_model)*y_real(1,t-l+1);
        end
        for l=i+2:nb+k_model+i
            y1(i+1,1)=y1(i+1,1)+matrixg_whole(k_model+i,l)*(u(1,t+i+1-l)-u(1,t+i-l));
        end
    end%init y1
 
    w_start=0;
    if k_model==1
        w_start=y_real(1,t);
    else
        for l=1:na+1
            w_start= w_start+matrix_f(l,k_model-1)*y_real(1,t-l+1);
        end
        for l=1:nb+k_model-1
            w_start= w_start+matrixg_whole(k_model-1,l)*(u(1,t-l)-u(1,t-l-1));
        end%init w_start
    end
    w(1,1)=GPC_alfa*w_start+(1-GPC_alfa)*y_set(1,t);
    for i=1:GPC_p-1
        w(i+1,1)=GPC_alfa*w(i,1)+(1-GPC_alfa)*y_set(1,t);
    end%init w

    u_delta = matrix_d *(w-y1);
    u(1,t)=u(1,t-1)+u_delta;%get u if you choose normal gpc 
    
end

subplot(2,1,1),plot(y_set,'r');hold on;plot(y_real,'b');legend('设定值','输出值');hold on;axis([0,simulationiter3-10,-inf,inf])
subplot(2,1,2),plot(u,'g');legend('控制律');hold on;axis([0,simulationiter3-10,-inf,inf])

function [A, B] = RLS(y, u)
a = [1   -1.5  0.5]';b = [0 1   0.6]';d = 3;  
na = length(a)-1;nb = length(b)-1;  
L = 1000;  
uk = zeros(d+nb,1);yk = zeros(na,1);  
u = randn(L,1);  

xi = sqrt(0.1)*randn(L,1); 

w = 4;
theta = [a(2:na+1);b];
thetae_1 = zeros(na+nb+1,1);
P = 10^w*eye(na+nb+1);
for k = 1:L
    phi = [-yk;uk(d:d+nb)]; 
    y(k) = phi'*theta+xi(k);
   
    K = P*phi/(1+phi'*P*phi);
    thetae(:,k) = thetae_1+K*(y(k)-phi'*thetae_1);
    P = (eye(na+nb+1)-K*phi')*P;
   
    thetae_1 = thetae(:,k);
   
    for i = d+nb:-1:2
        uk(i) = uk(i-1);
    end
    uk(1) = u(k);
   
    for i = na:-1:2
        yk(i) = yk(i-1);
    end
    yk(1) = y(k);
end

A = [1;thetae_1(1:2)]';
B=thetae_1(4:5)';
end