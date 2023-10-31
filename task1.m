%% 1)
clear ; clc;close all;
N = 1000;
n = 0:50;
 
u = -10 + (10+10)*rand(N,1);
e = randn(N,1);

figure 
plot(u,'k')
 
A = [1 -1];
B = [0 1.5];
m0 = idpoly(A,B);
u = iddata([],u);
e = iddata([],e);
y = sim(m0,[u e]);
z = [y,u];
m = arx(z,[1 1 1]);
 
A_est = m.A;
B_est = m.B;
m_est = idpoly(A_est,B_est);
y_est = sim(m_est,u);
 
figure 
plot(y,'r')
 
h = impz(B,A,n);
h_est = impz(B_est,A_est,n);
%% 2 a)
clear ; clc;close all;
N = 100;
n = 0:50;
 
u = zeros(N, 1);
u(1)=1;
e = zeros(N, 1);

figure 
plot(u,'k')
 
A = [1   -1.5347  0.7654];
B = [0   0.0654   -0.2060];
m0 = idpoly(A,B);
u = iddata([],u);
e = iddata([],e);
y = sim(m0,[u e]);
 
figure 
plot(y,'r')
%% 2 b)
clear ; clc;close all;
a = [1   -1.5347  0.7654]';b = [0.0654   -0.2060]';d = 3;  
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


