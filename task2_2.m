%% 1
num=1; 
den=[1 2 1]; 
t0=0.5; %采样周期
gs=tf(num,den); %连续时间过程传递函数
gd=c2d(gs,t0,'zoh') %零阶保持器离散传递函数

clear 
num=1; 
den=[1 1.5 0]; 
t0=0.5; %采样周期
gs=tf(num,den); %连续时间过程传递函数
gd=c2d(gs,t0,'zoh') %零阶保持器离散传递函数
%% 2
clear 
%极点配置间接自校正控制
clear all; close all;

a=[1 -1.4724 0.4724]; b=[0.09883 0.077];  c=[1 0.5]; d=5; %对象参数
Am=[1 -1.213 0.3679]; %期望闭环特征多项式
na=length(a)-1; nb=length(b)-1; nc=length(c)-1; %na、nb、nc为多项式A、B、C
nam=length(Am)-1; %Am阶次
nf=nb+d-1; ng=na-1;

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
yk=zeros(na,1); %输出初值
yrk=zeros(na,1); %期望输出初值
xik=zeros(nc,1); %白噪声初值
xiek=zeros(nc,1); %白噪声估计初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %期望输出
xi=sqrt(0.01)*randn(L,1); %白噪声序列

%RELS初值
thetae_1=0.001*ones(na+nb+1+nc,1);
P=10^6*eye(na+nb+1+nc);
lambda=1; %遗忘因子[0.9 1]
for k=1:L
    time(k)=k;
    y(k)=-a(2:na+1)*yk+b*uk(d:d+nb)+c*[xi(k);xik]; %采集输出数据
    
    %递推增广最小二乘法
    phie=[-yk(1:na);uk(d:d+nb);xiek];
    K=P*phie/(lambda+phie'*P*phie);
    thetae(:,k)=thetae_1+K*(y(k)-phie'*thetae_1);
    P=(eye(na+nb+1+nc)-K*phie')*P/lambda;
    
    xie=y(k)-phie'*thetae(:,k); %白噪声的估计值
    
    %提取辨识参数
    ae=[1 thetae(1:na,k)']; be=thetae(na+1:na+nb+1,k)'; ce=[1 thetae(na+nb+2:na+nb+1+nc,k)'];
    if nc>0
        if abs(ce(2))>0.8
            ce(2)=sign(ce(2))*0.8;
        end
    end
     %多项式B的分解
    br=roots(be); %求B的根
    b0=be(1); b1=1; %b0为B-；b1为B+
    Val=0.0; %通过修改临界值，确定B零点是否对消（零点绝对值小于临界值，则被抵消）
    for i=1:nb %分解B-、B+
        if abs(br(i))>=Val
            b0=conv(b0,[1 -br(i)]);
        else
            b1=conv(b1,[1 -br(i)]);
        end
    end
    
    Bm1=sum(Am)/sum(b0); %确定多项式Bm'
    
    %确定多项式A0
    %A0=ce; %可取A0=C
    na0=2*na-1-nam-(length(b1)-1); %观测器最低阶次
    A0=1;
    for i=1:na0
        A0=conv(A0,[1 0.5]); %生成观测器
    end

    %计算Diophantine方程，得到F、G、R
    [F1,G]=diophantine(ae,b0,d,A0,Am);
    F=conv(F1,b1); R=Bm1*A0; 
    nr=length(R)-1;

    u(k)=(-F(2:nf+1)*uk(1:nf)+R*[yr(k+d:-1:k+d-min(d,nr));yrk(1:nr-d)]-G*[y(k);yk(1:ng)])/F(1);%求控制量
    
    %更新数据
    thetae_1=thetae(:,k);
    
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=na:-1:2
        yk(i)=yk(i-1);
        yrk(i)=yrk(i-1);
    end
    yk(1)=y(k);
    yrk(1)=yr(k);
    for i=nc:-1:2
        xik(i)=xik(i-1);
        xiek(i)=xiek(i-1);
    end
    if nc>0
        xik(1)=xi(k);
        xiek(1)=xie;
    end
end
figure(1);
subplot(2,1,1);
plot(time,yr(1:L),'r:',time,y);
xlabel('k'); ylabel('y_r(k)、y(k)');
legend('y_r(k)','y(k)'); axis([0 L -20 20]);
subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)'); axis([0 L -40 40]);
%% 3
%极点配置直接自校正控制 
clear all; 
close all;

a=[1 -1.4724 0.4724]; b=[0.09883 0.077]; d=3; Am=[1 -1.213 0.3679]; %对象参数及期望闭环特征多项式
na=length(a)-1; nb=length(b)-1; nam=length(Am)-1; %na、nb、nam为多项式A、B、Am阶次
nf=nb+d-1; ng=na-1; %F、G的阶次

%确定多项式A0
na0=2*na-nam-nb-1; %观测器最低阶次
A0=1;
for i=1:na0
    A0=conv(A0,[1 0.3-i*0.1]);%生成观测器
end
AA=conv(A0,Am); naa=na0+nam; %A0*Am
nfg=max(naa,max(nf,ng)); %用于ufk、yuf更新
nr=na0; %R的阶次

L=400; %控制步数
uk=zeros(d+nb,1); %输入初值：uk(i)表示u(k-i)
ufk=zeros(d+nfg,1); %滤波输入初值
yk=zeros(max(na,d),1); %输出初值
yfk=zeros(d+nfg,1); %滤波输出初值
yrk=zeros(max(na,d),1); %期望输出初值
yr=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %期望输出
%RLS初值
thetae_1=0.001*ones(nf+ng+2,1);
P=10^6*eye(nf+ng+2);
lambda=1; %遗忘因子[0.9 1]
for k=1:L
    time(k)=k;
    y(k)=-a(2:na+1)*yk(1:na)+b*uk(d:d+nb); %采集输出数据
    ufk(d)=-AA(2:naa+1)*ufk(d+1:d+naa)+uk(d); %滤波输入输出
    yfk(d)=-AA(2:naa+1)*yfk(d+1:d+naa)+yk(d);
    
    %递推最小二乘法
    phie=[ufk(d:d+nf);yfk(d:d+ng)];
    K=P*phie/(lambda+phie'*P*phie);
    thetae(:,k)=thetae_1+K*(y(k)-phie'*thetae_1);
    P=(eye(nf+ng+2)-K*phie')*P/lambda;
    
    %提取辨识参数
    be0=thetae(1,k); thetaeb(:,k)=thetae(:,k)/be0;
    Fe=thetaeb(1:nf+1,k)'; Ge=thetaeb(nf+2:nf+ng+2,k)';
        
    Bm1=sum(Am)/be0; %Bm'
    R=Bm1*A0;
    
    u(k)=(-Fe(2:nf+1)*uk(1:nf)+R*[yr(k+d:-1:k+d-min(d,nr));yrk(1:nr-d)]-Ge*[y(k);yk(1:ng)])/Fe(1); %控制量
     %更新数据
    thetae_1=thetae(:,k);
    
    for i=d+nb:-1:2
        uk(i)=uk(i-1);
    end
    uk(1)=u(k);
    
    for i=max(na,d):-1:2
        yk(i)=yk(i-1);
        yrk(i)=yrk(i-1);
    end
    yk(1)=y(k);
    yrk(1)=yr(k);
    
    for i=d+nfg:-1:d+1
        ufk(i)=ufk(i-1);
        yfk(i)=yfk(i-1);
    end
end

figure(1);
subplot(2,1,1);
plot(time,yr(1:L),'r:',time,y);
xlabel('k'); ylabel('y_r(k)、y(k)');
legend('y_r(k)','y(k)'); axis([0 L -20 20]);

subplot(2,1,2);
plot(time,u);
xlabel('k'); ylabel('u(k)'); axis([0 L -5 5]);
    
function [F1,G]=diophantine(A,B,d,A0,Am)
%***********************************************************************
  %功能：Diophanine方程的求解
  %调用格式：[F1,G]=diophantine(A,B,d,A0,Am)
  %输入参数：多项式A、B系数向量、纯延迟d、多项式A0、Am系数向量（行向量）
  %输出参数：Diophanine方程的解F1、G（行向量）
%***********************************************************************  
dB=[zeros(1,d) B];
na=length(A)-1; nd=length(dB)-1;
T1=conv(A0,Am); nt=length(T1); T=[T1';zeros(na+nd-nt,1)];

%得到Sylvester 矩阵
AB=zeros(na+nd);
for i=1:na+1
    for j=1:nd
        AB(i+j-1,j)=A(i);
    end
end
for i=1:nd+1
    for j=1:na
        AB(i+j-1,j+nd)=dB(i);
    end
end
%得到F1,G
L=(AB)\T;
F1=[ L(1:nd)]';
G=[ L(nd+1:na+nd)]';
end