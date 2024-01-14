classdef Vehicle < handle
    %VEHICLE REMUS模型

    properties
        %% 本体参数
        xg; yg; zg;  % Center of Gravity wrt Origin at CB
        Ixx; Iyy; Izz;  % Moments of Inertia wrt Origin at CB
        rho; W; B; m;
        %% STD REMUS Non-Linear Maneuvering Coefficients: Forces
        Xuu; Xudot; Xwq; Xqq; Xvr; Xrr
        Yvv; Yrr; Yuv; Yvdot; Yrdot; Yur; Ywp; Ypq; Yuudr
        Zww; Zqq; Zuw; Zwdot; Zqdot; Zuq; Zvp; Zrp; Zuuds
        %% STD REMUS Non-Linear Maneuvering Coefficients: Moments
        Kpp; Kpdot;
        Mww; Mqq; Muw; Mwdot; Mqdot; Muq; Mvp; Mrp; Muuds
        Nvv; Nrr; Nuv; Nvdot; Nrdot; Nur; Nwp; Npq; Nuudr
        K;
        %% 状态量及初始值
        X=zeros(12, 1);  % X=[u v w p q r x y z phi theta psi]
    end

    methods
        function obj = Vehicle(p)
            %VEHICLE 构造此类的实例
            %   将profile中的参数复制到vehicle中
            %   copy data from profile to vehicle
            %% 本体参数
            obj.xg=p.xg; obj.yg=p.yg; obj.zg=p.zg;  % Center of Gravity wrt Origin at CB
            obj.Ixx=p.Ixx; obj.Iyy=p.Iyy; obj.Izz=p.Izz;  % Moments of Inertia wrt Origin at CB
            obj.rho=p.rho; obj.W=p.W; obj.B=p.B;
            %% STD REMUS Non-Linear Maneuvering Coefficients: Forces
            obj.Xuu=p.Xuu;  % Cross-flow Drag
            obj.Xudot=p.Xudot;  % Added Mass
            obj.Xwq=p.Xwq;  % Added Mass Cross-term
            obj.Xqq=p.Xqq;  % Added Mass Cross-term
            obj.Xvr=p.Xvr;  % Added Mass Cross-term
            obj.Xrr=p.Xrr;  % Added Mass Cross-term
            obj.Yvv=p.Yvv;  % Cross-flow Drag
            obj.Yrr=p.Yrr;  % Cross-flow Drag
            obj.Yuv=p.Yuv;  % Body Lift Force and Fin Lift
            obj.Yvdot=p.Yvdot;  % Added Mass
            obj.Yrdot=p.Yrdot;  % Added Mass
            obj.Yur=p.Yur;  % Added Mass Cross Term and Fin Lift
            obj.Ywp=p.Ywp;  % Added Mass Cross-term
            obj.Ypq=p.Ypq;  % Added Mass Cross-term
            obj.Yuudr=p.Yuudr;  % Fin Lift Force
            obj.Zww=p.Zww;  % Cross-flow Drag
            obj.Zqq=p.Zqq;  % Cross-flow Drag
            obj.Zuw=p.Zuw;  % Body Lift Force and Fin Lift
            obj.Zwdot=p.Zwdot;  % Added Mass
            obj.Zqdot=p.Zqdot;  % Added Mass
            obj.Zuq=p.Zuq;  % Added Mass Cross Term and Fin Lift
            obj.Zvp=p.Zvp;  % Added Mass Cross Term
            obj.Zrp=p.Zrp;  % Added Mass Cross Term
            obj.Zuuds=p.Zuuds;  % Fin Lift Force
            %% STD REMUS Non-Linear Maneuvering Coefficients: Moments
            obj.Kpp=p.Kpp;  % Rolling Resistance
            obj.Kpdot=p.Kpdot;  % Added Mass
            obj.Mww=p.Mww;  % Cross-flow Drag
            obj.Mqq=p.Mqq;  % Cross-flow Drag
            obj.Muw=p.Muw;  % Body and Fin Lift and Munk Moment
            obj.Mwdot=p.Mwdot;  % Added Mass
            obj.Mqdot=p.Mqdot;  % Added Mass
            obj.Muq=p.Muq;  % Added Mass Cross Term and Fin Lift
            obj.Mvp=p.Mvp;  % Added Mass Cross Term
            obj.Mrp=p.Mrp;  % Added Mass Cross Term
            obj.Muuds=p.Muuds;  % Fin Lift Moment
            obj.Nvv=p.Nvv;  % Cross-flow DragKprop
            obj.Nrr=p.Nrr;  % Cross-flow Drag
            obj.Nuv=p.Nuv;  % Body and Fin Lift and Munk Moment
            obj.Nvdot=p.Nvdot;  % Added Mass
            obj.Nrdot=p.Nrdot;  % Added Mass
            obj.Nur=p.Nur;  % Added Mass Cross Term and Fin Lift
            obj.Nwp=p.Nwp;  % Added Mass Cross Term
            obj.Npq=p.Npq;  % Added Mass Cross Term
            obj.Nuudr=p.Nuudr;  % Fin Lift Moment

            obj.m = obj.W / obj.rho;
            obj.K = [obj.m - obj.Xudot 0 0 0 obj.m * obj.zg -obj.m * obj.yg;
                   0 obj.m - obj.Yvdot 0 -obj.m * obj.zg 0 obj.m * obj.xg - obj.Yrdot;
                   0 0 obj.m - obj.Zwdot obj.m * obj.yg -obj.m * obj.xg - obj.Zqdot 0;
                   0 -obj.m * obj.zg obj.m * obj.yg obj.Ixx - obj.Kpdot 0 0;
                   obj.m * obj.zg 0 -obj.m * obj.xg - obj.Mwdot 0 obj.Iyy - obj.Mqdot 0;
                   -obj.m * obj.yg obj.m * obj.xg - obj.Nvdot 0 0 0 obj.Izz - obj.Nrdot];
        end

        function f = kinetic_equation(obj, U)
            %METHOD1 动力学与运动学
            %   kinetic_equation:convert X to Xdot
            dr=U(1);    ds=U(2);    Xprop=U(3); Kprop=U(4);
            J1=[cos(obj.X(12))*cos(obj.X(11))  -sin(obj.X(12))*cos(obj.X(10))+cos(obj.X(12))*sin(obj.X(11))*sin(obj.X(10))  sin(obj.X(12))*sin(obj.X(10))+cos(obj.X(12))*sin(obj.X(11))*cos(obj.X(10));
                sin(obj.X(12))*cos(obj.X(11))  cos(obj.X(12))*cos(obj.X(10))+sin(obj.X(12))*sin(obj.X(11))*sin(obj.X(10))  -cos(obj.X(12))*sin(obj.X(10))+sin(obj.X(12))*sin(obj.X(11))*cos(obj.X(10));
                -sin(obj.X(11))  cos(obj.X(11))*sin(obj.X(10))  cos(obj.X(11))*cos(obj.X(10))];
            J2=[1  sin(obj.X(10))*tan(obj.X(11))  cos(obj.X(10))*tan(obj.X(11));
                0  cos(obj.X(10))  -sin(obj.X(10));
                0  sin(obj.X(10))/cos(obj.X(11))  cos(obj.X(10))/cos(obj.X(11))];

            XHS=-(obj.W-obj.B)*sin(obj.X(11));
            YHS=(obj.W-obj.B)*cos(obj.X(11))*sin(obj.X(10));
            ZHS=(obj.W-obj.B)*cos(obj.X(11))*cos(obj.X(10));
            KHS=-obj.yg*obj.W*cos(obj.X(11))*cos(obj.X(10))-obj.zg*obj.W*cos(obj.X(11))*sin(obj.X(10));
            MHS=-obj.zg*obj.W*sin(obj.X(11))-obj.xg*obj.W*cos(obj.X(11))*cos(obj.X(10));
            NHS=-obj.xg*obj.W*cos(obj.X(11))*sin(obj.X(10))-obj.yg*obj.W*sin(obj.X(11));

            Xe=XHS+obj.Xuu*obj.X(1)*abs(obj.X(1))+...
                (obj.Xwq-obj.m)*obj.X(3)*obj.X(5)+(obj.Xqq+obj.m*obj.xg)*obj.X(5)^2+(obj.Xvr+obj.m)*obj.X(2)*obj.X(6)+(obj.Xrr+obj.m*obj.xg)*obj.X(6)^2-...
                obj.m*obj.yg*obj.X(4)*obj.X(5)-obj.m*obj.zg*obj.X(4)*obj.X(6)+Xprop;
            Ye=YHS+obj.Yvv*obj.X(2)*abs(obj.X(2))+obj.Yrr*obj.X(6)*abs(obj.X(6))+...
                obj.m*obj.yg*obj.X(6)^2+(obj.Yur-obj.m)*obj.X(1)*obj.X(6)+(obj.Ywp+obj.m)*obj.X(3)*obj.X(4)+(obj.Ypq-obj.m-obj.xg)*obj.X(4)*obj.X(5)+...
                obj.Yuv*obj.X(1)*obj.X(2)+obj.m*obj.yg*obj.X(4)^2+obj.m*obj.zg*obj.X(5)*obj.X(6)+obj.Yuudr*obj.X(1)^2*dr;
            Ze=ZHS+obj.Zww*obj.X(3)*abs(obj.X(3))+obj.Zqq*obj.X(5)*abs(obj.X(5))+...
                (obj.Zuq+obj.m)*obj.X(1)*obj.X(5)+(obj.Zvp-obj.m)*obj.X(2)*obj.X(4)+(obj.Zrp-obj.m*obj.xg)*obj.X(6)*obj.X(4)+obj.Zuw*obj.X(3)*obj.X(1)+...
                obj.m*obj.zg*(obj.X(4)^2+obj.X(5)^2)-obj.m*obj.yg*obj.X(6)*obj.X(5)+obj.Zuuds*obj.X(1)^2*ds;
            Ke=KHS+obj.Kpp*obj.X(4)*abs(obj.X(4))-...
                (obj.Izz-obj.Iyy)*obj.X(5)*obj.X(6)+obj.m*(obj.X(1)*obj.X(5)-obj.X(2)*obj.X(4))-obj.m*obj.zg*(obj.X(3)*obj.X(4)-obj.X(1)*obj.X(6))+Kprop;
            Me=MHS+obj.Mww*obj.X(3)*abs(obj.X(3))+obj.Mqq*obj.X(5)*abs(obj.X(5))+...
                (obj.Muq-obj.m*obj.xg)*obj.X(1)*obj.X(5)+(obj.Mvp+obj.m*obj.xg)*obj.X(2)*obj.X(4)+(obj.Mrp-(obj.Ixx-obj.Izz))*obj.X(6)*obj.X(4)+...
                obj.m*obj.zg*(obj.X(2)*obj.X(6)-obj.X(3)*obj.X(5))+obj.Muw*obj.X(1)*obj.X(3)+obj.Muuds*obj.X(1)^2*ds;
            Ne=NHS+obj.Nvv*obj.X(2)*abs(obj.X(2))+obj.Nrr*obj.X(6)*abs(obj.X(6))+...
                (obj.Nur-obj.m*obj.xg)*obj.X(1)*obj.X(6)+(obj.Nwp+obj.m*obj.xg)*obj.X(3)*obj.X(4)+(obj.Npq-(obj.Iyy-obj.Ixx))*obj.X(4)*obj.X(5)-...
                obj.m*obj.yg*(obj.X(2)*obj.X(6)-obj.X(3)*obj.X(5))+obj.Nuv*obj.X(1)*obj.X(2)+obj.Nuudr*obj.X(1)^2*dr;

            f=zeros(12, 1);
            a = obj.K;
            b = [Xe; Ye; Ze; Ke; Me; Ne];
            f(1:6)=obj.myGauss(a, b);
            f(7:9)=J1*obj.X(1:3);
            f(10:12)=J2*obj.X(4:6);
        end

        function euler(obj, dt, U)
            %METHOD1 euler's method
            obj.X=obj.X+kinetic_equation(obj, U).*dt;
        end

        function RK4(obj, dt, U)
            %METHOD1 Rongo Kuta method
            Xn=obj.X;
            K1=kinetic_equation(obj, U);
            obj.X=Xn+K1.*dt./2;
            K2=kinetic_equation(obj, U);
            obj.X=Xn+K2.*dt./2;
            K3=kinetic_equation(obj, U);
            obj.X=Xn+K3.*dt;
            K4=kinetic_equation(obj, U);
            obj.X=Xn+(K1+2*K2+2*K3+K4).*dt./6;
        end

        function [x] = myGauss(~, a, b) % 列主元高斯消去法
            len_n = length(a);
            x = zeros(len_n, 1);
            a = [a b];

            for k = 1:len_n - 1
                max = k;

                for i = k + 1:len_n

                    if a(i, k) > a(max, k)
                        max = i;
                    end

                end

                temp = a(k, k:len_n + 1);
                a(k, k:len_n + 1) = a(max, k:len_n + 1);
                a(max, k:len_n + 1) = temp;

                for i = k + 1:len_n
                    a(i, k) = -a(i, k) / a(k, k);
                    a(i, k + 1:len_n + 1) = a(i, k + 1:len_n + 1) + a(i, k) * a(k, k + 1:len_n + 1);
                end

            end

            x(len_n, 1) = a(len_n, len_n + 1) / a(len_n, len_n);

            for i = len_n - 1:-1:1
                sum = 0;

                for j = i + 1:len_n
                    sum = sum + x(j, 1) * a(i, j);
                end

                x(i, 1) = (a(i, len_n + 1) - sum) / a(i, i);
            end

        end
    end
end

