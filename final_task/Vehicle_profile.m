classdef Vehicle_profile
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明

    properties
        %% 本体参数
        xg=0.0; yg=0.0; zg=1.96e-2;  % Center of Gravity wrt Origin at CB
        Ixx=1.77e-1; Iyy=3.45; Izz=3.45;  % Moments of Inertia wrt Origin at CB
        rho=1.03e3; W=2.99e2; B=3.06e2;
        %% STD REMUS Non-Linear Maneuvering Coefficients: Forces
        Xuu=-1.62;  % Cross-flow Drag
        Xudot=-9.3e-1;  % Added Mass
        Xwq=-3.55e1;  % Added Mass Cross-term
        Xqq=-1.93;  % Added Mass Cross-term
        Xvr=3.55e1;  % Added Mass Cross-term
        Xrr=-1.93;  % Added Mass Cross-term
        Yvv=-1.31e2;  % Cross-flow Drag
        Yrr=6.32e-1;  % Cross-flow Drag
        Yuv=-2.86e1;  % Body Lift Force and Fin Lift
        Yvdot=-3.55e1;  % Added Mass
        Yrdot=1.93;  % Added Mass
        Yur=5.22;  % Added Mass Cross Term and Fin Lift
        Ywp=3.55e1;  % Added Mass Cross-term
        Ypq=1.93;  % Added Mass Cross-term
        Yuudr=9.64;  % Fin Lift Force
        Zww=-1.31e2;  % Cross-flow Drag
        Zqq=-6.32e-1;  % Cross-flow Drag
        Zuw=-2.86e1;  % Body Lift Force and Fin Lift
        Zwdot=-3.55e1;  % Added Mass
        Zqdot=-1.93;  % Added Mass
        Zuq=-5.22;  % Added Mass Cross Term and Fin Lift
        Zvp=-3.55e1;  % Added Mass Cross Term
        Zrp=1.93;  % Added Mass Cross Term
        Zuuds=-9.64;  % Fin Lift Force
        %% STD REMUS Non-Linear Maneuvering Coefficients: Moments
        Kpp=-1.3e-3;  % Rolling Resistance
        Kpdot=-1.41e-2;  % Added Mass
        Mww=3.18;  % Cross-flow Drag
        Mqq=-9.4;  % Cross-flow Drag
        Muw=2.4e1;  % Body and Fin Lift and Munk Moment
        Mwdot=-1.93;  % Added Mass
        Mqdot=-4.88;  % Added Mass
        Muq=-2;  % Added Mass Cross Term and Fin Lift
        Mvp=-1.93;  % Added Mass Cross Term
        Mrp=4.86;  % Added Mass Cross Term
        Muuds=-6.15;  % Fin Lift Moment
        Nvv=-3.18;  % Cross-flow DragKprop
        Nrr=-9.4;  % Cross-flow Drag
        Nuv=-2.4e1;  % Body and Fin Lift and Munk Moment
        Nvdot=1.93;  % Added Mass
        Nrdot=-4.88;  % Added Mass
        Nur=-2;  % Added Mass Cross Term and Fin Lift
        Nwp=-1.93;  % Added Mass Cross Term
        Npq=-4.86;  % Added Mass Cross Term
        Nuudr=-6.15;  % Fin Lift Moment
    end

end

