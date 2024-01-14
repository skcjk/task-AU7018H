classdef Controller < handle
    %CONTROLLER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        E=0;
    end
    
    methods
        
        function output = output(obj, X)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            e=-X(9);
            e2=X(3);
            obj.E=obj.E+e;
            c=0.05*e+3*e2+0.0001*obj.E;
            dt=0.01;
            output=[dt, -0.0175, c, 3.86, -5.43e-1];
            
        end
    end
end

