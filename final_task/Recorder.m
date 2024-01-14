classdef Recorder < handle
    %RECORDER 用于记录仿真结果
    
    properties
        result;
        initial_size = 2000;  % 初始数组大小
        increment = 1000;     % 增量
        current_index = 1;
        len_of_X;
        t=0;
    end
    
    methods
        function obj = Recorder(X)
            %RECORDER 构造此类的实例
            %   此处显示详细说明
            obj.len_of_X = size(X, 1)+1;
            obj.result = zeros(obj.len_of_X, obj.initial_size);
        end
        
        function check_resualt_size(obj)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            current_size = size(obj.result, 2);
            if obj.current_index > current_size
                new_resualt = zeros(obj.len_of_X, current_size+obj.increment);
                new_resualt(:,1:current_size)=obj.result;
                obj.result = new_resualt;
            end
        end
        
        function save_data(obj, X, dt)
            obj.check_resualt_size();
            obj.result(1:obj.len_of_X-1, obj.current_index)=X;
            obj.result(obj.len_of_X, obj.current_index)=obj.t;
            obj.t=obj.t+dt;
            obj.current_index = obj.current_index +1;
        end

    end
end

