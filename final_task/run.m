function result = run(total)
%TEST_F 此处显示有关此函数的摘要
%   此处显示详细说明
a=Vehicle(Vehicle_profile());
controller=Controller();
recorder=Recorder(a.X);
for i=1:total
    output = controller.output(a.X);
    recorder.save_data(a.X, output(1));
    a.RK4(output(1), output(2:5))
end
result = recorder.result;
end

