result = run(10000);  % 计算仿真结果，可将run编译为MEX以减少计算时间
%% 设置图窗
f=figure;
fig_Axes=gca(f);
fig_Axes.View=[-20, 30];
fig_Axes.ZDir="reverse";
hold on;
%% 
axis(fig_Axes,'equal');
ylim(fig_Axes, [-20, 20]);
zlim(fig_Axes, [-20, 20]);
fig_line=plot3(fig_Axes,result(7,:), result(8,:), result(9,:), 'r:','LineWidth', 1.5);  % 画路径