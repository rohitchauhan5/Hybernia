function [ ind_plat ] = find_plateau( y_var, x )
%FIND_PLATEAU ~ at what index does the graph start to plateu?
i = 1;

%while ((abs(y_var(length(y_var)) - y_var(i))/(x(length(y_var))-x(i)) > abs((y_var(length(y_var))-y_var(1))/4430)) & (abs(y_var(length(y_var)) - y_var(i))/(x(length(y_var))-x(i))) > 0.0005) %break when average slope, scaled to range, is less than some threshold
while (abs((y_var(length(y_var)) - y_var(i))/(x(length(y_var))-x(i))) > abs((y_var(length(y_var))-y_var(1))/(x(length(y_var))-x(1)))/28.92 & abs((y_var(length(y_var)) - y_var(i))/(x(length(y_var))-x(i))) > 0.0005) %break when average slope, scaled to range, is less than some threshold

%     disp('time')
%     x(i)
%     disp('slope')
%     abs(y_var(length(y_var)) - y_var(i))/(x(length(y_var))-x(i))
%     disp('comp')
%     abs((y_var(length(y_var))-y_var(1))/4430)
%     pause
    i=i+1;
    if i == length(y_var)
        break
    end
end

% using differentials (below) didn't work because plateaus may be too
% bumpy!
% dervs = diff(y_var);
% while all(abs(dervs(i:length(dervs),1))<0.1)==0
%     i = i+1;
%     if i == length(dervs)+1
%         break
%     end
% end
ind_plat = i;
end

