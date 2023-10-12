function xq = plot_graph_se_signficance(y,fs,colour,SigLineVal,SigLineRange,SigLineColour)

% y [(subjects)x(time_pts)]

x = ((1:size(y,2))/fs)*1000;
xq = ((1:0.5:size(y,2))/fs)*1000;

y_se = std(y)/((size(y,1))^0.5);
y = interp1(x,mean(y),xq);
y_se = interp1(x,y_se,xq);
y_l = y - y_se;
y_u = y + y_se;

fill([xq,fliplr(xq)],[y_l,fliplr(y_u)],'k','LineWidth',0.01,'EdgeColor','none'); hold on;
alpha(0.25);
plot(xq,y,'Color',colour,'LineWidth',3); hold on;
alpha(0.25);

for s = 1:length(SigLineRange)
    for xx = 1:length(SigLineRange{s})
        if SigLineRange{s}(xx)==1
            x = xq(2*SigLineRange{s}(xx)-1:2*SigLineRange{s}(xx));        
        elseif SigLineRange{s}(xx)==size(y,2)
            x = xq(2*SigLineRange{s}(xx)-2:2*SigLineRange{s}(xx)-1);
        else
            x = xq(2*SigLineRange{s}(xx)-2:2*SigLineRange{s}(xx));
        end
        plot(x,SigLineVal(s)*ones(1,length(x)),'Color',SigLineColour{s},'LineWidth',3); hold on;
    end
end

