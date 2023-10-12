function pred_stats_subperiods4plots(yRange,sigLine,filename)

pathfolder = 'C:\Users\mszgm1\OneDrive - The University of Nottingham\LEL 2013';
addpath(genpath(pathfolder));

folder = sprintf('%s//analyses 2022',pathfolder);
filename = sprintf('%s_subperiods',filename);
filename = sprintf('%s//%s.mat',folder,filename);
a = load(filename);

% {'Fp1', 'AF3', 'F7', 'F3', 'FC1', 'FC5', 'T7', 'C3', ...
%    1      2     3     4      5      6     7     8  
%  'CP1', 'CP5', 'P7', 'P3', 'Pz', 'PO3', 'O1', 'Oz', ...
%    9      10    11    12    13    14     15    16
%  'O2', 'PO4', 'P4', 'P8', 'CP6', 'CP2', 'C4', 'T8', ...
%   17     18    19    20     21    22     23    24
%  'FC6', 'FC2', 'F4', 'F8', 'AF4', 'Fp2', 'Fz', 'Cz'};
%   25     26     27    28    29     30     31    32

%electrodes = [1:32];
electrodes = [1:10,21:32];
r_real = []; r_pseudo = []; r_back = [];
r_real_shuffled = []; r_pseudo_shuffled = []; r_back_shuffled = [];
foldNum = 5;
subperiodNum = size(a.predstats,2);
for s = 1:20
    for p = 1:subperiodNum
        for f = 1:foldNum
            r_real(s,p,f) = mean(mean(Fisher(a.predstats{s,p}(f).real.r(:,electrodes))));
            r_pseudo(s,p,f) = mean(mean(Fisher(a.predstats{s,p}(f).pseudo.r(:,electrodes))));
            r_back(s,p,f) = mean(mean(Fisher(a.predstats{s,p}(f).back.r(:,electrodes))));
%             r_real_shuffled(s,p,f) = mean(mean(Fisher(a.predstats{s,p}(f).real_shuffled.r(:,electrodes))));
%             r_pseudo_shuffled(s,p,f) = mean(mean(Fisher(a.predstats{s,p}(f).pseudo_shuffled.r(:,electrodes))));
%             r_back_shuffled(s,p,f) = mean(mean(Fisher(a.predstats{s,p}(f).back_shuffled.r(:,electrodes))));
        end
    end
end

r_real_stats = mean(r_real,3); 
r_pseudo_stats = mean(r_pseudo,3); 
r_back_stats = mean(r_back,3);
% r_real_shuffled_stats = mean(r_real_shuffled,3); 
% r_pseudo_shuffled_stats = mean(r_pseudo_shuffled,3);
% r_back_shuffled_stats = mean(r_back_shuffled,3);

alpha1 = 0.05/(subperiodNum*3);
real_vs_pseudo = []; real_vs_back = []; pseudo_vs_back = [];
real_vs_shuffled = []; pseudo_vs_shuffled = []; back_vs_shuffled = [];
for p = 1:subperiodNum
%      [~,real_vs_pseudo(p)] = ttest(r_real_stats(:,p),r_pseudo_stats(:,p));
%      [~,real_vs_back(p)] = ttest(r_real_stats(:,p),r_back_stats(:,p));
%      [~,pseudo_vs_back(p)] = ttest(r_pseudo_stats(:,p),r_back_stats(:,p));
%     [~,real_vs_shuffled(p)] = ttest(r_real_stats(:,p),r_real_shuffled_stats(:,p));
%     [~,pseudo_vs_shuffled(p)] = ttest(r_pseudo_stats(:,p),r_pseudo_shuffled_stats(:,p));
%     [~,back_vs_shuffled(p)] = ttest(r_back_stats(:,p),r_back_shuffled_stats(:,p));
    real_vs_pseudo(p) = bootstrapping(r_real_stats(:,p)-r_pseudo_stats(:,p),10000);
    real_vs_back(p) = bootstrapping(r_real_stats(:,p)-r_back_stats(:,p),10000);
    pseudo_vs_back(p) = bootstrapping(r_pseudo_stats(:,p)-r_back_stats(:,p),10000);
%     real_vs_pseudo(p) = signrank(r_real_stats(:,p)-r_pseudo_stats(:,p));
%     real_vs_back(p) = signrank(r_real_stats(:,p)-r_back_stats(:,p));
%     pseudo_vs_back(p) = signrank(r_pseudo_stats(:,p)-r_back_stats(:,p));
%     real_vs_shuffled(p) = bootstrapping(r_real_stats(:,p)-r_real_shuffled_stats(:,p),10000);
%     pseudo_vs_shuffled(p) = bootstrapping(r_pseudo_stats(:,p)-r_pseudo_shuffled_stats(:,p),10000);
%     back_vs_shuffled(p) = bootstrapping(r_back_stats(:,p)-r_back_shuffled_stats(:,p),10000);

%     ci_real_vs_pseudo = bootci(10000,{@(x)(mean(x)),r_real_stats(:,p)-r_pseudo_stats(:,p)},'Alpha',alpha1,'type','bca');
%     if max(ci_real_vs_pseudo)*min(ci_real_vs_pseudo)>0
%         real_vs_pseudo(p) = 1;
%     else
%         real_vs_pseudo(p) = 0;
%     end
%     ci_real_vs_back = bootci(10000,{@(x)(mean(x)),r_real_stats(:,p)-r_back_stats(:,p)},'Alpha',alpha1,'type','bca');
%     if max(ci_real_vs_back)*min(ci_real_vs_back)>0
%         real_vs_back(p) = 1;
%     else
%         real_vs_back(p) = 0;
%     end
%     ci_pseudo_vs_back = bootci(10000,{@(x)(mean(x)),r_pseudo_stats(:,p)-r_back_stats(:,p)},'Alpha',alpha1,'type','bca');
%     if max(ci_pseudo_vs_back)*min(ci_pseudo_vs_back)>0
%         pseudo_vs_back(p) = 1;
%     else
%         pseudo_vs_back(p) = 0;
%     end
end
p = [real_vs_pseudo,real_vs_back,pseudo_vs_back];
[~,~,~,adj_p]=fdr_bh(p);
[adj_p(1:length(adj_p)/3);adj_p(length(adj_p)/3+1:2*length(adj_p)/3);adj_p(2*length(adj_p)/3+1:end)]
% p = [real_vs_shuffled,pseudo_vs_shuffled,back_vs_shuffled];
%[real_vs_pseudo;real_vs_back;pseudo_vs_back]


x = ((1:subperiodNum)/64)*1000;
xq = ((1:0.5:subperiodNum)/64)*1000;
x_real_vs_pseudo = find(adj_p(1:length(adj_p)/3)<0.05);
x_real_vs_back = find(adj_p(length(adj_p)/3+1:2*length(adj_p)/3)<0.05); 
x_pseudo_vs_back = find(adj_p(2*length(adj_p)/3+1:end)<0.05);
% se
r_real_stats_se = (std(r_real_stats))/(s^0.5); r_real_stats_se = interp1(x,r_real_stats_se,xq);
r_pseudo_stats_se = (std(r_pseudo_stats))/(s^0.5); r_pseudo_stats_se = interp1(x,r_pseudo_stats_se,xq);
r_back_stats_se = (std(r_back_stats))/(s^0.5); r_back_stats_se = interp1(x,r_back_stats_se,xq);
% mean
r_real_stats = interp1(x,mean(r_real_stats),xq);
r_pseudo_stats = interp1(x,mean(r_pseudo_stats),xq);
r_back_stats = interp1(x,mean(r_back_stats),xq);
% upper and lower bounds
r_real_l = r_real_stats - r_real_stats_se; r_real_u = r_real_stats + r_real_stats_se;
r_pseudo_l = r_pseudo_stats - r_pseudo_stats_se; r_pseudo_u = r_pseudo_stats + r_pseudo_stats_se;
r_back_l = r_back_stats - r_back_stats_se; r_back_u = r_back_stats + r_back_stats_se;

%figure(1); clf;
fill([xq,fliplr(xq)],[r_real_l,fliplr(r_real_u)],'r','LineWidth',0.01,'EdgeColor','none'); hold on;
fill([xq,fliplr(xq)],[r_pseudo_l,fliplr(r_pseudo_u)],'b','LineWidth',0.01,'EdgeColor','none'); hold on;
fill([xq,fliplr(xq)],[r_back_l,fliplr(r_back_u)],'g','LineWidth',0.01,'EdgeColor','none'); hold on;
alpha(0.25);
plot(xq,r_real_stats,'Color',[1 0 0 0.5],'LineWidth',3); hold on;
plot(xq,r_pseudo_stats,'Color',[0 0 1 0.5],'LineWidth',3); hold on;
plot(xq,r_back_stats,'Color',[0 1 0 0.5],'LineWidth',3); hold on;
alpha(0.2);

for xx = 1:length(x_real_vs_pseudo)
    if 2*x_real_vs_pseudo(xx)==2
        x = xq(2*x_real_vs_pseudo(xx)-1:2*x_real_vs_pseudo(xx));        
    elseif 2*x_real_vs_pseudo(xx)>62
        x = xq(2*x_real_vs_pseudo(xx)-2:2*x_real_vs_pseudo(xx)-1);
    else
        x = xq(2*x_real_vs_pseudo(xx)-2:2*x_real_vs_pseudo(xx));
    end
    plot(x,sigLine(1)*ones(1,length(x)),'Color','#6E2C00','LineWidth',3); hold on;
end
for xx = 1:length(x_real_vs_back)
    if 2*x_real_vs_back(xx)==2
        x = xq(2*x_real_vs_back(xx)-1:2*x_real_vs_back(xx));
    elseif 2*x_real_vs_back(xx)>62
        x = xq(2*x_real_vs_back(xx)-2:2*x_real_vs_back(xx)-1);
    else
        x = xq(2*x_real_vs_back(xx)-2:2*x_real_vs_back(xx));
    end
    plot(x,sigLine(2)*ones(1,length(x)),'Color','#DB7400','LineWidth',3); hold on;
end
for xx = 1:length(x_pseudo_vs_back)
    if 2*x_pseudo_vs_back(xx)==2
        x = xq(2*x_pseudo_vs_back(xx)-1:2*x_pseudo_vs_back(xx));
    elseif 2*x_pseudo_vs_back(xx)>62
        x = xq(2*x_pseudo_vs_back(xx)-2:2*x_pseudo_vs_back(xx)-1);
    else
        x = xq(2*x_pseudo_vs_back(xx)-2:2*x_pseudo_vs_back(xx));
    end
    plot(x,sigLine(3)*ones(1,length(x)),'Color','#EDBB99','LineWidth',3);
end
ylim(yRange);
% hold on;
% plot(x,mean(r_real_shuffled_stats),'--r'); hold on;
% plot(x,mean(r_pseudo_shuffled_stats),'--b'); hold on;
% plot(x,mean(r_back_shuffled_stats),'--g');
xlabel('Time lag (ms)');
%ylabel('Reconstruction Accuracy (Fisher-z)');

ax = gca;
ax.Box = 'off';
ax.FontSize = 12;
ax.LineWidth = 1.5;

% mean_r_real = mean(r_real_stats); err_real = std(r_real_stats)/(length(r_real_stats))^0.5;
% mean_r_pseudo = mean(r_pseudo_stats); err_pseudo = std(r_pseudo_stats)/(length(r_pseudo_stats))^0.5;
% mean_r_back = mean(r_back_stats); err_back = std(r_back_stats)/(length(r_back_stats))^0.5;
% mean_r = [mean_r_real,mean_r_pseudo,mean_r_back];
% err = [err_real,err_pseudo,err_back];
% figure(2); clf;
% b = bar(mean_r,'LineWidth',4); 
% b.FaceColor = 'flat';
% b.CData(1,:) = [1 0 0];
% b.CData(2,:) = [0 0 1];
% b.CData(3,:) = [0 1 0];
% hold on;
% errorbar([1,2,3],mean_r,err,'k','linestyle','none','LineWidth',4);




function Fisher_r = Fisher(r)
Fisher_r = 0.5*log((1+r)./(1-r));
end

end