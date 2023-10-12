function pred_stats_z_scores

pathfolder = 'C:\Users\mszgm1\OneDrive - The University of Nottingham\LEL 2013';
addpath(genpath(pathfolder));

folder = sprintf('%s//analyses 2022',pathfolder);
filename = 'mTRF_encoding_Theta_RealPseudoBack_SpectrogramPhoneticsPhonemes_0-500ms_eegraw_64HzFs';
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
r_real_shuffled_std = []; r_pseudo_shuffled_std = []; r_back_shuffled_std = []; 
foldNum = 5;
for s = 1:length(a.predstats)
    for f = 1:foldNum
        r_real(s,f,:,1) = mean(Fisher(a.predstats{s}(f).real.r));
        r_pseudo(s,f,:,1) = mean(Fisher(a.predstats{s}(f).pseudo.r));
        r_back(s,f,:,1) = mean(Fisher(a.predstats{s}(f).back.r));
        r_real(s,f,:,2) = mean(Fisher(a.predstats{s}(f).real_shuffled.r));
        r_pseudo(s,f,:,2) = mean(Fisher(a.predstats{s}(f).pseudo_shuffled.r));
        r_back(s,f,:,2) = mean(Fisher(a.predstats{s}(f).back_shuffled.r));
        r_real_shuffled_std(s,f,:) = std(Fisher(a.predstats{s}(f).real_shuffled.r));
        r_pseudo_shuffled_std(s,f,:) = std(Fisher(a.predstats{s}(f).pseudo_shuffled.r));
        r_back_shuffled_std(s,f,:) = std(Fisher(a.predstats{s}(f).back_shuffled.r));
    end
end
% [(subjects)x(foldNum)x(electrodes)x(matched/shuffled)]

% r_real_z = (r_real(:,:,:,1)-r_real(:,:,:,2))./r_real_shuffled_std;
% r_pseudo_z = (r_pseudo(:,:,:,1)-r_pseudo(:,:,:,2))./r_pseudo_shuffled_std;
% r_back_z = (r_back(:,:,:,1)-r_back(:,:,:,2))./r_back_shuffled_std;
r_real_z = r_real(:,:,:,2);
r_pseudo_z = r_pseudo(:,:,:,2);
r_back_z = r_back(:,:,:,2);
r_real_stats = mean(mean(r_real_z(:,:,electrodes),3),2);
r_pseudo_stats = mean(mean(r_pseudo_z(:,:,electrodes),3),2);
r_back_stats = mean(mean(r_back_z(:,:,electrodes),3),2);

r_real_topoplot = squeeze(mean(mean(r_real_z),2));
r_pseudo_topoplot = squeeze(mean(mean(r_pseudo_z),2));
r_back_topoplot = squeeze(mean(mean(r_back_z),2));
size(r_real_topoplot)
cmin = -0.6; cmax = 0.6;
figure(1); clf;
subplot(1,3,1);
topoplot(r_real_topoplot, 'channel_loc.elp', 'maplimits', [cmin cmax]); 
colormap(jet(128)); caxis([cmin cmax]); colorbar;
subplot(1,3,2);
topoplot(r_pseudo_topoplot, 'channel_loc.elp', 'maplimits', [cmin cmax]); 
colormap(jet(128)); caxis([cmin cmax]); colorbar;
subplot(1,3,3);
topoplot(r_back_topoplot, 'channel_loc.elp', 'maplimits', [cmin cmax]); 
colormap(jet(128)); caxis([cmin cmax]); colorbar;

p_real_vs_pseudo_bts = bootstrapping(r_real_stats-r_pseudo_stats,10000);
p_real_vs_back_bts = bootstrapping(r_real_stats-r_back_stats,10000);
p_pseudo_vs_back_bts = bootstrapping(r_pseudo_stats-r_back_stats,10000);


[~,~,~,adj_p_bts]=fdr_bh([p_real_vs_pseudo_bts,p_real_vs_back_bts,p_pseudo_vs_back_bts]);
[mean(r_real_stats),mean(r_pseudo_stats),mean(r_back_stats);adj_p_bts]

mean_r_real = mean(r_real_stats); err_real = std(r_real_stats)/(length(r_real_stats))^0.5;
mean_r_pseudo = mean(r_pseudo_stats); err_pseudo = std(r_pseudo_stats)/(length(r_pseudo_stats))^0.5;
mean_r_back = mean(r_back_stats); err_back = std(r_back_stats)/(length(r_back_stats))^0.5;
mean_r = [mean_r_real,mean_r_pseudo,mean_r_back];
err = [err_real,err_pseudo,err_back];

figure(2); clf;
b = bar(mean_r,'LineWidth',4); 
b.FaceColor = 'flat';
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 0 1];
b.CData(3,:) = [0 1 0];
hold on;
errorbar([1,2,3],mean_r,err,'k','linestyle','none','LineWidth',4);

% for s = 1:length(a.predstats)
%     for e = 1:length(electrodes)
%         r_real(e,s) = mean(Fisher(a.predstats{s}(e).real.r(:,e)))';
%         r_pseudo(e,s) = mean(Fisher(a.predstats{s}(e).pseudo.r(:,e)))';
%         r_back(e,s) = mean(Fisher(a.predstats{s}(e).back.r(:,e)))';
%     end
% end
% 
% mean(mean(r_real))
% mean(mean(r_pseudo))
% mean(mean(r_back))
% [~,p_real_vs_pseudo] = ttest(mean(r_real),mean(r_pseudo));
% [~,p_real_vs_back] = ttest(mean(r_real),mean(r_back));
% [~,p_pseudo_vs_back] = ttest(mean(r_pseudo),mean(r_back));
% p_real_vs_pseudo
% p_real_vs_back
% p_pseudo_vs_back


function Fisher_r = Fisher(r)
Fisher_r = 0.5*log((1+r)./(1-r));
end

end