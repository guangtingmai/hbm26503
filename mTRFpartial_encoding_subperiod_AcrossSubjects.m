function mTRFpartial_encoding_subperiod_AcrossSubjects(comparison,freqName,targetName,partialName)

folder = 'C:\Users\mszgm1\OneDrive - The University of Nottingham\LEL 2013';
addpath(genpath(folder));
subfolder = sprintf('%s//analyses 2022',folder);
name1 = sprintf('mTRF_encoding_%s_RealPseudoBack_%s(target)_%s(partialed)2',freqName,targetName,partialName);
name2 = '0-500ms_eegraw_64HzFs';

loadname = sprintf('%s//%s_%s.mat',subfolder,name1,name2);
a = load(loadname);
eeg_partial = a.eeg_partial;
targetStimFea = a.targetStimFea;
bestlambda = a.bestlambda;

subjectID1 = '03_20am';  subjectID2 = '03_20eve'; subjectID3 = '03_22pm';   subjectID4 = '03_22eve'; subjectID5 = '03_24am';  
subjectID6 = '03_24pm';  subjectID7 = '03_24eve'; subjectID8 = '03_26pm';   subjectID9 = '03_29am';  subjectID10 = '03_29pm'; 
subjectID11 = '03_30am'; subjectID12 = '03_30pm'; subjectID13 = '03_31am';  subjectID14 = '03_31pm'; subjectID15 = '04_02eve'; 
subjectID16 = '04_06am'; subjectID17 = '04_06pm'; subjectID18 = '04_06eve'; subjectID19 = '04_07am'; subjectID20 = '04_07pm';
subjectID = {subjectID1, subjectID2, subjectID3, subjectID4, subjectID5, subjectID6, subjectID7, subjectID8, subjectID9, subjectID10, ...
     subjectID11, subjectID12, subjectID13, subjectID14, subjectID15, subjectID16, subjectID17, subjectID18, subjectID19, subjectID20};


downsamplefactor = 4;
fs = 256/downsamplefactor;

predstats = {};
dirTRF = 1;
subperiods = {};
% subperiods{1} = [0,50]; subperiods{2} = [50,100]; subperiods{3} = [100,150]; subperiods{4} = [150,200]; subperiods{5} = [200,250]; 
% subperiods{6} = [250,300]; subperiods{7} = [300,350]; subperiods{8} = [350,400]; subperiods{9} = [400,450]; subperiods{10} = [450,500];
for n = 1:0.5/(1/fs)
    subperiods{n} = [(n/fs)*1000-15,(n/fs)*1000+15];
end

%savefolder = 'C:\Users\mszgm1\OneDrive - University College London\LEL data analyses';
shuffle_or_not = 0;
for s = 16:length(subjectID)
    [pred,predstats(s,:)] = mTRFpartial_subperiod_SingleSubject...
      (bestlambda{s}, eeg_partial{s}, targetStimFea{s}, comparison, fs, dirTRF, subperiods, shuffle_or_not);
    savename = sprintf('%s//%s_subperiods_pred_%s.mat',subfolder,name1,subjectID{s});
    save(savename,'pred');
    s
end

% savename = sprintf('%s//%s_subperiods.mat',subfolder,name1);
% save(savename,'predstats');