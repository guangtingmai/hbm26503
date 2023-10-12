function mTRF_encoding_subperiod_AcrossSubjects(freqRange, comparison)

folder = 'C:\Users\mszgm1\OneDrive - The University of Nottingham\LEL 2013';
addpath(genpath(folder));
subfolder = sprintf('%s//analyses 2022',folder);

switch comparison
    case 1 % real vs pseudo vs back
        Comparison = 'RealPseudoBack';
    otherwise % real vs pseudo
        Comparison = 'RealPseudo';
end
switch freqRange
    case 1 % Delta
        FreqRange = 'Delta';
    case 2 % Theta
        FreqRange = 'Theta';
    otherwise % DeltaTheta
        FreqRange = 'DeltaTheta';
end
stimIdx = [4,7,8];
name1 = sprintf('mTRF_encoding_%s_%s_SpectrogramPhoneticsPhonemes',FreqRange,Comparison);
name2 = '0-500ms_eegraw_64HzFs';

loadname = sprintf('%s//%s_%s.mat',subfolder,name1,name2);
a = load(loadname);
bestlambda = a.bestlambda;

subjectID1 = '03_20am';  subjectID2 = '03_20eve'; subjectID3 = '03_22pm';   subjectID4 = '03_22eve'; subjectID5 = '03_24am';  
subjectID6 = '03_24pm';  subjectID7 = '03_24eve'; subjectID8 = '03_26pm';   subjectID9 = '03_29am';  subjectID10 = '03_29pm'; 
subjectID11 = '03_30am'; subjectID12 = '03_30pm'; subjectID13 = '03_31am';  subjectID14 = '03_31pm'; subjectID15 = '04_02eve'; 
subjectID16 = '04_06am'; subjectID17 = '04_06pm'; subjectID18 = '04_06eve'; subjectID19 = '04_07am'; subjectID20 = '04_07pm';
subjectID = {subjectID1, subjectID2, subjectID3, subjectID4, subjectID5, subjectID6, subjectID7, subjectID8, subjectID9, subjectID10, ...
     subjectID11, subjectID12, subjectID13, subjectID14, subjectID15, subjectID16, subjectID17, subjectID18, subjectID19, subjectID20};

downsamplefactor = 4;
fs = 256/downsamplefactor;

% pred = {};
% predstats = {};
dirTRF = 1;
subperiods = {};
% subperiod{1} = [0,50]; subperiod{2} = [50,100]; subperiod{3} = [100,150]; subperiod{4} = [150,200]; subperiod{5} = [200,250]; 
% subperiod{6} = [250,300]; subperiod{7} = [300,350]; subperiod{8} = [350,400]; subperiod{9} = [400,450]; subperiod{10} = [450,500];
for n = 1:0.5/(1/fs)
    subperiods{n} = [(n/fs)*1000-15,(n/fs)*1000+15];
end
for s = 18:length(subjectID)
    [pred,predstats] = mTRF_subperiod_SingleSubject...
      (subfolder, subjectID{s}, bestlambda{s}, comparison, freqRange, stimIdx, fs, downsamplefactor, dirTRF, subperiods, 'raw');
    savename = sprintf('%s//%s_subperiods_pred_%s.mat',subfolder,name1,subjectID{s});
    save(savename,'pred');
    s
end


end
