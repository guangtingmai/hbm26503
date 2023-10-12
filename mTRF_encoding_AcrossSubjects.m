function mTRF_encoding_AcrossSubjects(freqRange, comparison, stimIdx, ica_or_raw, trialpercent)

pathfolder = 'C:\Users\mszgm1\OneDrive - The University of Nottingham\LEL 2013';
addpath(genpath(pathfolder));

subjectID1 = '03_20am';  subjectID2 = '03_20eve'; subjectID3 = '03_22pm';   subjectID4 = '03_22eve'; subjectID5 = '03_24am';  
subjectID6 = '03_24pm';  subjectID7 = '03_24eve'; subjectID8 = '03_26pm';   subjectID9 = '03_29am';  subjectID10 = '03_29pm'; 
subjectID11 = '03_30am'; subjectID12 = '03_30pm'; subjectID13 = '03_31am';  subjectID14 = '03_31pm'; subjectID15 = '04_02eve'; 
subjectID16 = '04_06am'; subjectID17 = '04_06pm'; subjectID18 = '04_06eve'; subjectID19 = '04_07am'; subjectID20 = '04_07pm';
subjectID = {subjectID1, subjectID2, subjectID3, subjectID4, subjectID5, subjectID6, subjectID7, subjectID8, subjectID9, subjectID10, ...
     subjectID11, subjectID12, subjectID13, subjectID14, subjectID15, subjectID16, subjectID17, subjectID18, subjectID19, subjectID20};

folder = sprintf('%s//analyses 2022',pathfolder);

downsamplefactor = 4;
new_fs = 256/downsamplefactor;

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

model = {};
bestlambda = {};
pred = {};
predstats = {};
dirTRF = 1; % encoding
for s = 1:length(subjectID)
    [model{s},bestlambda{s},pred{s},predstats{s},stimFeaName,tmin,tmax] = mTRF_SingleSubject...
        (folder, subjectID{s}, freqRange, comparison, dirTRF, stimIdx, ica_or_raw, downsamplefactor,trialpercent);
    s
end

savename = sprintf('%s//mTRF_encoding_%s_%s_%s_%d-%dms_eeg%s_%dHzFs_TrialPercent%d.mat',folder,FreqRange,...
    Comparison,stimFeaName,tmin,tmax,ica_or_raw,new_fs,trialpercent);
save(savename,'model','bestlambda','pred','predstats');