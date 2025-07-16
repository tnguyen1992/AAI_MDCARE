% %%% DCARE Pipeline
%%%%%%%%%%%%%%%%%% ARC Project %%%%%%%%%%%%%%%%%%%%%%%
% convert raw data to homer2 nirs format

% wo sind unsere Daten?
srcPath = 'P:\projects\DCARE\DCARE\rawData\';

% wohin sollen die Daten gespeichert werden?
% desPath = 'P:\projects\DCARE\DCARE\procData\nirsData\';

%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d']);
  end

  %% Konvertierung
  
  for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/']));
  
%   Sub1DesFile = strcat(desPath, sprintf(['DCARE_%02d_sub1'], ...
%                       i), '.nirs');
%   Sub2DesFile = strcat(desPath, sprintf(['DCARE_%02d_sub2'], ...
%                       i), '.nirs');
% 
% % -------------------------------------------------------------------------
% % Load SD file
% % -------------------------------------------------------------------------
% SDfile = 'P:\projects\DCARE\DCARE\MATLAB\CARE.SD';
% load(SDfile, '-mat', 'SD');

% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

%   Sub1_wl1File = strcat(Sub1SrcDir, sprintf(['DCARE_%02d'], i), '.wl1');
%   Sub1_wl2File = strcat(Sub1SrcDir, sprintf(['DCARE_%02d'], i), '.wl2');
%   Sub1_hdrFile = strcat(Sub1SrcDir, sprintf(['DCARE_%02d'], i), '.hdr');
%  
%   F1=[Sub1_wl1File; Sub1_wl2File; Sub1_hdrFile];
%   spm_fnirs_read_nirscout(F1);

  
  Sub2_wl1File = strcat(Sub2SrcDir, sprintf(['DCARE_%02d'], i), '.wl1');
  Sub2_wl2File = strcat(Sub2SrcDir, sprintf(['DCARE_%02d'], i), '.wl2');
  Sub2_hdrFile = strcat(Sub2SrcDir, sprintf(['DCARE_%02d'], i), '.hdr');
  
  F2=[Sub2_wl1File; Sub2_wl2File; Sub2_hdrFile];
  spm_fnirs_read_nirscout(F2);
  
  

  end
  
  

%% convert nirs file to concentrations
clear all
% wo sind unsere Daten?
srcPath = 'P:\projects\DCARE\DCARE\rawData\';


%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d']);
  end

  %% Konvertierung
  
  for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_NIRSFile = strcat(Sub1SrcDir, sprintf(['NIRS.mat']));
%   Sub1_chFile = strcat(Sub1SrcDir, sprintf(['ch_config.txt']));

  Sub2_NIRSFile = strcat(Sub2SrcDir, sprintf(['NIRS.mat']));
%   Sub2_chFile = strcat(Sub2SrcDir, sprintf(['ch_config.txt']));

%   FC1=[Sub1_NIRSFile; Sub1_chFile];
  
%     load(Sub1_NIRSFile);
%     age=5;
%     spm_fnirs_convert_ui(Sub1_NIRSFile, P, age)
    
%   FC2=[Sub2_NIRSFile; Sub2_chFile];
%     clear P
    load(Sub2_NIRSFile);
    age=36;
    spm_fnirs_convert_ui(Sub2_NIRSFile, P, age)

  end
  
  %% convert nirs file to concentrations
clear all
% wo sind unsere Daten?
srcPath = 'P:\projects\DCARE\DCARE\rawData\';


%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d']);
  end

  %% Preprocessing
  
  for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_NIRSFile = strcat(Sub1SrcDir, sprintf(['NIRS.mat']));
%   Sub1_chFile = strcat(Sub1SrcDir, sprintf(['ch_config.txt']));

  Sub2_NIRSFile = strcat(Sub2SrcDir, sprintf(['NIRS.mat']));
%   Sub2_chFile = strcat(Sub2SrcDir, sprintf(['ch_config.txt']));

%   FC1=[Sub1_NIRSFile; Sub1_chFile];
  
%     spm_fnirs_temporalpreproc_ui(Sub1_NIRSFile)
    
%   FC2=[Sub2_NIRSFile; Sub2_chFile];
    spm_fnirs_temporalpreproc_ui(Sub2_NIRSFile)

  end
  
  
  clear all

%% 1st level estimation  
% wo sind unsere Daten?
srcPath = 'P:\projects\DCARE\DCARE\rawData\';


%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d']);
  end

  %% 1st level params
  
  for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_NIRSFile = strcat(Sub1SrcDir, sprintf(['NIRS.mat']));
%   Sub1_chFile = strcat(Sub1SrcDir, sprintf(['ch_config.txt']));

  Sub2_NIRSFile = strcat(Sub2SrcDir, sprintf(['NIRS.mat']));
%   Sub2_chFile = strcat(Sub2SrcDir, sprintf(['ch_config.txt']));

%   FC1=[Sub1_NIRSFile; Sub1_chFile];
  
    spm_fnirs_specify1st_ui(Sub1_NIRSFile)
    
%   FC2=[Sub2_NIRSFile; Sub2_chFile];
    spm_fnirs_specify1st_ui(Sub2_NIRSFile)

  end

  %% SPM estimation
  
  for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/HbO/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/HbO/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_SPMFile = strcat(Sub1SrcDir, sprintf(['SPM.mat']));
  Sub2_SPMFile = strcat(Sub2SrcDir, sprintf(['SPM.mat']));
  
    spm_fnirs_estimate_ui(Sub1_SPMFile)
    spm_fnirs_estimate_ui(Sub2_SPMFile)

  end
  
  
  
%% Extract SPM information
% wo sind unsere Daten?
srcPath = 'P:\projects\DCARE\DCARE\rawData\';


%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d']);
  end
  
  %%
  
    for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/HbO/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/HbO/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_SPMFile = strcat(Sub1SrcDir, sprintf(['SPM.mat']));
  Sub2_SPMFile = strcat(Sub2SrcDir, sprintf(['SPM.mat']));
      
  beta_id(i,1)=i;

    load(Sub1_SPMFile)
    beta_coop(i,1:16)=SPM.beta(2,:);
    beta_ind(i,1:16)=SPM.beta(3,:);
    beta_rest(i,1:16)=SPM.beta(4,:);
    beta_preform(i,1:16)=SPM.beta(6,:);
    
    clear SPM
    load(Sub2_SPMFile)
    beta_coop(i,17:32)=SPM.beta(2,:);
    beta_ind(i,17:32)=SPM.beta(3,:);
    beta_rest(i,17:32)=SPM.beta(4,:);
    beta_preform(i,17:32)=SPM.beta(6,:);
    
    end
  
%% WTC Analyse
% wo sind unsere Daten?
srcPath = 'P:\projects\DCARE\DCARE\rawData\';


%% welche Probanden sind in dem Ordner?
sourceList    = dir([srcPath, 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d']);
  end


 %% Wtc calculation per condition
  for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_NIRSFile = strcat(Sub1SrcDir, sprintf(['NIRS.mat']));
  Sub2_NIRSFile = strcat(Sub2SrcDir, sprintf(['NIRS.mat']));
  eventFile = load(strcat(Sub1SrcDir, sprintf(['multiple_conditions.mat'])));
  load(Sub1_NIRSFile)
  Y = reshape(spm_vec(rmfield(Y, 'od')), [P.ns P.nch 3]); 
  Y = spm_fnirs_preproc(Y, P); 
  Y = Y(:, :, 1);
  hbo_sub1 = spm_fnirs_filter(Y, P, 7.8125); 
  load(Sub2_NIRSFile)
  Y = reshape(spm_vec(rmfield(Y, 'od')), [P.ns P.nch 3]); 
  Y = spm_fnirs_preproc(Y, P); 
  Y = Y(:, :, 1);
  hbo_sub2 = spm_fnirs_filter(Y, P, 7.8125); 
  
  
% Calculate wavelet transform coherence, cross-correlations and control analyses
% 1) cut into trials
% events
fs=P.fs;
coop=eventFile.onsets{1,1}(1,:)*fs;
ind=eventFile.onsets{1,2}(1,:)*fs;
rest=eventFile.onsets{1,3}(1,:)*fs;
pf=eventFile.onsets{1,5}(1,:)*fs;

% format timepoints x channels x trial
hbo_sub1_coop(:,:,1) = hbo_sub1(coop(1,1):coop(1,1)+936,:);
hbo_sub1_coop(:,:,2) = hbo_sub1(coop(1,2):coop(1,2)+936,:);
hbo_sub1_ind(:,:,1) = hbo_sub1(ind(1,1):ind(1,1)+936,:);
hbo_sub1_ind(:,:,2) = hbo_sub1(ind(1,2):ind(1,2)+936,:);
hbo_sub1_rest(:,:,1) = hbo_sub1(rest(1,1):rest(1,1)+624,:);
hbo_sub1_rest(:,:,2) = hbo_sub1(rest(1,2):rest(1,2)+624,:);
hbo_sub1_rest(:,:,3) = hbo_sub1(rest(1,3):rest(1,3)+624,:);
hbo_sub1_pf(:,:) = hbo_sub1(pf(1,1):pf(1,1)+2343,:);

hbo_sub2_coop(:,:,1) = hbo_sub2(coop(1,1):coop(1,1)+936,:);
hbo_sub2_coop(:,:,2) = hbo_sub2(coop(1,2):coop(1,2)+936,:);
hbo_sub2_ind(:,:,1) = hbo_sub2(ind(1,1):ind(1,1)+936,:);
hbo_sub2_ind(:,:,2) = hbo_sub2(ind(1,2):ind(1,2)+936,:);
hbo_sub2_rest(:,:,1) = hbo_sub2(rest(1,1):rest(1,1)+624,:);
hbo_sub2_rest(:,:,2) = hbo_sub2(rest(1,2):rest(1,2)+624,:);
hbo_sub2_rest(:,:,3) = hbo_sub2(rest(1,3):rest(1,3)+624,:);
hbo_sub2_pf(:,:) = hbo_sub2(pf(1,1):pf(1,1)+2343,:);

% create t variable for WTC
t = 0:1/fs:length(hbo_sub1_coop)/fs - 1/fs;
t_coop = t';
t = 0:1/fs:length(hbo_sub1_ind)/fs - 1/fs;
t_ind = t';
t = 0:1/fs:length(hbo_sub1_rest)/fs - 1/fs;
t_rest = t';
t = 0:1/fs:length(hbo_sub1_pf)/fs - 1/fs;
t_pf = t';

% 2) calc. coupling
% cooperation
% find period
  sigPart1 = [t_coop, hbo_sub1_coop(:,1,1)];
  sigPart2 = [t_coop, hbo_sub2_coop(:,1,1)];
  try
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 10, 1, 'first');
  pnoi(2) = find(period < 50, 1, 'last');

for trial=1:2
for ch=1:16
    sigPart1 = [t_coop, hbo_sub1_coop(:,ch,trial)];
    sigPart2 = [t_coop, hbo_sub2_coop(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  

      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end

    
   coherence(ch, trial)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end

    catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
    end
% individual
% find period
  sigPart1 = [t_ind, hbo_sub1_ind(:,1,1)];
  sigPart2 = [t_ind, hbo_sub2_ind(:,1,1)];
  try
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 4, 1, 'first');
  pnoi(2) = find(period < 20, 1, 'last');

for trial=1:2
for ch=1:16
    sigPart1 = [t_ind, hbo_sub1_ind(:,ch,trial)];
    sigPart2 = [t_ind, hbo_sub2_ind(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  

      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end

    
   coherence(ch, trial,2)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end
catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
    end
% rest
% find period
  sigPart1 = [t_rest, hbo_sub1_rest(:,1,1)];
  sigPart2 = [t_rest, hbo_sub2_rest(:,1,1)];
  try
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 4, 1, 'first');
  pnoi(2) = find(period < 20, 1, 'last');

for trial=1:3
for ch=1:16
    sigPart1 = [t_rest, hbo_sub1_rest(:,ch,trial)];
    sigPart2 = [t_rest, hbo_sub2_rest(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end
    
   coherence(ch, trial,3)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end
catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
    end
% rest
% find period
  sigPart1 = [t_pf, hbo_sub1_pf(:,1,1)];
  sigPart2 = [t_pf, hbo_sub2_pf(:,1,1)];
  try
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 4, 1, 'first');
  pnoi(2) = find(period < 20, 1, 'last');

for trial=1
for ch=1:16
    sigPart1 = [t_pf, hbo_sub1_pf(:,ch,trial)];
    sigPart2 = [t_pf, hbo_sub2_pf(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end
    
   coherence(ch, trial,4)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end
catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
  end
  save(strcat(srcFolder, sprintf(['WTC.mat'])),'coherence');
  end
  

  
%% 3) random pair analysis and/or shuffling

rp_id = repmat(numOfPart,1,15);
rp_id_extra = rp_id(1,1:10);
rp_id(1,991:1000) = rp_id_extra;
rp_id = rp_id(randperm(numel(rp_id)));

%% 
  for i = 586:1:size(rp_id,2)
      
  srcFolder   = strcat(srcPath, sprintf(['DCARE_%02d/'], rp_id(1,i)));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/']));
  pos = randi(length(rp_id));
  card = rp_id(1,pos);
  srcFolder2   = strcat(srcPath, sprintf(['DCARE_%02d/'], card));
  Sub2SrcDir = strcat(srcFolder2, sprintf(['Subject2/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_NIRSFile = strcat(Sub1SrcDir, sprintf(['NIRS.mat']));
  Sub2_NIRSFile = strcat(Sub2SrcDir, sprintf(['NIRS.mat']));
  eventFile = load(strcat(Sub1SrcDir, sprintf(['multiple_conditions.mat'])));

  load(Sub1_NIRSFile)
  Y = reshape(spm_vec(rmfield(Y, 'od')), [P.ns P.nch 3]); 
  Y = spm_fnirs_preproc(Y, P); 
  Y = Y(:, :, 1);
  hbo_sub1 = spm_fnirs_filter(Y, P, 7.8125); 
  load(Sub2_NIRSFile)
  Y = reshape(spm_vec(rmfield(Y, 'od')), [P.ns P.nch 3]); 
  Y = spm_fnirs_preproc(Y, P); 
  Y = Y(:, :, 1);
  hbo_sub2 = spm_fnirs_filter(Y, P, 7.8125); 
  
  
% Calculate wavelet transform coherence, cross-correlations and control analyses
% 1) cut into trials
% events
fs=P.fs;
coop=eventFile.onsets{1,1}(1,:)*fs;
ind=eventFile.onsets{1,2}(1,:)*fs;
rest=eventFile.onsets{1,3}(1,:)*fs;
pf=eventFile.onsets{1,5}(1,:)*fs;

% format timepoints x channels x trial
hbo_sub1_coop(:,:,1) = hbo_sub1(coop(1,1):coop(1,1)+936,:);
hbo_sub1_coop(:,:,2) = hbo_sub1(coop(1,2):coop(1,2)+936,:);
hbo_sub1_ind(:,:,1) = hbo_sub1(ind(1,1):ind(1,1)+936,:);
hbo_sub1_ind(:,:,2) = hbo_sub1(ind(1,2):ind(1,2)+936,:);
hbo_sub1_rest(:,:,1) = hbo_sub1(rest(1,1):rest(1,1)+624,:);
hbo_sub1_rest(:,:,2) = hbo_sub1(rest(1,2):rest(1,2)+624,:);
hbo_sub1_rest(:,:,3) = hbo_sub1(rest(1,3):rest(1,3)+624,:);
hbo_sub1_pf(:,:) = hbo_sub1(pf(1,1):pf(1,1)+2343,:);

eventFile2 = load(strcat(Sub1SrcDir, sprintf(['multiple_conditions.mat'])));
coop=eventFile.onsets{1,1}(1,:)*fs;
ind=eventFile.onsets{1,2}(1,:)*fs;
rest=eventFile.onsets{1,3}(1,:)*fs;
pf=eventFile.onsets{1,5}(1,:)*fs;

hbo_sub2_coop(:,:,1) = hbo_sub2(coop(1,1):coop(1,1)+936,:);
hbo_sub2_coop(:,:,2) = hbo_sub2(coop(1,2):coop(1,2)+936,:);
hbo_sub2_ind(:,:,1) = hbo_sub2(ind(1,1):ind(1,1)+936,:);
hbo_sub2_ind(:,:,2) = hbo_sub2(ind(1,2):ind(1,2)+936,:);
hbo_sub2_rest(:,:,1) = hbo_sub2(rest(1,1):rest(1,1)+624,:);
hbo_sub2_rest(:,:,2) = hbo_sub2(rest(1,2):rest(1,2)+624,:);
hbo_sub2_rest(:,:,3) = hbo_sub2(rest(1,3):rest(1,3)+624,:);
try
hbo_sub2_pf(:,:) = hbo_sub2(pf(1,1):pf(1,1)+2343,:);
catch
    fprintf('Preschool Form skipped!')
end

% create t variable for WTC
t = 0:1/fs:length(hbo_sub1_coop)/fs - 1/fs;
t_coop = t';
t = 0:1/fs:length(hbo_sub1_ind)/fs - 1/fs;
t_ind = t';
t = 0:1/fs:length(hbo_sub1_rest)/fs - 1/fs;
t_rest = t';
t = 0:1/fs:length(hbo_sub1_pf)/fs - 1/fs;
t_pf = t';

% 2) calc. coupling
% cooperation
% find period
  sigPart1 = [t_coop, hbo_sub1_coop(:,1,1)];
  sigPart2 = [t_coop, hbo_sub2_coop(:,1,1)];
  try
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 10, 1, 'first');
  pnoi(2) = find(period < 50, 1, 'last');

for trial=1:2
for ch=1:16
    sigPart1 = [t_coop, hbo_sub1_coop(:,ch,trial)];
    sigPart2 = [t_coop, hbo_sub2_coop(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  

      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end

    
   coherence(ch, trial)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end

    catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
    end
% individual
% find period
  sigPart1 = [t_ind, hbo_sub1_ind(:,1,1)];
  sigPart2 = [t_ind, hbo_sub2_ind(:,1,1)];
  try
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 4, 1, 'first');
  pnoi(2) = find(period < 20, 1, 'last');

for trial=1:2
for ch=1:16
    sigPart1 = [t_ind, hbo_sub1_ind(:,ch,trial)];
    sigPart2 = [t_ind, hbo_sub2_ind(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  

      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end

    
   coherence(ch, trial,2)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end
catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
    end
% rest
% find period
  sigPart1 = [t_rest, hbo_sub1_rest(:,1,1)];
  sigPart2 = [t_rest, hbo_sub2_rest(:,1,1)];
  try
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 4, 1, 'first');
  pnoi(2) = find(period < 20, 1, 'last');

for trial=1:3
for ch=1:16
    sigPart1 = [t_rest, hbo_sub1_rest(:,ch,trial)];
    sigPart2 = [t_rest, hbo_sub2_rest(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end
    
   coherence(ch, trial,3)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end
catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
    end
% rest
% find period
try
  sigPart1 = [t_pf, hbo_sub1_pf(:,1,1)];
  sigPart2 = [t_pf, hbo_sub2_pf(:,1,1)];
  
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > 4, 1, 'first');
  pnoi(2) = find(period < 20, 1, 'last');

for trial=1
for ch=1:16
    sigPart1 = [t_pf, hbo_sub1_pf(:,ch,trial)];
    sigPart2 = [t_pf, hbo_sub2_pf(:,ch,trial)];
    [Rsq{ch, trial}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
      for j=1:1:length(coi)
        Rsq{ch, trial}(period >= coi(j), j) = NaN;
      end
    
   coherence(ch, trial,4)=nanmean(nanmean(Rsq{ch, trial}(pnoi(1):pnoi(2),:)));
end
end
catch
    fprintf('Inconsistent data in coherence calculation step for Dyad %d, skipped.\n', i);
  end
  save(sprintf(['P:/projects/DCARE/DCARE/MATLAB/procData/rpaData/rpa_%02d.mat'],i),'coherence');
  end
  
  
  

  
    
    