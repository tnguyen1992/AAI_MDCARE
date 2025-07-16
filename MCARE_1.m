% %%% DCARE Pipeline
% convert raw data to homer2 nirs format

% wo sind unsere Daten?
srcPath = 'P:\teamdata\hoehl\projects\DCARE\MCARE\MCARE_fNIRS_rawData\';

% wohin sollen die Daten gespeichert werden?
desPath = 'P:\teamdata\hoehl\projects\DCARE\MCARE\MCARE_fNIRS_procData\01_nirs\';

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
  
  Sub1DesFile = strcat(desPath, sprintf(['DCARE_%02d_sub1'], ...
                      i), '.nirs');
  Sub2DesFile = strcat(desPath, sprintf(['DCARE_%02d_sub2'], ...
                      i), '.nirs');

% -------------------------------------------------------------------------
% Load SD file
% -------------------------------------------------------------------------
SDfile = 'P:\projects\DCARE\DCARE\MATLAB\CARE.SD';
load(SDfile, '-mat', 'SD');

% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_wl1File = strcat(Sub1SrcDir, sprintf(['DCARE_%02d'], i), '.wl1');
  Sub1_wl2File = strcat(Sub1SrcDir, sprintf(['DCARE_%02d'], i), '.wl2');
  Sub1_hdrFile = strcat(Sub1SrcDir, sprintf(['DCARE_%02d'], i), '.hdr');
 
  Sub2_wl1File = strcat(Sub2SrcDir, sprintf(['DCARE_%02d'], i), '.wl1');
  Sub2_wl2File = strcat(Sub2SrcDir, sprintf(['DCARE_%02d'], i), '.wl2');
  Sub2_hdrFile = strcat(Sub2SrcDir, sprintf(['DCARE_%02d'], i), '.hdr');
  

% -------------------------------------------------------------------------
% Convert and export data
% -------------------------------------------------------------------------
convertData(Sub1DesFile, Sub1_wl1File, Sub1_wl2File, Sub1_hdrFile, SD,...
            i);
convertData(Sub2DesFile, Sub2_wl1File, Sub2_wl2File, Sub2_hdrFile, SD,...
            i);

end

% -------------------------------------------------------------------------
% SUBFUNCTION data convertion
% -------------------------------------------------------------------------
function convertData (desFile, wl1File, wl2File, hdrFile, SD, num)
wl1 = load(wl1File);                                                        % load .wl1 file
wl2 = load(wl2File);                                                        % load .wl2 file

d = [wl1 wl2];                                                              % d matrix from .wl1 and .wl2 files

fid = fopen(hdrFile);
tmp = textscan(fid,'%s','delimiter','\n');                                  % this just reads every line
hdr_str = tmp{1};
fclose(fid);

keyword = 'Sources=';                                                       % find number of sources
tmp = hdr_str{strncmp(hdr_str, keyword, length(keyword))};
NIRxSources = str2double(tmp(length(keyword)+1:end));

keyword = 'Detectors=';                                                     % find number of detectors
tmp = hdr_str{strncmp(hdr_str, keyword, length(keyword))};
NIRxDetectors = str2double(tmp(length(keyword)+1:end));

if NIRxSources < SD.nSrcs || NIRxDetectors < SD.nDets                       % Compare number of sources and detectors to SD file
   error('The number of sources and detectors in the NIRx files does not match your SD file...');
end

keyword = 'SamplingRate=';                                                  % find Sample rate
tmp = hdr_str{strncmp(hdr_str, keyword, 13)};
fs = str2double(tmp(length(keyword)+1:end));

% find Active Source-Detector pairs
keyword = 'S-D-Mask="#';
ind = find(strncmp(hdr_str, keyword, length(keyword))) + 1;
ind2 = find(strncmp(hdr_str(ind+1:end), '#', 1)) - 1;
ind2 = ind + ind2(1);
sd_ind = cell2mat(cellfun(@str2num, hdr_str(ind:ind2), 'UniformOutput', 0));
sd_ind = sd_ind';
sd_ind = logical([sd_ind(:);sd_ind(:)]);
d = d(:, sd_ind);

% find NaN values in the recorded data -> channels should be pruned as 'bad'
for i=1:size(d,2)
    if nonzeros(isnan(d(:,i)))
        SD.MeasListAct(i) = 0;
    end
end

% find event markers and build s vector
keyword = 'Events="#';
ind = find(strncmp(hdr_str, keyword, length(keyword))) + 1;
ind2 = find(strncmp(hdr_str(ind+1:end), '#', 1)) - 1;
ind2 = ind + ind2(1);
events = cell2mat(cellfun(@str2num, hdr_str(ind:ind2), 'UniformOutput', 0));
events = events(:,2:3);
markertypes = unique(events(:,1));
s = zeros(length(d),length(markertypes));
for i = 1:length(markertypes)
    s(events(events(:,1) == markertypes(i), 2), i) = 1;
end

% create t, aux varibles
aux = ones(length(d),1);                                                    %#ok<NASGU>
t = 0:1/fs:length(d)/fs - 1/fs;
t = t';                                                                     %#ok<NASGU>

fprintf('Saving NIRS file: %s...\n', desFile);
save(desFile, 'd', 's', 't', 'aux', 'SD');
fprintf('Data stored!\n\n');

end
  
