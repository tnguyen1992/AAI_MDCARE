%% CARE Hyperscanning Pipeline
% Cross-correlation analysis

% -------------------------------------------------------------------------
% Path settings
% -------------------------------------------------------------------------
desPath = 'D:/InvTie/';      % memory space for processed data
addpath('C:\Users\tnguyen\OneDrive - Fondazione Istituto Italiano Tecnologia\Documenti\MATLAB\toolboxes\wavelet-coherence-master')
% -------------------------------------------------------------------------
% Dyad list
% -------------------------------------------------------------------------
sourceList    = dir([strcat(desPath), ...
    strcat('*02a_preproc*', '.mat')]);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
numOfSources  = length(sourceList);
numOfPart     = zeros(1, numOfSources);

for i=1:1:numOfSources
    %     numOfPart(i)  = sscanf(sourceList{i}, ...
    %         strcat( '_d%d_02a_preproc*', ...
    %         '.mat'));
    numOfPart(i) = str2double(sourceList{i}(8:9));
end



%% -------------------------------------------------------------------------
% RPA
% -------------------------------------------------------------------------
poi = [10 50];
numOfPermutations = 139; %(N-1)

for i = 1:1:numOfSources
    for n=1:1:numOfPermutations

        fprintf('<strong>Dyad %s</strong>\n', sourceList{1, i});
        fprintf('<strong>Permutation %d</strong>\n', n);

        load(sourceList{1, i});

        hboSub1 = data_preproc.sub1.hbo;

        % Load a different file for each permutation
        permutationFileIndex = mod(n, length(sourceList)) + 1;  % Change this line based on your requirements
        load(sourceList{1, permutationFileIndex});
        hboSub2 = data_preproc.sub2.hbo;

        % Adjust both files to the shortest length
        minLen = min(length(hboSub1), length(hboSub2));
        hboSub1 = hboSub1(1:minLen,1:16);
        hboSub2 = hboSub2(1:minLen,1:16);

        t = (0:(1/data_preproc.sub1.fs):((size(hboSub1, 1) - 1) / ...
            data_preproc.sub1.fs))';

        numOfChan = size(hboSub1, 2);

        % -------------------------------------------------------------------------
        % Basic variable
        % Determine events
        % -------------------------------------------------------------------------

        durCollaboration  = round(120 * ...                % duration collaboration condition
            data_preproc.sub1.fs - 1);
        durIndividual     = round(120 * ...                 % duration individual condition
            data_preproc.sub1.fs - 1);
        durBaseline       = round(80 * ...                  % duration baseline condition
            data_preproc.sub1.fs - 1);

        sMatrix = data_preproc.sub1.s;

        evtCollaboration  = find(sMatrix(:, 1) > 0);
        evtIndividual     = find(sMatrix(:, 2) > 0);
        evtBaseline       = find(sMatrix(:, 3) > 0);
        % -------------------------------------------------------------------------
        % Estimate periods of interest
        % -------------------------------------------------------------------------
        pnoi = zeros(2,1);
        ii = numOfChan;
        while (isnan(hboSub1(1, ii)) || isnan(hboSub2(1, ii)))                        % check if 16 th channel was not rejected in both subjects during preprocessing
            ii = ii - 1;                                                                % if 16th channel was rejected in at least on subject
            if ii == 0                                                                 % search for next channel which was not rejected
                break;
            end
        end
        if ii ~= 0
            sigPart1 = [t, hboSub1(:,ii)];
            sigPart2 = [t, hboSub2(:,ii)];
            [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0, 'AR1', 'auto');
            pnoi(1) = find(period > poi(1), 1, 'first');
            pnoi(2) = find(period < poi(2), 1, 'last');
        else
            period = NaN;                                                             % if all channel were rejected, the value period cannot be extimated and will be therefore set to NaN
        end

        % -------------------------------------------------------------------------
        % Allocate memory
        % -------------------------------------------------------------------------
        Rsq{numOfChan} = [];
        Rsq(:) = {NaN(length(period), length(t))};

        % -------------------------------------------------------------------------
        % Calculate Coherence increase between conditions for every channel of the
        % dyad
        % -------------------------------------------------------------------------
        for ii=1:1:numOfChan
            if ~isnan(hboSub1(1, ii)) && ~isnan(hboSub2(1, ii))                         % check if this channel was not rejected in both subjects during preprocessing
                sigPart1 = [t, hboSub1(:,ii)];
                sigPart2 = [t, hboSub2(:,ii)];
                [Rsq{ii}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence

                for j=1:1:length(coi)
                    Rsq{ii}(period >= coi(j), j) = NaN;
                end



                % calculate mean activation in frequency band of interest
                % collaboration condition
                for j=1:1:length(evtCollaboration)
                    try
                        meanCohCollab(j)  = nanmean(nanmean(Rsq{ii}(pnoi(1):pnoi(2), ...
                            evtCollaboration(j):evtCollaboration(j) + ...
                            durCollaboration)));
                    catch
                        meanCohCollab(j)  = nanmean(nanmean(Rsq{ii}(pnoi(1):pnoi(2), ...
                            evtCollaboration(j):end)));
                    end

                end

                % individual condition
                for j=1:1:length(evtIndividual)
                    try
                        meanCohIndiv(j)   = nanmean(nanmean(Rsq{ii}(pnoi(1):pnoi(2), ...
                            evtIndividual(j):evtIndividual(j) + ...
                            durIndividual)));
                    catch
                        meanCohIndiv(j)   = nanmean(nanmean(Rsq{ii}(pnoi(1):pnoi(2), ...
                            evtIndividual(j):end)));
                    end

                end

                % baseline
                for j=1:1:length(evtBaseline)
                    try
                        meanCohBase(j)    = nanmean(nanmean(Rsq{ii}(pnoi(1):pnoi(2), ...
                            evtBaseline(j):evtBaseline(j) + ...
                            durBaseline)));
                    catch
                        meanCohBase(j)    = nanmean(nanmean(Rsq{ii}(pnoi(1):pnoi(2), ...
                            evtBaseline(j):end)));
                    end
                end




            end

            try
                collaboration(ii)  = nanmean(meanCohCollab);
                % average mean coherences over trials
            catch
                collaboration(ii) = NaN;
            end
            try
                individual(ii)     = nanmean(meanCohIndiv);
            catch
                individual(ii) = NaN;
            end
            try
                baseline(ii)       = nanmean(meanCohBase);
            catch
                baseline(ii)       = NaN;
            end
        end
        coherences(1:3, 1:16,n) = [collaboration; individual; baseline];

    end

    coherences_m=nanmean(coherences,3);

    % put results into the output data structure
    data_wtc.hbo.coherence_m       = coherences_m;
    data_wtc.hbo.coherence       = coherences;
    %     data_wtc.Rsq                  = Rsq;
    data_wtc.paramStrings         = {'Cooperation', 'Individual', 'Rest'};
    data_wtc.channel              = 1:1:size(hboSub1, 2);
    data_wtc.t                    = t;
    data_wtc.cfg.period           = period;
    data_wtc.cfg.poi              = poi;

    % save  values

    % Extract the filename without the extension
    [~, filename, ~] = fileparts(sourceList{1, i});

    % Remove the "-mat" suffix
    filenameWithoutMat = strrep(filename, '-mat', '');

    cfg             = [];
    cfg.desFolder   = strcat(desPath, 'RPA/');
    cfg.filename    = strcat(filenameWithoutMat,'_rpa.mat');

    file_path = strcat(cfg.desFolder, cfg.filename,  ...
        '.mat');


    save(file_path,'data_wtc');
    fprintf('Data stored!\n\n');
    clear data_wtc data_preproc coherences coherences_m meanCohCollab meanCohBase meanCohIndiv Rsq

end




%% %% Extract rpa for revision 2
clear all
clc


%% Data extraction
sourceList    = dir(['D:/InvTie/RPA/', 'MCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['MCARE_d%d_02a_*']);
  end

%% WTC 
i=1;
for id = numOfPart

  RPAFile = strcat('D:/InvTie/RPA/', sprintf(['MCARE_d%02d_02a_preproc_002_rpa.mat.mat'],id));
  load(RPAFile)
  

	for condition=1:1:3
        for channel=1:1:16
                rpa(i,1)=id;
                rpa(i,2)=condition;
                rpa(i,3)=channel;
                rpa(i,4)=data_wtc.hbo.coherence_m(condition, channel);
                i=i+1;
        end
    end

end

%% Data extraction
sourceList    = dir(['D:/InvTie/RPA/', 'DCARE_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ['DCARE_d%d_02a_*']);
  end

%% WTC 
for id = numOfPart

  RPAFile = strcat('D:/InvTie/RPA/', sprintf(['DCARE_d%02d_02a_preproc_001_rpa.mat.mat'],id));
  load(RPAFile)
  

	for condition=1:1:3
        for channel=1:1:16
                rpa(i,1)=id+90;
                rpa(i,2)=condition;
                rpa(i,3)=channel;
                rpa(i,4)=data_wtc.hbo.coherence_m(condition, channel);
                i=i+1;
        end
    end

end

dlmwrite('MDCARE_AAI_rpa.csv', rpa)