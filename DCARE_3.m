%DCARE Pipeline GLM
%%%%%%%%%%%%%%%

srcPath = 'C:\Users\Trinh Nguyen\Documents\MATLAB\DCARE\hmrData\';                        % raw data location
desPath = 'C:\Users\Trinh Nguyen\Documents\MATLAB\DCARE\glmData\';                  % processed data location
% gsePath = '\\fs.univie.ac.at\homedirs\nguyenq22\Documents\Projekte\InControl\procData\gen\';                    % path to CARE.SD


% prefix='CT';   
%% Scan for all subjects
if ~exist('numOfPart', 'var')                                               % estimate number of participants in raw data folder
  sourceList    = dir([srcPath, 'DCARE_*_sub1.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['DCARE_%d_sub1']);
  end
end


%%
for i = numOfPart
  fprintf('<strong>Dyad %d</strong>\n', i);
  
  % load preprocessed data
  srcFolder   = strcat(srcPath);
  filename    = sprintf(['DCARE_%02d_sub2'], i);
  
  fprintf('Load preprocessed data...\n');
  load([srcFolder filename]);
  
%%

% -------------------------------------------------------------------------
% Basic variables
% -------------------------------------------------------------------------
colRest     = 3;
colCoop     = 1;
colInd      = 2;
colPF       = 5;
    
fs=7.8125;

colAll            = colCoop|colInd|colRest|colPF;

durCoop  = round(120 * ...                % duration cooperation condition
                                  fs - 1);
durInd     = round(120 * ...                 % duration individual condition
                                  fs - 1);
durRest       = round(80 * ...                  % duration rest condition
                                  fs - 1);
durPF       = round(300 * ...                  % duration preschool form condition
                                  fs - 1);

% -------------------------------------------------------------------------
% Adapt the s matrix
% -------------------------------------------------------------------------
sMatrix = s;

evtCoop  = find(sMatrix(:, colCoop) > 0);
evtInd     = find(sMatrix(:, colInd) > 0);
evtRest       = find(sMatrix(:, colRest) > 0);
evtPF       = find(sMatrix(:, colPF) > 0);

for j = evtCoop'
    sMatrix(j:j+durCoop, colCoop) = 1; 
end

for j = evtInd'
    sMatrix(j:j+durInd, colInd) = 1;
end

for j = evtRest'
    sMatrix(j:j+durRest, colRest) = 1;
end

for j = evtPF'
    sMatrix(j:j+durPF, colPF) = 1;
end

% eventMarkers      = data_preproc.eventMarkers(colAll);
sMatrix_c           = sMatrix(:, colCoop);
% Convolution with hrf function
xBF = spm_get_bf(struct('dt',0.0018,'name','hrf'));
U.u=sMatrix_c;
U.name = {'reg'};
convreg_c = spm_Volterra(U, xBF.bf);

sMatrix_i           = sMatrix(:, colInd);
% Convolution with hrf function
xBF = spm_get_bf(struct('dt',0.0018,'name','hrf'));
U.u=sMatrix_i;
U.name = {'reg'};
convreg_i = spm_Volterra(U, xBF.bf);

sMatrix_r           = sMatrix(:, colRest);
% Convolution with hrf function
xBF = spm_get_bf(struct('dt',0.0018,'name','hrf'));
U.u=sMatrix_r;
U.name = {'reg'};
convreg_r = spm_Volterra(U, xBF.bf);

sMatrix_pf           = sMatrix(:, colPF);
% Convolution with hrf function
xBF = spm_get_bf(struct('dt',0.0018,'name','hrf'));
U.u=sMatrix_pf;
U.name = {'reg'};
convreg_pf = spm_Volterra(U, xBF.bf);

convreg=[convreg_c convreg_i convreg_r convreg_pf];  

% -------------------------------------------------------------------------
% Adapt the s matrix
% -------------------------------------------------------------------------
fprintf('<strong>Conduct generalized linear model regression for all channels...</strong>\n');
% y           =dc;
% sMatrix           =sMatrix(:, colAll);
% Aaux        =aux;
% trange      =[0 120];
% glmSolveMethod = 1;
% idxBasis    =2;
% paramsBasis =[0.1 3.0 10.0 1.8 3.0 10.0];
% rhoSD_ssThresh =0;
% driftOrder=3;
% flagSSmethod=[];
% flagMotionCorrect=0;
% 
% [yavg, yavgstd, tHRF, nTrials, ynew, yresid, ysum2, beta, yR] = hmrDeconvHRF_DriftSS(y, s, t, SD, Aaux, tIncAuto,...
%     trange, glmSolveMethod, idxBasis, paramsBasis, rhoSD_ssThresh, flagSSmethod, driftOrder, flagMotionCorrect );
% 
% data_glm.yavg=yavg;
% data_glm.yavgstd=yavgstd;
% data_glm.tHRF=tHRF;
% data_glm.nTrials=nTrials;
% data_glm.ynew=ynew;
% data_glm.yresid=yresid;
% data_glm.ysum2=ysum2;
% data_glm.beta=beta;
% data_glm.yR=yR;
data.hbo=hbo;
data.fs=fs;

data_glm = execGLM(convreg, data);

%% save beta values of glm regression
  cfg             = [];
  desFolder   = strcat(desPath);
  filename    = sprintf(['CARE_%02d_glm_sub2'],i);
  
  file_path = strcat(desFolder, filename, ...
                     '.mat');

  fprintf('The generalized linear model coefficients of dyad %d will be saved in:\n', i); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'data_glm');
  fprintf('Data stored!\n\n');
  clear data_glm data_preproc 
end

%% GLM Function definition
function data_out = execGLM(s, data_in)
    % build output matrix
    beta = zeros(16,5);                                 
  
  for channel = 1:1:size(data_in.hbo, 2)
    % conduct generalized linear model regression
    % beta estimates for a generalized linear regression of the responses 
    % in data_in.hbo(:, channel) on the predictors in the sMatrix
    if ~isnan(data_in.hbo(1, channel))                                      % check if channel was not rejected during preprocessing
      beta(channel,:) = glmfit(s, data_in.hbo(:, channel));
    else
      beta(channel,:) = NaN;
    end
    
    
  end
  
  % put results into a structure
%   data_out.eventMarkers = evtMark;
  data_out.s            = s;
  data_out.hbo          = data_in.hbo;
  data_out.time         = (1:1:size(data_in.hbo, 1)) / data_in.fs;
  data_out.fsample      = data_in.fs;
  data_out.channel      = 1:1:size(data_in.hbo, 2);
  data_out.beta         = beta(:, 2:end);                                    % for the existing conditions only the columns 2:end are relevant  
fprintf('<strong>Generalized linear model was computed and saved...</strong>\n');


end