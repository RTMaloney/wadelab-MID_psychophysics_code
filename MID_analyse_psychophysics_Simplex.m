% Generate a mean set of data across all 4 subjects, and model the parameters in the same way as individual subjects.
% NOTE: this function loads the PROCESSED data files ie the mean parameters for each condition (not the raw values).
% Remember that in all matrices, frequencies are given in ROWS and amplitudes in COLUMNS.
% When plotting the matrices we transpose the matrix to put frequency across the x-axis.

% *** Remember that the order of the matrices (within the cell arrays) is CD coherence, CD contrast
% IOVD coherence, IOVD contrast. ***

% In this version we use the Simplex algorithm (fminsearch) to perform the modelling
% (not the backslash operator).

%% Load the data

% Each subject code, to load in turn.
SubjCode = {'AS', 'JA', 'MF', 'RM'};
% We will DOWNSAMPLE the data of RM and JA so it is in the same space as everyone else.
% This will allow us to average across the 4 subjects before doing any modelling or upsampling.

for ii = 1:4 %loop across the 4 subjects.
    
    % Load the processed FULL cue data, to fit the model to:
    FULLcue{ii} = load (fullfile('data', [SubjCode{ii} '_Proc_thresh_FULLCUE_data.mat']));
    % Load the processed CD & IOVD data:
    MIDcues{ii} = load (fullfile('data', [SubjCode{ii} '_Proc_thresh_data.mat']));
    
    % Outliers in  JA data
    % Need to fix the 2 outliers in JA's IOVD coherence data before we do anything.
    % They are in MeanSensData{3}(3,3) or (21) and MeanSensData{3}(end)
    % We will do a linear interpolation of these values using the matlab function 'filloutliers'.
    % (Previously we used the median of all other values).
    if strcmp(SubjCode{ii}, 'JA')
        MIDcues{ii}.MeanSensData{3} = filloutliers(MIDcues{ii}.MeanSensData{3}, 'linear');
        %OutlierIdx = true(length(MIDcues{ii}.MeanSensData{3}(:)),1);
        %OutlierIdx(21) = false;
        %OutlierIdx(end) = false;
        %MIDcues{ii}.MeanSensData{3}(~OutlierIdx) = nan;
        %MIDcues{ii}.MeanSensData{3} = fillmissing(MIDcues{ii}.MeanSensData{3}, 'linear');
        %MIDcues{ii}.MeanSensData{3}(~OutlierIdx) = median(MIDcues{ii}.MeanSensData{3}(OutlierIdx));
    end
    
    % Transpose the data.
    % At this point we will transpose all the data matrices because we later plot frequency along the x-axis
    % Thus we want frequency to be along the x-axis (put it running along the COLUMNS)
    % and amplitude to be arranged along the ROWS.
    % This will avoid confusion.
    FULLcue{ii}.MeanSensData = cellfun(@transpose, FULLcue{ii}.MeanSensData, 'UniformOutput', false);
    CD_IOVD{ii}.Sens = cellfun(@transpose, MIDcues{ii}.MeanSensData, 'UniformOutput', false);
    CD_IOVD{ii}.Betas = cellfun(@transpose, MIDcues{ii}.MeanBetaData, 'UniformOutput', false);
    %CD_IOVD{ii}.ThreshErr = cellfun(@transpose, MIDcues{ii}.MeanThreshData_err, 'UniformOutput', false);
    
    % NOW, subjects RM and JA have sampled the CD & IOVD space more widely (9 freq * 7 ampl) than AS and MF (5 * 5).
    % FULL cue was always sampled at 5 frequencies * 5 amplitudes.
    % Here we will downsample their data to fit in the same space as MF and AS.
    
    % Establish the parameters we want in the end:
    desired_frequency = round(logspace(log10(0.5), log10(8), 5) * 100) / 100; % Same for all conditions
    desired_amplitudes_CD = round(logspace(log10(1.67), log10(16.67),5) * 100) / 100;
    desired_amplitudes_IOVD = round(logspace(log10(1.67), log10(167),5) * 100) / 100;
    
    % For RM and JA, we have data matrices of 9 rows (freq) * 7 columns (ampl).
    if strcmp(SubjCode{ii}, 'RM') || strcmp(SubjCode{ii}, 'JA')
        
        sourceFreq = round(logspace(log10(0.5), log10(8), 9) * 100) / 100;
        % For CD, we work out the source location
        sourceCDAmp = round(logspace(log10(1.67), log10(16.67),7) * 100) / 100;
        % and for IOVD, we also work out the source location
        sourceIOVDAmp = round(logspace(log10(1.67), log10(167),7) * 100) / 100;
        
        % We want to downsample/resample the data, so we set up the grids needed for the interpolation here.
        % Remember that frequency is the same for all.
        [d_fr, d_am_CD] = meshgrid(desired_frequency, desired_amplitudes_CD);  % grid of the desired parameters: CD
        [~, d_am_IOVD] = meshgrid(desired_frequency, desired_amplitudes_IOVD); % grid of the desired parameters: IOVD
        [~, CD_source_am] = meshgrid (sourceFreq, sourceCDAmp);                % grid of original CD parameters
        [source_fr, IOVD_source_am] = meshgrid (sourceFreq, sourceIOVDAmp);    % grid of original IOVD parameters
        
        % Now loop across the 2 threshold types: coherence and contrast.
        % Remember that resulting order of matrices are: CD coh, CD cont, (1&2) IOVD coh, IOVD cont (3&4).
        for signal = 1:2
            % Downsample the betas and sensitivities (inverse alphas) here.
            CD_IOVD{ii}.Betas{signal} = interp2(source_fr, CD_source_am, CD_IOVD{ii}.Betas{signal}, d_fr, d_am_CD);         %Downsampled betas for CD
            CD_IOVD{ii}.Betas{signal+2} = interp2(source_fr, IOVD_source_am, CD_IOVD{ii}.Betas{signal+2}, d_fr, d_am_IOVD); %Downsampled betas for IOVD
            CD_IOVD{ii}.Sens{signal} = interp2(source_fr, CD_source_am, CD_IOVD{ii}.Sens{signal}, d_fr, d_am_CD);           %Downsampled alphas for CD
            CD_IOVD{ii}.Sens{signal+2} = interp2(source_fr, IOVD_source_am, CD_IOVD{ii}.Sens{signal+2}, d_fr, d_am_IOVD);   %Downsampled alphas for IOVD
        end
    end
    
end % end of loop across the 4 subjects.

% So now we have the data matrices downsampled (where appropriate) and sorted, all in the same space
% for each of the 4 observers.
% We need to generate an average matrix now across these 4 observers and do the same modelling as performed on individual subjects.
% The averages are done in an inefficient and dirty way perhaps; BUT at least they avoid
% any indexing confusion where the conditions get mixed up.

%% Compute the mean matrices across subjects:  FULL cue data

tmp = [];
tmp2 = [];
for ii=1:4 %loop across subjects
    tmp(:,:,ii) = FULLcue{ii}.MeanSensData{1}; % for coherence
    tmp2(:,:,ii) = FULLcue{ii}.MeanSensData{2}; % for contrast
end
% Compute the mean matrix:
%(this works for a mean matrix and is the same as doing [1+2+3+4]/4).
FULLcue_av.MeanSensData{1} = mean(tmp,3);  % Cell {1} = coherence
FULLcue_av.MeanSensData{2} = mean(tmp2,3); % Cell {2} = contrast

%% Now, average the sensitivities:
% Need to loop across subjects (4), but also across conditions (4)
tmp = [];
tmp2 = [];
tmp3 = [];
tmp4 = [];
for subject = 1:4 % loop across the 4 subjects within each condition.
    tmp(:,:,subject) = CD_IOVD{subject}.Sens{1}; % For CD, coh
    tmp2(:,:,subject) = CD_IOVD{subject}.Sens{2}; % For CD, cont
    tmp3(:,:,subject) = CD_IOVD{subject}.Sens{3}; % For IOVD, coh
    tmp4(:,:,subject) = CD_IOVD{subject}.Sens{4}; % For IOVD, cont
end

% Now compute mean sensitivity across subjects for each condition.
CD_IOVD_av.Sens{1} = mean(tmp,3); % CD, coh
CD_IOVD_av.Sens{2} = mean(tmp2,3); % CD, cont
CD_IOVD_av.Sens{3} = mean(tmp3,3); % IOVD, coh
CD_IOVD_av.Sens{4} = mean(tmp4,3); % IOVD, cont

%% Now, do the same for the betas:
tmp = [];
tmp2 = [];
tmp3 = [];
tmp4 = [];
for subject = 1:4 % loop across the 4 subjects within each condition.
    tmp(:,:,subject) = CD_IOVD{subject}.Betas{1}; % For CD, coh
    tmp2(:,:,subject) = CD_IOVD{subject}.Betas{2}; % For CD, cont
    tmp3(:,:,subject) = CD_IOVD{subject}.Betas{3}; % For IOVD, coh
    tmp4(:,:,subject) = CD_IOVD{subject}.Betas{4}; % For IOVD, cont
end

% Now compute mean betas across subjects for each condition.
CD_IOVD_av.Betas{1} = mean(tmp,3);  % CD, coh
CD_IOVD_av.Betas{2} = mean(tmp2,3); % CD, cont
CD_IOVD_av.Betas{3} = mean(tmp3,3); % IOVD, coh
CD_IOVD_av.Betas{4} = mean(tmp4,3); % IOVD, cont

%%
% So now we have our mean matrices across subjects: full cue (coh/cont) CD (coh/cont) and IOVD (coh/cont)
% Remember that they have been transposed: arranged so that frequency is along the x-axis (put it running along the COLUMNS)
% and amplitude is along the y-axis (ROWS).

% We can perform the modelling, but first we need to resample the CD data (both averages and single subject)
% to be in same space as IOVD and FULL.

%% Resample CD data to be in same amplitude space as IOVD and FULL.
% Here, we will first add ones to the space where there's no data for CD.
% This is because we can assume there is zero sensitivity for CD at ranges > 16.67 arcmin.
% Then we can DOWNSAMPLE the data; hopefully giving us a smoother estimate.

% Remember that now, the desired and current source frequencies are the same.
sourceFreq = desired_frequency;
% Same goes for the CD source amplitudes.
sourceCDAmp = [desired_amplitudes_CD, 35.98, 77.51, 167]; % Add the hypothetical amplitudes to CD
new_CD_desired_amp = desired_amplitudes_IOVD; % same existing amplitude range as FULL and IOVD

% We want to resample the data, so we set up the grids needed for the interpolation here.
% Remember that frequency is the same for all.
[d_fr, d_am_CD] = meshgrid(desired_frequency, new_CD_desired_amp);  % grid of the desired parameters: CD
[CD_new_src_freq, CD_source_am] = meshgrid (sourceFreq, sourceCDAmp);             % grid of original CD parameters

%%
% Now loop across the 4 subjects and the 2 threshold types: coherence and contrast.
% Remember that resulting order of matrices are: CD coh, CD cont, (1&2)
for subject = 1:4
    for signal = 1:2
        
        % First, add the 1s to the empty space:
        CD_IOVD{subject}.Betas{signal} = [CD_IOVD{subject}.Betas{signal}; ones(3,5)]; %(those go on the bottom; later we the matrix upside down for plotting).
        CD_IOVD{subject}.Sens{signal} = [CD_IOVD{subject}.Sens{signal}; ones(3,5)]; 
        
        % Downsample the betas and sensitivities (inverse alphas) here.
        CD_IOVD{subject}.Betas{signal} = interp2(CD_new_src_freq, CD_source_am, CD_IOVD{subject}.Betas{signal}, d_fr, d_am_CD); % Downsampled betas for CD
        CD_IOVD{subject}.Sens{signal} = interp2(CD_new_src_freq, CD_source_am, CD_IOVD{subject}.Sens{signal}, d_fr, d_am_CD);   % Downsampled inverse alphas for CD
        
    end
end

%%
% Now do the same for the average (CD coh and cont only):
for signal = 1:2
    
    % Add the one's to the hypothetical space now:
    CD_IOVD_av.Betas{signal} = [CD_IOVD_av.Betas{signal}; ones(3,5)];
    CD_IOVD_av.Sens{signal} = [CD_IOVD_av.Sens{signal}; ones(3,5)];
    
    % Now perform the interpolation:
    CD_IOVD_av.Betas{signal} = interp2(CD_new_src_freq, CD_source_am, CD_IOVD_av.Betas{signal}, d_fr, d_am_CD);
    CD_IOVD_av.Sens{signal} = interp2(CD_new_src_freq, CD_source_am, CD_IOVD_av.Sens{signal}, d_fr, d_am_CD);
    
end

% We should now have everything we need to perform the modelling.
%% Perform the model fits: average data

% Default settings for the Simplex:
nfuncevals = 20000;
defopts = optimset ('fminsearch');
options = optimset (defopts, 'Display', 'iter', 'MaxFunEvals', nfuncevals, 'MaxIter', 8000, 'PlotFcns', @optimplotfval);

%%
for signal = 1:2 % loop across coherence and contrast
    
    % Initial guess of parameters
    maxparams = [realmax realmax]';  %Maximum allowable parameters
    params = [rand rand]'; %numbers for initial guess of starting parameters (replace these with best fitting paramters if initial search fails to minimize)
    
    data = FULLcue_av.MeanSensData{signal}(:);                              % FULL cue data, to be fit
    inv_alpha = [CD_IOVD_av.Sens{signal}(:), CD_IOVD_av.Sens{signal+2}(:)]; % Inverse alphas, CD and IOVD
    betas = [CD_IOVD_av.Betas{signal}(:), CD_IOVD_av.Betas{signal+2}(:)];   % Betas, CD and IOVD
    
    % Find the coefficients (the CD & IOVD weights).
    % This is model 1, without the betas.
    [w_m1_av_data(:, signal), SumSquaredError, dummy, output] = fminsearch (@(params) findBF_MID_simple_model(params, maxparams, data, inv_alpha), params, options);

    % Now solve model 2, which includes the beta values.
    [w_m2_av_data(:, signal), SumSquaredError, dummy, output] = fminsearch (@(params) findBF_MID(params, maxparams, data, betas, inv_alpha), params, options);

    % Now the moment of truth: Is the RMS error for the beta-weighted data
    % different to that of the non-beta-weighted data? Does the visual system
    % weight by SNR?
    
    % Set aside the model predictions computed according to the fitted weights.
    % We multiply (and sum) the inverse alphas for each cue by their respective weight through matrix multiplication.
    predictions(:,1) = [CD_IOVD_av.Sens{signal}(:), CD_IOVD_av.Sens{signal+2}(:)]*w_m1_av_data(:,signal); % Model predictions according to model 1
    predictions(:,2) = [CD_IOVD_av.Sens{signal}(:).*CD_IOVD_av.Betas{signal}(:), ...
        CD_IOVD_av.Sens{signal+2}(:).*CD_IOVD_av.Betas{signal+2}(:)]*w_m2_av_data(:,signal);              % Model predictions according to model 2
    
    % Now we could do a chi^2 goodness of fit between the data and the predictions
    % NMSE costs vary between -Inf (bad fit) to 1 (perfect fit). If the cost function is equal to zero, then x is no better than a straight line at matching xref.
    GofFit_av_data(signal,1) = goodnessOfFit(predictions(:,1), data, 'NMSE'); % Col 1 = Model 1, coh, cont
    GofFit_av_data(signal,2) = goodnessOfFit(predictions(:,2), data, 'NMSE'); % Col 2 = Model 2, coh, cont
    
end

%% Perform model fits: each subject
% Here we will model the data for each individual subject, just as done with the average above.
GofFit = [];
for subject = 1:4 % loop across subjects
    for signal = 1:2 % loop across coherence and contrast
        
        % Initial guess of parameters
        maxparams = [realmax realmax]';  %Maximum allowable parameters
        params = [rand rand]'; %numbers for initial guess of starting parameters (replace these with best fitting paramters if initial search fails to minimize)
        
        data = FULLcue{subject}.MeanSensData{signal}(:);                              % FULL cue data, to be fit
        inv_alpha = [CD_IOVD{subject}.Sens{signal}(:), CD_IOVD{subject}.Sens{signal+2}(:)]; % Inverse alphas, CD and IOVD
        betas = [CD_IOVD{subject}.Betas{signal}(:), CD_IOVD{subject}.Betas{signal+2}(:)];   % Betas, CD and IOVD
        
        % Find the coefficients (the CD & IOVD weights).
        % This is model 1, without the betas.
        [w_m1(:, signal, subject), SumSquaredError, dummy, output] = fminsearch (@(params) findBF_MID_simple_model(params, maxparams, data, inv_alpha), params, options);
        
        % Now solve model 2, which includes the beta values.
        [w_m2(:, signal, subject), SumSquaredError, dummy, output] = fminsearch (@(params) findBF_MID(params, maxparams, data, betas, inv_alpha), params, options);
        
        % Now the moment of truth: Is the RMS error for the beta-weighted data
        % different to that of the non-beta-weighted data? Does the visual system
        % weight by SNR?
        
        % Set aside the model predictions computed according to the fitted weights.
        % We multiply (and sum) the inverse alphas for each cue by their respective weight through matrix multiplication.
        predictions(:,1) = [CD_IOVD{subject}.Sens{signal}(:), CD_IOVD{subject}.Sens{signal+2}(:)]*w_m1(:,signal,subject); % Model predictions according to model 1
        predictions(:,2) = [CD_IOVD{subject}.Sens{signal}(:) .* CD_IOVD{subject}.Betas{signal}(:), ...
            CD_IOVD{subject}.Sens{signal+2}(:) .* CD_IOVD{subject}.Betas{signal+2}(:)]*w_m2(:,signal,subject);              % Model predictions according to model 2
        
        % Now we could do a chi^2 goodness of fit between the data and the predictions
        % NMSE costs vary between -Inf (bad fit) to 1 (perfect fit). If the cost function is equal to zero, then x is no better than a straight line at matching xref.
        GofFit(signal,1,subject) = goodnessOfFit(predictions(:,1), data, 'NMSE'); % Col 1 = Model 1, coh, cont
        GofFit(signal,2,subject) = goodnessOfFit(predictions(:,2), data, 'NMSE'); % Col 1 = Model 1, coh, cont
        
    end
end

% Now that we have the model coefficients, we can use these to plot the data.
% These will be interpolated to a wider parameter space (for presentation only)
% and the coefficients used on the interpolated data to display predictions.
% But this will be done in a different script: MID_psychophysics_plot_matrices.m

%%
% Save the modelling results, as well as the averages across subjects
save('MID_psychophysics_SIMPLEX_modelling_results.mat', 'GofFit_av_data', 'w_m1_av_data', 'w_m2_av_data', ...
    'GofFit', 'w_m1', 'w_m2', 'FULLcue_av', 'CD_IOVD_av');
