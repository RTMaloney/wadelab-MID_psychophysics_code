% MID_psychophysics_plot_matrices.m

% Loads and interpolates MID psychophysics data for display only.
% Uses the parameters established with the Simplex mode of model fitting to
% display the predicted data for the FULL cue.
% The matrices are for the manuscript describing this work.

% Load modelling coefficients and average matrices across the 4 observers:
load 'MID_psychophysics_SIMPLEX_modelling_results.mat'

%% Set up the fixed variables:
frequency = round(logspace(log10(0.5), log10(8), 9) * 100) / 100;

CDamplitudes = round(logspace(log10(1.67), log10(16.67),7) * 100) / 100;
IOVDamplitudes = round(logspace(log10(1.67), log10(167),7) * 100) / 100;

%Let's just fix a few of the overlapping values so they are the same in the 2 cues:
IOVDamplitudes(4) = 16.67;  %was 16.7
CDamplitudes(5) = 7.75; % was 7.74

% Set up labels for the Y-axis.
% We want them doubled to represent 'binocular' disparity, and slightly adjusted to look better:
BinocAmpl = num2cell(round(fliplr(IOVDamplitudes * 2))); % double for binoc disparity/displacement
BinocAmpl{1} = 333;
BinocAmpl{4} = 33.3;
BinocAmpl{end-1} = 7.2;
BinocAmpl{end} = 3.33;

%% Plot matrix of velocities for methods
% First, let's make a plot of the MAX velocities tested,
% as displayed in the ECVP 2016 poster:

% Make a matrix of the velocities so we can do a contour plot
% to illustrate the different velocities (for the methods).
% NOTE: we are calculating the MONOCULAR velocities:
a = logspace(log10(1.67), log10(167), 20) / 60; %Div by 60 to convert to degrees
f = logspace(log10(0.5), log10(8), 18);
Max_vel = [];
for ii = 1:length(a)
    for jj = 1: length(f)
        
        % Find velocity of the hot spot:
        prd = 1./f(jj); %in sec
        angFreq = 2 * pi * f(jj); % remember that a sine wave is always some multiple of 2*pi, because it is a periodic function.
        RefreshRate = 120; %absolute screen refresh rate
        %the effective frame rate PER EYE: since each eye is stimulated on successive video frames
        PeyeFR = RefreshRate/2; %Per eye f.r. always half the absolute f.r.
        %Determine the shift in PHASE of the wave, if any. =0 if no shift
        PhaseShiftCycles = 0; %fraction of a cycle you want to shift it by, eg 1/4, 1/3, 1/2 etc etc
        PhaseShift = PhaseShiftCycles * 2 * pi; %as with angular frequency, the phase is always some multiple of 2*pi (ie it is specified in radians).
        %The no. of frames for a full cycle
        FrmsFullCyc = round(PeyeFR*prd); %must have integer no. of frames
        
        %Length of the sine wave:
        t = linspace(0, prd, FrmsFullCyc); %One complete cycle, IN SEC, in steps of per-eye frames
        
        % Velocity is calculated according to the equation given here:
        % web('http://labman.phys.utk.edu/phys135/modules/m9/oscillations.htm')
        %Compute the velocity:
        velocity = -angFreq * a(ii) * sin(angFreq * t + PhaseShift);
        
        %Take the max velocity. Note that the min velocity should be the same value, only negative.
        % Divide by 60 to get the results in deg/s rather than arcmin/s
        Max_vel(ii,jj) = max(velocity);
    end
end
figure(1)
contourf(Max_vel,linspace(0, max(max(Max_vel)),25))
colorbar('Ticks',[1, 25, 50, 75, 100, 125], ...
    'TickDirection', 'out');
set(gca, 'YTickLabel', fliplr(BinocAmpl),...
    'YTick', linspace(1,20,7), ...
    'YScale', 'linear', ...
    'TickDir', 'out', 'Box', 'off', ...
    'XTickLabel',num2cell(frequency), ...
    'XTick', linspace(1,18,9), ...
    'XScale', 'linear', ...
    'FontSize',14);

%Set rendering and save the velocity figure as eps:
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(1), 'Figures/2019 matrices/MID_methods_velocities.eps','epsc');
close

%% Let's plot the AVERAGE data matrices, for CD and IOVD.
% So there are 4 plots here: CD coherence, CD contrast, IOVD coh, IOVD cont.
% IMPORTANT: the data are stored *in that order*, so when we loop across signal (ie coh or cont)
% we add 2 to the signal (1,2) to give us the IOVD indices (3,4)

% Average sensitivity data are in CD_IOVD_av.Sens{}, and are each 5*5 matrices.
% We also interpolate the beta values, to use them (multiplied by the model coefficients)
% when estimating the FULL cue results based on CD and IOVD.

% NOTE: in all plots, we need to set 0.9 as the min value plotted in the matrix (though they are never <1.0)
% This is to get the colorbar to fit in the right place to have the NaNs values at the end.

%Set the colormap:
cm = colormap; close 
cm = [0 0 0; cm];

% Set color ranges for matrices: defines colorbars
GridRange = [0.9 3; 0.9 6]; % smaller range for Coherence compared to Contrast
% Using the same range for all conditions loses the effect in CD. (colorbar is too narrow)

% Interpolate to 18 since that is a multiple of 9, the no. of frequencies tested,
InterpFreq = logspace(log10(0.5),log10(8),18);
% If we interpolate to an odd number, 19, we get the more sensible logarithmic range (easier for plotting)
InterpAmpl = logspace(log10(1.67), log10(167), 20);

sourceFreq = round(logspace(log10(0.5), log10(8), 5) * 100) / 100;
% For CD, we work out the source location
sourceCDAmp = round(logspace(log10(1.67), log10(16.67),5) * 100) / 100;
% and for IOVD, we also work out the source location
sourceIOVDAmp = round(logspace(log10(1.67), log10(167),5) * 100) / 100;

% We want to resample the data, so we set up the grids needed for the interpolation here.
% Remember that frequency is the same for all.
[d_fr, d_am_CD] = meshgrid(InterpFreq, InterpAmpl);                 % grid of the desired parameters: CD
[~, d_am_IOVD] = meshgrid(InterpFreq, InterpAmpl);                  % grid of the desired parameters: IOVD
[~, CD_source_am] = meshgrid (sourceFreq, sourceCDAmp);             % grid of original CD parameters
[source_fr, IOVD_source_am] = meshgrid (sourceFreq, sourceIOVDAmp); % grid of original IOVD parameters

%%
for signal = 1:2
    
    % ------------------- %
    %     Do CD first
    % ------------------- %
    Sens = CD_IOVD_av.Sens{signal}; % CD sensitivity
    fgrid_CD{signal} = interp2(source_fr, CD_source_am, Sens, d_fr, d_am_CD);
    
    Betas = CD_IOVD_av.Betas{signal}; % CD betas.
    fgrid_CD_betas{signal} = interp2(source_fr, CD_source_am, Betas, d_fr, d_am_CD);
    
    figure
    %imagesc(flipud(fgrid_CD{signal}))
    imagesc(flipud(fgrid_CD{signal}), GridRange(signal,:))
    colormap(cm)
    colorbar('TickDirection', 'out')
    %     colorbar('Ticks',0.5:0.5:3, ...
    %         'TickLabels',{'','1','','2','','3'}, ...
    %         'TickDirection', 'out')
    
    % Set up the axes. This is tricky, because the variables are logarithmically spaced,
    % but the matrix is 'linear', so we specify the values logarithmically, but trick the scale into being linear
    % Things get confusing because we have 10 unique amplitudes that are NOT spaced logarithmically all together (though they are for
    % an individual cue eg CD or IOVD).
    % This is why we use the IOVD amplitudes to label the y-axis on both cues.
    %set(gca, 'YTickLabel', num2cell(fliplr(IOVDamplitudes)),...
    set(gca, 'YTickLabel', BinocAmpl,...
        'YTick', linspace(1,20,7), ...
        'YScale', 'linear', ...
        'TickDir', 'out', 'Box', 'off', ...
        'XTickLabel',num2cell(frequency), ...
        'XTick', linspace(1,18,9), ...
        'XScale', 'linear', ...
        'FontSize',14);
    %'XTick',linspace(0.5,20,9)); %,...
    %'FontSize',14)
    %'XTickLabel',num2cell([0 frequency]),...
    %     title([SubjCode, ' IOVD, coherence sensitivity'])
    %     ylabel('Max lateral displacement (arcmin)')
    %     xlabel('Temporal frequency (Hz)')
    
    %plot contour lines over data:
    %hold on
    %contour(flipud(fgrid_CD{signal}),5,'Color','w','LineStyle', '--', 'LineWidth', 0.3)
    
    % ------------------- %
    %     Now do IOVD
    % ------------------- %
    Sens = CD_IOVD_av.Sens{signal+2}; % IOVD sensitivity
    fgrid_IOVD{signal} = interp2(source_fr, IOVD_source_am, Sens, d_fr, d_am_IOVD);
    
    Betas = CD_IOVD_av.Betas{signal+2}; % IOVD betas
    fgrid_IOVD_betas{signal} = interp2(source_fr, IOVD_source_am, Betas, d_fr, d_am_IOVD);
    
    figure
    %imagesc(flipud(fgrid_IOVD{signal}))
    imagesc(flipud(fgrid_IOVD{signal}), GridRange(signal,:))
    colormap(cm)
    colorbar('TickDirection', 'out')
    %     colorbar('Ticks',0.5:0.5:7, ...
    %         'TickLabels',{'NA','1','','2','','3','','4','','5','','6','','7'}, ...
    %         'TickDirection', 'out')
    
    set(gca, 'YTickLabel', BinocAmpl,...
        'YTick', linspace(1,20,7), ...
        'YScale', 'linear', ...
        'TickDir', 'out', 'Box', 'off', ...
        'XTickLabel',num2cell(frequency), ...
        'XTick', linspace(1,18,9), ...
        'XScale', 'linear', ...
        'FontSize',14);
    
    %plot contour lines over data:
    %hold on
    %contour(flipud(fgrid_IOVD{signal}),5,'Color','w','LineStyle', '--', 'LineWidth', 0.3)
    
end

%Set rendering and save all the figures as eps:
figure(1)
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(1), 'Figures/2019 matrices/mean_matrix_CD_coh.eps','epsc');

figure(2)
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(2), 'Figures/2019 matrices/mean_matrix_IOVD_coh.eps','epsc');

figure(3)
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(3), 'Figures/2019 matrices/mean_matrix_CD_cont.eps','epsc');

figure(4)
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(4), 'Figures/2019 matrices/mean_matrix_IOVD_cont.eps','epsc');


%% Now plot the FULL cue average matrices.

% All of the desired and source frequencies and amplitudes for FULL are the same as above;
% except that the range of amplitudes for FULL is the same as IOVD: so we just re-use that here.

% Set color ranges for matrices: defines colorbars
GridRange = [0.9 5; 0.9 6]; % smaller range for Coherence compared to Contrast

% We want to resample the data, so we set up the grids needed for the interpolation here.
% Remember that frequency is the same for all.
[d_fr, d_am_FULL] = meshgrid(InterpFreq, InterpAmpl);               % grid of the desired parameters
[source_fr, FULL_source_am] = meshgrid (sourceFreq, sourceIOVDAmp); % grid of original parameters

for signal = 1:2
    
    Sens = FULLcue_av.MeanSensData{signal}; % FULL cue sensitivity
    fgrid_FULL{signal} = interp2(source_fr, FULL_source_am, Sens, d_fr, d_am_FULL);
    
    figure
    %imagesc(flipud(fgrid_FULL{signal}))
    imagesc(flipud(fgrid_FULL{signal}), GridRange(signal,:))
    colormap(cm)
    colorbar('TickDirection', 'out')
    % colorbar('Ticks',0.5:0.5:6, 'TickLabels',{'NA','1','','2','','3','','4','','5', '', '6'},'TickDirection', 'out')
    
    % Set up the axes. This is tricky, because the variables are logarithmically spaced,
    % but the matrix is 'linear', so we specify the values logarithmically, but trick the scale into being linear
    % Things get confusing because we have 10 unique amplitudes that are NOT spaced logarithmically all together (though they are for
    % an individual cue eg CD or IOVD).
    % This is why we use the IOVD amplitudes to label the y-axis on both cues.
    %set(gca, 'YTickLabel', num2cell(fliplr(IOVDamplitudes)),...
    set(gca, 'YTickLabel', BinocAmpl,...
        'YTick', linspace(1,20,7), ...
        'YScale', 'linear', ...
        'TickDir', 'out', 'Box', 'off', ...
        'XTickLabel',num2cell(frequency), ...
        'XTick', linspace(1,18,9), ...
        'XScale', 'linear', ...
        'FontSize',14);
    
    %plot contour lines over data:
    %hold on
    %contour(flipud(fgrid_FULL{signal}),5,'Color','w','LineStyle', '--', 'LineWidth', 0.3)    
end

% Set rendering and save the figures:
figure(5)
colorbar('Ticks',0.5:0.5:5, 'TickLabels',{'','1','','2','','3','','4','','5'},'TickDirection', 'out')
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(5), 'Figures/2019 matrices/mean_matrix_FULL_coh.eps','epsc');

figure(6)
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(6), 'Figures/2019 matrices/mean_matrix_FULL_cont.eps','epsc');

%% Plot estimated FULL cue data

% Here we want to take the model coefficients established earlier, and use those as well as
% the mean CD and IOVD data to generate predicted FULL cue matrices.
% These are the results of MODEL 2, which includes the betas as weights.

nansumdim = 2; % Dimension along which to sum the CD and IOVD estimates in nansum
GridRange = [0.9 5; 0.9 6];
for signal = 1:2
    
    %estimate the full cue data using the betas, inverse alphas and the weights from the model:
    % ie: CD alpha * CD beta * CD weight + IOVD alpha * IOVD beta * IOVD weight.
    % Need to use nansum here to ensure IOVD correctly added in space where no CD data exists
    fgrid_est{signal} = nansum([fgrid_CD{signal}(:).*fgrid_CD_betas{signal}(:) * w_m2_av_data(1,signal), ...
        fgrid_IOVD{signal}(:).*fgrid_IOVD_betas{signal}(:) * w_m2_av_data(2,signal)], nansumdim);
    
    % reshape the results back into a grid:
    fgrid_est{signal} = reshape(fgrid_est{signal}, 20,18);
    
    figure
    %imagesc(flipud(fgrid_est{signal}))
    imagesc(flipud(fgrid_est{signal}), GridRange(signal,:))
    colormap(cm)
    colorbar('TickDirection', 'out')
    % colorbar('Ticks',0.5:0.5:6, 'TickLabels',{'NA','1','','2','','3','','4','','5', '', '6'},'TickDirection', 'out')
    
    set(gca, 'YTickLabel', BinocAmpl,...
        'YTick', linspace(1,20,7), ...
        'YScale', 'linear', ...
        'TickDir', 'out', 'Box', 'off', ...
        'XTickLabel',num2cell(frequency), ...
        'XTick', linspace(1,18,9), ...
        'XScale', 'linear', ...
        'FontSize',14);
    
    %plot contour lines over data:
    %hold on
    %contour(flipud(fgrid_est{signal}),5,'Color','w','LineStyle', '--', 'LineWidth', 0.3) 
end

% Set rendering and save the figures:
figure(7)
colorbar('Ticks',0.5:0.5:5, 'TickLabels',{'','1','','2','','3','','4','','5'},'TickDirection', 'out')
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(7), 'Figures/2019 matrices/mean_matrix_ESTIM_coh.eps','epsc');

figure(8)
set(gcf,'renderer','painters')  % painters seems the best renderer for saving eps images
saveas(figure(8), 'Figures/2019 matrices/mean_matrix_ESTIM_cont.eps','epsc');

%% Identify peak performance velocities:
% --------------------------------- %
% find the hot spots for both cues:
% --------------------------------- %
% (in case we need them)
f=[]; a=[];
for signal=1:2
    
    % CD hot spots:
    [r, c] = find(fgrid_CD{signal}==max(max(fgrid_CD{signal})));
    %Rows are amplitude:
    a(signal,1) = InterpAmpl(r);
    f(signal,1) = InterpFreq(c);
    
    % IOVD hot spots:
    [r, c] = find(fgrid_IOVD{signal}==max(max(fgrid_IOVD{signal})));
    %Rows are amplitude:
    a(signal,2) = InterpAmpl(r);
    f(signal,2) = InterpFreq(c);
    
    %Find velocity of the hot spot:
    prd(signal,:) = 1./f(signal,:); %in sec
    angFreq(signal,:) = 2 * pi * f(signal,:); % remember that a sine wave is always some multiple of 2*pi, because it is a periodic function.
    RefreshRate = 120; %absolute screen refresh rate
    %the effective frame rate PER EYE: since each eye is stimulated on successive video frames
    PeyeFR = RefreshRate/2; %Per eye f.r. always half the absolute f.r.
    %Determine the shift in PHASE of the wave, if any. =0 if no shift
    PhaseShiftCycles = 0; %fraction of a cycle you want to shift it by, eg 1/4, 1/3, 1/2 etc etc
    PhaseShift = PhaseShiftCycles * 2 * pi; %as with angular frequency, the phase is always some multiple of 2*pi (ie it is specified in radians).
    %The no. of frames for a full cycle
    FrmsFullCyc(signal,:) = round(PeyeFR*prd(signal,:)); %must have integer no. of frames
    
    %Length of the sine wave:
    t_CD = linspace(0, prd(signal,2), FrmsFullCyc(signal,1)); %One complete cycle, IN SEC, in steps of per-eye frames
    t_IOVD = linspace(0, prd(signal,2), FrmsFullCyc(signal,2)); %One complete cycle, IN SEC, in steps of per-eye frames
    
    %Compute the velocity:
    velocity_CD = -angFreq(signal,1) * a(signal,1) * sin(angFreq(signal,1) * t_CD + PhaseShift);
    velocity_IOVD = -angFreq(signal,2) * a(signal,2) * sin(angFreq(signal,2) * t_IOVD + PhaseShift);
    
    %Take the max velocity. Note that the min velocity should be the same value, only negative.
    %Divide by 60 to get the results in deg/s rather than arcmin/s
    Hot_max_vels(signal,1) = max(velocity_CD)/60
    Hot_max_vels(signal,2) = max(velocity_IOVD)/60
end


