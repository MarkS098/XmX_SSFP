clc; close all; clearvars;

% Filepath and directory numbers
data_dir = '/home/mark/NMR Data/XmX processing/LF_Acetylacetone_30_11_25_DMSO_NaOH_293K'; % main experiment directory
data_mat_name = 'peak_B_acac'; % File name for the .mat file

% Processing parameters
save_data_mat = true;
plot_spect = true;
normalize = false;
file_type = 'fid'; % options: abs, fid, 1r, 1i
dir_num = [111,10:109]; % experiment folder names, first file is reference
mean_val = zeros(1,numel(dir_num)); % pre-allocating mean FID value array
TR_vals = zeros(1,numel(dir_num)); % pre-allocating TR array
SNR = zeros(1,numel(dir_num));

n_skip = 4; % number of points to remove from the start of the FID
bounds = [20,30]; % maxima calculation boundaries
noise_bound = [120,160]; % in case a peak is mixed with the noise we take the mean value of the noise in its expected region
min_height = 4e8; % minimum intensity for a point in the spectrum to be considered a peak
min_seperation = 4; % minimum seperation for points to be considered peaks

for j = 1:numel(dir_num)
    
    % Reading the data
    directory = [data_dir, filesep, num2str(dir_num(j))];
    switch file_type
        case '1r'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1r']);
            data = fread(file_ID,'int32');
            fclose(file_ID);
            
        case '1i'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1i']);
            data = fread(file_ID,'int32');
            fclose(file_ID);
        
        case 'fid' 
            data = read_bruker_data(directory);
            data = data(n_skip+1:end);
            spect = abs(fftshift(fft(data)));
         
        case 'abs'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1i']);
            data_imag = fread(file_ID,'int32');
            fclose(file_ID);

            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1r']);
            data_real = fread(file_ID,'int32');
            fclose(file_ID);

            data = data_real + 1i*data_imag;
            data = abs(data);

    end


    % Read acquisition and processing parameters
    parameters = read_ssfp_acqus(directory);
    proc_parameters = read_ssfp_procs(directory);

    % Acquisition parameters
    TD = parameters.n_acqps;
    SFO1 = parameters.carrier_frq;
    O1 = parameters.offset_frq; O1p = O1/(SFO1*1e-6);
    SWH = parameters.sp_width;
    TR_vals(j) = parameters.rep_time*1e3;

    % Processing parameters
    SI = proc_parameters.SI;
    NC_proc = proc_parameters.nc_proc;
    n_points = TD;

    if (strcmp(file_type,'1r') || strcmp(file_type,'1i')|| strcmp(file_type,'abs'))
        n_skip = 0;
        n_points = SI;
        
        % Getting the spectrum via FFT and scaling
        spect = flip(data)/(2^NC_proc);
    end

    % Calculating frequency axis
    Hz_axis = linspace(-SWH/2,SWH/2,n_points - n_skip) + O1; % We remove the corresponding number of points from the axis for plotting
    ppm_axis = Hz_axis/(SFO1*1e-6); % Converting the frequency axis to ppm

    % Taking a section from the spectrum for maxima calculations
    sub_axis = ppm_axis(ppm_axis > min(bounds) & ppm_axis < max(bounds));
    sub_spect = spect(ppm_axis > min(bounds) & ppm_axis < max(bounds));
    
    % Determining the maxima of the spectrum
    [pks,locs] = findpeaks(sub_spect,sub_axis,"MinPeakHeight",min_height,"MinPeakDistance",min_seperation);

    if j==1
        num_peaks = numel(pks);
    end
    

    % Sort the peak intensities into arrays for plotting
    % if the off-resonance peak disappears take the mean of the noise in
    if numel(locs) < 2
        peaks(j) = mean(spect(ppm_axis > min(noise_bound) & ppm_axis < max(noise_bound))); % off-resonance

    else
        peaks(j) = spect(find(ppm_axis == min(locs))); % off-resonance
    end

    % Calculating the SNR for both peaks
    noise_std = std(noise_bound, 1);
    noise_avg = mean(noise_bound);

    SNR(j) = (peaks(j) - noise_avg)/noise_std; 

    if plot_spect == true
        % Plotting the spectrum with removed points for every dataset
        figure()
        hold on
        title(['TR = ',num2str(TR_vals(j)),' ms'])
        plot(ppm_axis,spect)
        xlabel('\delta^{13} (ppm)')
        set(gca,'XDir','reverse')
    end
end

% Plotting on-resonance vs off-resonance peaks intensity values
ref_val = max([peaks]); % maximum value for normalization
peaks(peaks == ref_val) = [];
SNR(1) = [];
TR_vals(1) = [];

if normalize == true
    peaks = peaks/ref_val;
end


figure()
hold on
plot(TR_vals, peaks,'--ro','LineWidth',1)
xlabel('TR (ms)')
ylabel('Normalized Intensity')
legend('on-resonance')

% Saving the .mat file in the data directory
if save_data_mat == true
    save([data_dir,filesep,data_mat_name],'peaks','TR_vals')
end