clc; close all; clearvars;

% Filepath and directory numbers
data_dir = '/home/mark/NMR Data/XmX processing/XmX data/LF_DimethylAcetamide_2W_11_09_25'; % main experiment directory
data_mat_name = 'LF_DimethylAcetamide_31_08_25 25K-1r'; % File name for the .mat file

% Processing parameters
save_data_mat = true;
plot_spect = false;
file_type = 'fid'; % options: fid, 1r, 1i
dir_num = [300,195:207]; % experiment folder names, first file is reference
mean_val = zeros(1,numel(dir_num)); % pre-allocating mean FID value array
TR_vals = zeros(1,numel(dir_num)); % pre-allocating TR array
SNR_A = zeros(1,numel(dir_num));
SNR_B = zeros(1,numel(dir_num));

n_skip = 4; % number of points to remove from the start of the FID
bounds = [33,38]; % maxima calculation boundaries
noise_bound = [40,45]; % in case a peak is mixed with the noise we take the mean value of the noise in its expected region
min_height = 100; % minimum intensity for a point in the spectrum to be considered a peak
min_seperation = 2.5; % minimum seperation for points to be considered peaks

for j = 1:numel(dir_num)
    
    % Reading the data
    directory = [data_dir, filesep, num2str(dir_num(j))];
    switch file_type
        case '1r'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1r']);
            Data = fread(file_ID,'int32');
            fclose(file_ID);
            
        case '1i'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1i']);
            Data = fread(file_ID,'int32');
            fclose(file_ID);
        
        case 'fid' 
            Data = read_bruker_data(directory);
            Data = Data(n_skip+1:end);
            spect = abs(fftshift(fft(Data)));
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

    if (strcmp(file_type,'1r') || strcmp(file_type,'1i'))
        n_skip = 0;
        n_points = SI;
        
        % Getting the spectrum via FFT and scaling
        spect = flip(Data)/(2^NC_proc);
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
        peak_B(j) = mean(spect(ppm_axis > min(noise_bound) & ppm_axis < max(noise_bound))); % off-resonance

    else
        peak_B(j) = spect(find(ppm_axis == min(locs))); % off-resonance
    end

    peak_A(j) = spect(find(ppm_axis == max(locs))); % on resonance

    % Calculating the SNR for both peaks
    noise_std = std(noise_bound, 1);
    noise_avg = mean(noise_bound);

    SNR_A(j) = (peak_A(j) - noise_avg)/noise_std; 
    SNR_B(j) = (peak_B(j) - noise_avg)/noise_std; 

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
ref_val = max([peak_B]); % maximum value for normalization
peak_B(peak_B == ref_val) = [];
SNR_B(1) = [];
TR_vals(1) = [];

figure()
hold on
plot(TR_vals, peak_B/ref_val,'--ro','LineWidth',1)
xlabel('TR (ms)')
ylabel('Normalized Intensity')
legend('on-resonance')

figure()
hold on
plot(TR_vals, SNR_B,'--o','LineWidth',1)
xlabel('TR (ms)')
ylabel('SNR')
legend('on-resonance')

% Saving the .mat file in the data directory
if save_data_mat == true
    save([data_dir,filesep,data_mat_name],'peak_A','peak_B','TR_vals')
end