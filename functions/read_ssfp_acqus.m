% Reads Bruker's acqus file and extracts experimental parameters
% relevant for SSFP reconstruction. Syntax:
%
%              parameters=read_ssfp_acqus(directory)
%
% Parameters:
%
%    directory - path to the directory containing 
%                the acqus file
%
% Outputs the following fields of the parameters structure:
%
%    .isotope      - working isotope, character string
%
%    .carrier_frq  - carrier frequency for the
%                    working isotope, Hz
%
%    .offset_frq   - transmitter offset, Hz
%
%    .sp_width     - spectral width, Hz
%
%    .flip_angle   - pulse flip angle, degrees
%
%    .n_phase_incs - number of SSFP phase increments
%
%    .n_acqps      - number of complex points ac-
%                    quired in each SFFP readout
%
%    .n_scans      - number of SSFP signal avera-
%                    ging scans as [dummy acquired]
%
%    .pulse_dur    - pulse duration, seconds
%
%    .dead_pp      - pre-pulse dead time, seconds
%
%    .dead_ap      - after-pulse dead time, seconds
%
% mark.shif@weizmann.ac.il
% ilya.kuprov@weizmann.ac.il

function parameters=read_ssfp_acqus(directory)

% Read acquisition parameters file
file_name=[directory filesep 'acqus'];
disp(['Bruker parameter file:      ' file_name]);
acqus=readlines(file_name);

% Isotope
isotope=acqus(startsWith(acqus,'##$NUC1= '));
if isempty(isotope)
    error('working isotope not found.');
else
    parameters.isotope=char(extractBetween(isotope,'<','>'));
    disp(['Working isotope:            ' parameters.isotope]);
end

% Carrier frequency
carrier_frq_string=acqus(startsWith(acqus,'##$SFO1='));
if contains(carrier_frq_string,'.')
    parameters.carrier_frq=1e6*str2double(regexp(carrier_frq_string,'\d+(\.\d*)','match'));
else
    parameters.carrier_frq=1e6*str2double(regexp(carrier_frq_string,'\d+','match'));
end
if isempty(parameters.carrier_frq)
    error('carrier frequency not found.');
else
    disp(['Carrier frequency:          ' num2str(parameters.carrier_frq) ' Hz']);
end

% Transmitter offset
offset_freq_string=acqus(startsWith(acqus,'##$O1='));
if contains(offset_freq_string,'.')
    parameters.offset_frq=str2double(regexp(offset_freq_string,'\d+(\.\d*)','match'));
else
    parameters.offset_frq=str2double(regexp(offset_freq_string,'\d+','match'));
end
if isempty(parameters.offset_frq)
    error('transmitter offset not found.');
else
    disp(['Transmitter offset:         ' num2str(parameters.offset_frq) ' Hz']);
end

% Spectral width
sp_width_string=acqus(startsWith(acqus,'##$SW_h='));
if contains(sp_width_string,'.')
    parameters.sp_width=str2double(regexp(sp_width_string,'\d+(\.\d*)','match'));
else
    parameters.sp_width=str2double(regexp(sp_width_string,'\d+','match'));
end
if isempty(parameters.sp_width)
    error('spectral width not found.');
else
    disp(['Spectral width:             ' num2str(parameters.sp_width) ' Hz']);
end

% Pulse flip angle
constants_index=find(acqus=='##$CNST= (0..63)');
if isempty(constants_index)||(~isscalar(constants_index))
    error('flip angle not found.');
else
    constants=split(acqus(constants_index + 1));
    parameters.flip_angle=str2double(constants(11));
    disp(['Pulse flip angle:           ' num2str(parameters.flip_angle) ' degrees']);
end

% Pre-pulse dead time
delays_index=find(acqus=='##$D= (0..63)');
if isempty(delays_index)||(~isscalar(delays_index))
    error('relaxation delay not found.');
else
    delay_lengths=split(acqus(delays_index + 1));
    parameters.dead_pp=str2double(delay_lengths(2));
    disp(['Pre-pulse dead time:        ' num2str(parameters.dead_pp) ' seconds']);
end

% Pulse duration
pulse_index=find(acqus=='##$P= (0..63)');
if isempty(pulse_index)||(~isscalar(pulse_index))
    error('pulse duration not found.');
else
    pulse_lengths=split(acqus(pulse_index + 1));
    parameters.pulse_dur=1e-6*str2double(pulse_lengths(4));
    disp(['Pulse duration:             ' num2str(parameters.pulse_dur) ' seconds']);
end

% After-pulse dead time
dead_ap=acqus(startsWith(acqus,'##$DE='));
parameters.dead_ap=1e-6*str2double(regexp(dead_ap,'\d+','match'));
if isempty(parameters.dead_ap)
    error('after-pulse dead time not found.');
else
    disp(['After-pulse dead time:      ' num2str(parameters.dead_ap) ' seconds']);
end

% Number of dummy scans
d_scans_string=acqus(startsWith(acqus,'##$DS='));
d_scans=str2double(regexp(d_scans_string,'\d+','match'));
if isempty(d_scans)
    error('number of dummy scans not found.');
else
    parameters.n_scans(1)=d_scans;
    disp(['Number of dummy scans:      ' int2str(parameters.n_scans(1))]);
end

% Number of acquired scans
a_scans_string=acqus(startsWith(acqus,'##$NS='));
a_scans=str2double(regexp(a_scans_string,'\d+','match'));
if isempty(a_scans)
    error('number of acquired scans not found.');
else
    parameters.n_scans(2)=a_scans;
    disp(['Number of acquired scans:   ' int2str(parameters.n_scans(2))]);
end

% Number of complex points
n_acqps_string=acqus(startsWith(acqus,'##$TD='));
parameters.n_acqps=str2double(regexp(n_acqps_string,'\d+','match'))/2;
if isempty(parameters.n_acqps)
    error('number of complex points not found.');
else
    disp(['Number of complex points:   ' int2str(parameters.n_acqps)]);
end

% Reading TR from the acqus file
D_idx = find(acqus == '##$D= (0..63)');
delays = split(acqus(D_idx + 1));
parameters.rep_time = str2double(delays(11));

if isempty(parameters.rep_time)
    error('repetition time not found.');
else
    disp(['repetition time:   ' num2str(parameters.rep_time) ' seconds']);
end

end

% My work always tried to unite the true with the beau-
% tiful; but when I had to choose one or the other, I 
% usually chose the beautiful.
%
% Hermann Weyl

