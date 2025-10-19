% Reads Bruker's ser or fid file and extracts the FID
% Parameters:
%
%    directory - path to the directory containing 
%                the ser or fid file
%
% Output:
%    
%    fid - a vector array containing the FID for the experiment 
%
%

function [fid] = read_bruker_data(directory)

filenames = {'fid','ser'};
    file_ID = -1;
    for i = 1:numel(filenames)
        path = fullfile(directory, filenames{i});
        if exist(path,'file')
            file_ID = fopen(path,'r');
            break;
        end
    end

    if file_ID == -1
        error('fid or ser file not found in: %s', directory);
    end

data = fread(file_ID,'int32');
fclose(file_ID);
fid = data(1:2:end) + 1i*data(2:2:end);

end