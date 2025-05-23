function VPmonkeysessionList = VPmonkey_importSessionlist %(filename, dataLines)
%IMPORTFILE Import data from a text file
%  VPMONKEYSESSIONLIST = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  VPMONKEYSESSIONLIST = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  VPmonkeysessionList = importfile("/home/rick/neuro/iannettilab/Projects/VP monkey/VPmonkey_sessionList.csv", [4, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 27-Jun-2022 14:25:18

%% Input handling

filename = [getRoot '/VPmonkey/VPmonkey_sessionList.csv'];

% If dataLines is not specified, define defaults
% if nargin < 2
    dataLines = [4, Inf];
% end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 21);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["session", "subjects", "date", "headfixed", "starttime", "stoptime", "video_setup", "video_recording", "vis_stimulus", "vis_offset", "vis_photocell", "vis_comments", "som_intensity", "som_offset", "som_twitch", "aud_ampgain", "aud_offset", "aud_beep_attenuation", "aud_noise_attenuation", "music_experiment", "comments"];
opts.VariableTypes = ["string", "categorical", "datetime", "char", "datetime", "datetime", "categorical", "categorical", "categorical", "double", "double", "char", "char", "double", "char", "double", "double", "double", "double", "double", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["headfixed", "vis_comments", "som_intensity", "som_twitch"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["subjects", "headfixed", "video_setup", "video_recording", "vis_stimulus", "vis_comments", "som_intensity", "som_twitch", "comments"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "date", "InputFormat", "dd/MM/yy");
opts = setvaropts(opts, "starttime", "InputFormat", "HH:mm");
opts = setvaropts(opts, "stoptime", "InputFormat", "HH:mm");
% opts = setvaropts(opts, "session", "TrimNonNumeric", true);
% opts = setvaropts(opts, "session", "ThousandsSeparator", ",");

% Import the data
VPmonkeysessionList = readtable(filename, opts);

end