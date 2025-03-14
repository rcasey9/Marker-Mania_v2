clc; close all; clear;

directory = pwd;
addpath([directory '\include'])
addpath([directory '\MocapUtils'])
addpath([directory '\MoCapTools-master'])
addpath([directory '\MoCapTools-master\lib'])
addpath([directory '\MoCapTools-master\extlib'])
addpath([directory '\MoCapTools-master\extlib\btk'])

ViconFolder = [];
files = dir('C:\Program Files\Vicon');
for ii = length(files):-1:1
    folder = files(ii).name;
    if contains(folder,'Nexus2')
        ViconFolder = folder;
        break
    end
end
if ~isempty(ViconFolder)
    addpath(['C:\Program Files\Vicon\' ViconFolder '\SDK\MATLAB'])
end
files = dir('C:\Program Files (x86)\Vicon');
for ii = length(files):-1:1
    folder = files(ii).name;
    if contains(folder,'Nexus2')
        ViconFolder = folder;
        break
    end
end
if ~isempty(ViconFolder)
    addpath(['C:\Program Files (x86)\Vicon\' ViconFolder '\SDK\MATLAB'])
else
    warning('Vicon MATLAB SDK not found. Add it to your MATLAB path manually')
end



files = dir('Pipelines');
for ii = length(files):-1:1
    file = files(ii).name;
    if contains(file,'.Pipeline')
        copyfile([directory '\Pipelines\' file], ['C:\Users\Public\Documents\Vicon\Nexus2.x\Configurations\Pipelines\' file])
    end
end