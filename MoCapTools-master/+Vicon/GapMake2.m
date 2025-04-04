function [errorTable, markerData, newGapBoolean] = GapMake2(markerData, ikXml, varargin)
% GapMake: Runs inverse kinematics on a struct of marker data to determine
% where there are high errors, and deletes those frames with high error
% [errorTable, markerData, newGapBoolean] = GapMake(markerData, ikXml, varargin)
%   In:
%       markerData - a struct of (gap-filled) marker data, such as one
%         generated by Vicon.ExtractMarkers
%       ikXml - path to IK setup file (.xml)
%   Optional Inputs:
%       VerboseLevel - 0 (minimal, default), 1 (normal), 2 (debug mode)
%       ErrorThresholdLow - threshold of IK error above which bad marker
%         data will be deleted
%       ErrorThresholdHigh - threshold of IK error above which marker data
%         will be considered bad
%       
%   Out:
%       errorTable - a table containing IK error for every marker at every
%         frame
%       markerData - updated struct of marker data with gaps created
%       newGapBoolean - true if there were gaps created, false otherwise
%
%   See also: Vicon.GapFill, Vicon.IterativeGapFilling.

% Any region where the error exceeds lowThresh will be flagged as having
% high error, but only those regions where the error ever exceeds
% highThresh will actually be deleted. So if the error exceeds lowThresh
% but not highThresh, nothing will happen in that region. If the error
% exceeds both, then the entire region exceeding lowThresh will be deleted.
% Maybe these should change over iterations?

% Have different thresholds for joints vs. skin markers?

% Have another input for a doWeCare vector? 

p = inputParser;
p.addParameter('VerboseLevel',0);
p.addParameter('ErrorThresholdLow',0.04);
p.addParameter('ErrorThresholdHigh',0.06);
 
p.parse(varargin{:});
 
verboseLevel = p.Results.VerboseLevel;
lowThresh = p.Results.ErrorThresholdLow;
highThresh = p.Results.ErrorThresholdHigh;


%% Identify and delete bad frames

% Remove frames that have absolutely no data
trcTable = Osim.interpret(markerData, 'TRC');
noDataFrames = all(isnan(trcTable{:, 2:end}), 2);
if any(noDataFrames)
    warning('Removing frames with absolutely no data');
    trcTable=trcTable(~noDataFrames,:);
end

% run inverse kinematics to get errors

ikwarnings = evalc('[~, errorTable, ~] = Osim.IK(trcTable, ikXml);');
if verboseLevel == 2
    fprintf('%s\n', ikwarnings);
    % visualize the error data
    figure
    Osim.viewIKError(errorTable);
    figure;
    % visualize with a surface plot
    Osim.viewIKError(errorTable, 1);
end

fprintf('   Max Error of Gap Region: %f\n',max(max(errorTable{:,2:end})));

% get what frames need to be deleted for each marker based on which parts
% of the trial are important (determined from fpFile) and which parts have
% high error (determined by the other three inputs)

% change this to include a doWeCare vector? 
[markerNames, frames] = Osim.getFramesToDelete2(errorTable); %, fpFile);

if verboseLevel >= 1
    fprintf('   Num Gaps Made: %d\n',length(frames));
end
% markerNames is a cell array like {'L_ASIS', 'R_ASIS'}
% frames is a cell array of vectors like {[1:500], [200:300, 400:450]}
if isempty(frames)
    fprintf('   No regions with exeedingly high error!\n');
    newGapBoolean = false;
    return;
end
newGapBoolean = true;
% for each marker, replace its bad frames with NaN
markerData = deleteFrames(markerData, markerNames, frames);
% markerData can be gap filled again
end
