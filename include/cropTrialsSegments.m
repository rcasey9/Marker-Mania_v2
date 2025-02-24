function  markerStructLists = cropTrialsSegments(markerStruct,cropLength)
% To crop the trials to small ones and later combine them together

markerStructLists = struct();

tempMarker = fieldnames(markerStruct);
tempMarker = tempMarker{1};
markerHeader = markerStruct.(tempMarker).Header;
totalLength = length(markerHeader);

numSegments = ceil(totalLength/cropLength);

for i = 1:numSegments
    i_start = (i-1)*cropLength+1;
    if i ~= numSegments
        i_end = i*cropLength;
    else
        i_end = totalLength;
    end

    SegName = ['Seg_',num2str(i)];
    markerStructLists.(SegName) = markerStruct;
    markerNames = fieldnames(markerStructLists.(SegName));

    for j = 1:length(markerNames)
        markerName = markerNames{j};
        markerTable = markerStructLists.(SegName).(markerName);
        markerTable = markerTable(i_start:i_end,:);
        markerStructLists.(SegName).(markerName) = markerTable;
    end
end




end