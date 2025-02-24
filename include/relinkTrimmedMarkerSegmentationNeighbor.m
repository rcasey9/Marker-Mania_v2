
function [markerDict,markersetToFill,relinkedFlag] = relinkTrimmedMarkerSegmentationNeighbor(markerDict,markerDictRef,clusters,refmarkerset,cluster_jump_threshold,debugFlag,verbose)
if verbose
disp('%%%%%Relink Trimmed Marker Segments Using Neighbor Frame%%%%%')
end
markerStructnames = keys(markerDict);
markerStructname = markerStructnames{1};
Frames = markerDict({markerStructname});
Frames = Frames{:};
startFrameOffset = Frames(1,1);
totalFrames = size(Frames,1);

markerSet = markerStructnames;
CmarkerSet = markerSet(contains(markerSet,'C_'));
markerSegDict = segmentMarkers(markerDict,verbose);

relinkedFlag = 0;
foundFlag = 0;
gaps = {1,2,3,4,5,6,7,8,9,10};
segRangeMax = max(cell2mat(gaps));
reverseStr = '';
markersetToFill = {};


for mm = 1:length(refmarkerset) % loop through marker set
    currentMarker = refmarkerset{mm};
    if ~isKey(markerDict,{currentMarker})
        continue
    end
    if verbose
    disp(['Searching for: ', currentMarker])
    end
    markerSeg = markerSegDict({currentMarker});
    markerSeg = markerSeg{:};
    reverseStr = '';
    for ii = 1:length(markerSeg)
        segLoc = markerSeg(ii);
        foundFlag = 1;
        msg = ['Segments To Go: ' num2str(length(markerSeg) - ii) newline];
        if verbose
        msg = sprintf('Segments To Go: %d\n', length(markerSeg) - ii);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
        while foundFlag
            foundFlag = 0;
            data = getMarkerCoordinates(markerDict,currentMarker,1:totalFrames)';
            markerSet = keys(markerDict);
            temp_names=markerSet(contains(markerSet,'C_'));
            for jj = 1:length(temp_names)
                currentFound = 0;
                cmp_marker = temp_names{jj};
                cmpX = getMarkerCoordinates(markerDict,cmp_marker,1:totalFrames)';
                currentX = getMarkerCoordinates(markerDict,currentMarker,1:totalFrames)';
                
                hasIntersection = any(~isnan(cmpX(:,1)) & ~isnan(currentX(:,1)));
                if hasIntersection || ~isKey(markerDict,{cmp_marker}) || strcmp(currentMarker,cmp_marker)
                    continue
                end




                for gg = 1:length(gaps)
                    gap = gaps{gg};
                    if mod(ii,2)
                        checkLoc = segLoc-gap;
                        segRange = segLoc:segLoc+segRangeMax;
                    else
                        checkLoc = segLoc+gap;
                        segRange = segLoc-segRangeMax:segLoc;
                    end
                    if checkLoc < 1 || checkLoc > totalFrames || ~isnan(data(checkLoc,1))
                        continue
                    end
                    if any(segRange < 1)
                        segRange = 1:segLoc;
                    elseif any(segRange > totalFrames)
                        segRange = segLoc:totalFrames;
                    end
                    CmarkerSeg = markerSegDict({cmp_marker});
                    CmarkerSeg = CmarkerSeg{:};
                    cmp_data = getMarkerCoordinates(markerDict,cmp_marker,checkLoc)';           
                    if ~any(isnan(cmp_data)) && any(ismember(checkLoc,CmarkerSeg))
                        currentMarkerCoordinate = data(segLoc,:);
                        if any(segRange<1) || any(segRange>totalFrames) || any(isnan(data(segRange,1))) || length(segRange)<2
                            referenceVel = 1;
                        else
                            referenceVel = mean(sum(diff(data(segRange,:)).^2,2).^0.5);
                        end
                        referenceDiff = (sum((currentMarkerCoordinate-cmp_data).^2)).^0.5;
                        referenceDiffAxis = currentMarkerCoordinate-cmp_data;
                        referenceVelAxis = mean((diff(fillmissing(data(segRange,:),'linear',1))),'omitnan');
                        if any(isnan(referenceVelAxis)) || length(segRange) == 1
%                             disp('          No Reference Vel, Set to default value')
                            referenceVelAxis = [1,1,1];
                        end
%                         if referenceDiff <= (referenceVel*gap+5)

                        if mod(ii,2)
                            thresholdCheck = -referenceDiffAxis + referenceVelAxis*gap;
                        else
                            thresholdCheck = referenceDiffAxis + referenceVelAxis*gap;
                        end

                        if all(abs(thresholdCheck) < 20)
%                             check = all(abs(referenceDiffAxis-referenceVelAxis*gap) < 10);
%                             check = abs(referenceDiffAxis-referenceVelAxis*gap) < 10;
                            cmpX = getMarkerCoordinates(markerDict,cmp_marker,1:totalFrames)';
                            cmp_data_range = find(~isnan(cmpX(:,1)));
                            starting = cmp_data_range(1);
                            ending = cmp_data_range(end);
                            cmp_data_range = cmp_data_range(1):cmp_data_range(end);
                            
                            msg = ['     NeighborFrame Found(',num2str(thresholdCheck,3),') : ', currentMarker,' with ',cmp_marker, ' at frames: ',num2str(cmp_data_range(1)+startFrameOffset-1),' - ',num2str(cmp_data_range(end)+startFrameOffset-1)];
                            if verbose
                            disp(msg)
                            end
                            markerDict = replaceCurrentMarkersWithCMarkers(markerDict,currentMarker,cmp_marker,cmp_data_range);

                            jumped = 0;
                            if jumped
                                markerStruct = markerStructBefore;
                                if verbose
                                disp('     Relink Incorrectly, Revert Back')
                                end
                            else
                                markerDict({cmp_marker}) = [];
                                relinkedFlag = 1;
                                markersetToFill = [markersetToFill(:)',currentMarker];
                                foundFlag = 1;
                                currentFound = 1;
                                
                                % For Debug
                                if debugFlag
                                    currentMarkerCord = markerDict({currentMarker});
                                    currentMarkerCord = currentMarkerCord{:};

                                    figure(1)
                                    clf
                                    sgtitle(currentMarker)
                                    subplot(3,1,1)
                                    plot(currentMarkerCord(:,1),currentMarkerCord(:,2),'b.-')
                                    hold on
                                    subplot(3,1,2)
                                    plot(currentMarkerCord(:,1),currentMarkerCord(:,3),'b.-')
                                    hold on
                                    subplot(3,1,3)
                                    plot(currentMarkerCord(:,1),currentMarkerCord(:,4),'b.-')
                                    hold on
                                    subplot(3,1,1)
                                    plot(currentMarkerCord(cmp_data_range,1),currentMarkerCord(cmp_data_range,2),'r.-')
                                    subplot(3,1,2)
                                    plot(currentMarkerCord(cmp_data_range,1),currentMarkerCord(cmp_data_range,3),'r.-')
                                    subplot(3,1,3)
                                    plot(currentMarkerCord(cmp_data_range,1),currentMarkerCord(cmp_data_range,4),'r.-')
                                end
                                break;
                            end
                            
                        end
                    end
                end

                if currentFound
                    break
                end
            end
            if mod(ii,2) && foundFlag
                segLoc = cmp_data_range(1);
            elseif ~mod(ii,2) && foundFlag
                segLoc = cmp_data_range(end);
            end
            if segLoc < 1 || segLoc > totalFrames
                foundFlag = 0;
            end
        end
        if ~contains(msg,'Segment')
            msg = newline;
            if verbose
            fprintf(msg);
            end
        end
        if verbose
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end

% if debugFlag
%     saveas(gcf,['debugFig\',currentMarker,'_Neighbor.png'])
% end
end
markersetToFill = unique(markersetToFill);
end