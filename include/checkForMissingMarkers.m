function [missing] = checkForMissingMarkers(markerStruct, markerSet,verbose)
missing = false;
for mm = 1:length(markerSet)
    currentMarker = markerSet{mm};
    if contains(currentMarker(1:2),'C_')==0
        if any([any(isnan(markerStruct.(currentMarker).x)), any(isnan(markerStruct.(currentMarker).y)),any(isnan(markerStruct.(currentMarker).z))])
        missing = true;
        if verbose
        disp(['    MISSING:',currentMarker])
        end
        end
    end
end

end