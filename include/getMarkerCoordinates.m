function coordinate = getMarkerCoordinates(markerDict,markerName,frameNums)
coordinate = markerDict({markerName});

coordinate = coordinate{:};
coordinate = [coordinate(frameNums,2)'; coordinate(frameNums,3)'; coordinate(frameNums,4)'];

end