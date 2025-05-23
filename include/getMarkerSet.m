function [clusters, clusters_jump_threshold, markerStructRef] = getMarkerSet(folderPath,viconPath,custom_jump_thresholds)
fprintf('\n \n  %%%%%% GETTING MARKERSET %%%%%% \n \n');
 if ~exist(folderPath, 'dir')
       error(['Folder Does Not Exist: ' newline folderPath newline])
end
[~,status_result] = system('tasklist /FI "imagename eq nexus.exe" /fo table /nh');
%Check if Vicon is Running
if ~contains(status_result, 'Nexus.exe')
    %Open Vicon if it isn't running
    system([viconPath ' &'])
    pause(20)
end
vicon = ViconNexus();

files = dir(folderPath);
L = length(files);
index = false(1, L);
for k = 1:L
    M = length(files(k).name);
    if M > 4 && contains(files(k).name, '.c3d')
        index(k) = true;
    end
end
files = files(index);
L = length(files);
index = false(1, L);
for k = 1:L
    M = length(files(k).name);
    if (contains(files(k).name, 'static') || contains(files(k).name, 'Static') || contains(files(k).name, 'STATIC'))
        index(k) = true;
    end
end
staticFiles = files(index);
if isempty(staticFiles)
    error(['No Static Trial Found in Folder: ' newline folderPath newline])
end

File = staticFiles(1).name;
filename = [folderPath '\' File(1:length(File)-4)];
disp(['Using: ' File(1:length(File)-4)])
doing_vicon_operations = true;
while doing_vicon_operations
try   
vicon.OpenTrial(filename, 60);
vicon.RunPipeline('ExportC3D', '', 200);
vicon.SaveTrial(60);
[ names, ~, active ] = vicon.GetSubjectInfo();
subject = names{active};
disp(['Found Subject: ' subject])
markerStructRef = Vicon.ExtractMarkers([filename,'.c3d']);

segments = vicon.GetSegmentNames(subject);
clusters = {};
for qq = 1:length(segments)
    segment = segments{qq};
    [parent, ~, markers] = vicon.GetSegmentDetails(subject, segment);
    jj = 1;
    %%Check if marker has been detached or isn't labeled
    while jj <= length(markers)
        marker = markers{jj};
        try
            check = markerStructRef.(marker).x(1);
        catch 
            markers(jj) = [];
            continue
        end
        if isnan(check)
             markers(jj) = [];
             continue
        end
        jj = jj+1;
    end
    if length(markers) < 4
        diff = 4-length(markers);
        [~, ~, parentMarkers] = vicon.GetSegmentDetails(subject, parent);
        clusterCenter = [0;0;0];
        markerString = [];
        for ll = 1:length(markers)
            marker = markers{ll};
            markerString = [markerString marker ','];
            clusterCenter = clusterCenter + [markerStructRef.(marker).x(1);markerStructRef.(marker).y(1);markerStructRef.(marker).z(1)]./length(markers);
        end
        dists = [];
        for pm = 1:length(parentMarkers)
            parentMarker = parentMarkers{pm};
            dists(end+1) = sqrt((clusterCenter(1)-markerStructRef.(parentMarker).x(1)).^2 + (clusterCenter(2)-markerStructRef.(parentMarker).y(1)).^2 + (clusterCenter(3)-markerStructRef.(parentMarker).z(1)).^2);
        end
        [~,inds] = sort(dists,'ascend');
        for jj = 1:diff
            ind = inds(jj);
          markers{end + 1} = parentMarkers{ind};
          disp(['Adding ' parentMarkers{ind} ' to cluster: ' markerString])
        end  
    end
    clusters{qq} = markers;
    clusters_jump_threshold{qq} = {20};
    for cl = 1:length(custom_jump_thresholds)
        if strcmpi(segment,custom_jump_thresholds{cl}{1})
            clusters_jump_threshold{qq} = custom_jump_thresholds{cl}(2);
        end
    end
end
vicon.CloseTrial(60);
catch

    createEndnoteFilter(folderPath,filename)
    warning('Problem communicating with Vicon... Attempting to reconnect')
    checkReopenVicon(viconPath);
    Vicon_Openned = false;
    while ~Vicon_Openned
    try
    vicon = ViconNexus();
    catch 
        pause(2)
        continue
    end
    Vicon_Openned = true;
    end
    continue

end
doing_vicon_operations = false;
end
end
