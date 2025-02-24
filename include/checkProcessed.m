function missing = checkProcessed(clusters,markerDict,markerStructRef,cf,lp)
markerDictRef = markerStruct2dict(markerStructRef);
missing = false;
try
for cc = 1:length(clusters)
    cluster = clusters{cc};
    for mm = 1:length(cluster)
        marker = cluster(mm);
        maxErrThresh = 160;
        errThresh = 130;
        if contains(marker{1},caseInsensitivePattern('THI')) % Special case for thigh markers (prone to higher errors)
            maxErrThresh = 190;
            errThresh = 140;
        end
        
        ind = strcmp(marker,cluster);
        otherMarkers = cluster(~ind);
        
        arr = lookup(markerDict,marker);
        arr = arr{1}(:,2:4);
        arrRef = lookup(markerDictRef,marker);
        arrRef = arrRef{1};
        arrRef = [arrRef(1,2:4)]';
        donorArr = zeros(length(arr),3,length(otherMarkers));
        targArr = zeros(3,length(otherMarkers));
        for om = 1:length(otherMarkers)
            otherMarker = otherMarkers(om);
            temp = lookup(markerDict,otherMarker);
            donorArr(:,:,om) = temp{1}(:,2:4);
            temp = lookup(markerDictRef,otherMarker);
            targArr(:,om) = [temp{1}(1,2:4)]';
        end
        mErr = NaN(length(arr),1);
        
        for i = 1:length(arr)
        donorTraj = [];
        targetTraj = [];
            if ~isnan(arr(i,2))
            for j = 1:length(otherMarkers)
                if ~isnan(donorArr(i,1,j))
                    donorTraj(1:3,end+1) = squeeze(donorArr(i,1:3,j));
                    targetTraj(1:3,end+1) = targArr(1:3,j);

                end
                
            end
            if size(donorTraj,2) >= 3

            p = absor(targetTraj,donorTraj);
                            
            M = p.R * arrRef + p.t;
            mErr(i) = sqrt((M(1) - arr(i,1)).^2 + (M(2) - arr(i,2)).^2 + (M(3) - arr(i,3)).^2);
            

            end

            end
            
        end
        fs = cf; %low pass filter
        fn=(fs/2); %nyquist freq. (1/2 the sample rate)
        fc=lp; %cutoff freq.
        [b,a]=butter(4,fc/fn); %4th order butterworth filter
        filtMErr = filtfilt(b,a,mErr);
        if max(mErr) > maxErrThresh
            disp([marker{1} ' has too much peak marker error: ' num2str(max(mErr))])
            missing = true;
        elseif rmse(zeros(length(filtMErr),1),filtMErr) > errThresh
            missing = true;
            disp([marker{1} ' has too much total marker error: ' num2str(rmse(zeros(length(filtMErr),1),filtMErr))])
        end
        mErr = sort(mErr);
        if max(diff(mErr)) > 50
            missing = true;
            warning([marker{1} ' has large discontinuities'])
            
        end
        
    end
end
catch
    warning('Problem with Trial')
    missing = true;

end

end