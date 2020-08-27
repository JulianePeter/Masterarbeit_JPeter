%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEMI-AUTOMATIC NEEDLE SEGMENTATION VERSION 2 
% for HDR-Brachytherapy of the Prostate in TRUS data at UKSH Kiel.
%
% This algorithm is to be used on TRUS data that have not been segmented
% manually. It automatically segments the implant needles in a TRUS dataset 
% after manual ROI selection. The resulting coordinates of all segmented 
% needles are stored in the array finalNeedlePointsCoordinates.
%
% Interaction required by the user:
% 1) The path to the TRUS dataset to load needs to be specified in line 52.
% 2) An elliptical region of interest must be drawn on the centralmost 
% axial slice presented to the user covering the prostate and the needles. 
% The addition of a small rim (a few millimeters) around the visible
% needles is usually sufficient. The short semiaxis should divide the 
% prostate into two virtually symmetric parts. The ellipse dimensions can 
% be adapted and the ellipse can be shifted and rotated after creation. The 
% ROI selection is confirmed by double clicking on it.
%
% Additional .m-files required:
% loadTRUS
% findLocalMaxima2D
% mergeIdenticalMaxima
% imageAndNeedles2rgb
% label2rgb3d
% superimposeImageAndSegmentation
% getCurrentColumnPosition
% getIntensities
% coordinates2voxels
% voxels2coordinates
% imshow3Dfull (modified by Juliane Peter)
% imshow3Dfullcolour (imshow3Dfull modified by Juliane Peter)
%
% Further information on the segmentation algorithm can be found in the
% documents 'Hinweise_Segmentierungsalgorithmus_Implantationsnadeln.pdf' 
% and 'Masterarbeit_JPeter.pdf'.
%
% This code was written in Matlab R2020a.
%
% Author: Juliane Peter
% E-mail: juliane-peter@web.de
% July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
% Start timer.
%tic
%% Loading of Patient Dataset.
% Load a patient dataset in DICOM format, getting the TRUS image data,
% DICOM information.
% The user is required to insert the path to the DICOM data here:
[image,info] = loadTRUS('');

% Squeeze image dimensions.
image = squeeze(image);

% Store original image before further modification.
originalImage = image;

% Define colourmap.
cmap = colormap('jet');

% Add empty slice if the number of image slices is even to ensure to work on
% an odd number of slices comprising the whole image stack.
sliceAppended = false;
if mod(size(image,3),2) == 0
    image = zeros(size(originalImage,1),size(originalImage,2),size(originalImage,3)+1);
    image(:,:,1:1:size(originalImage,3)) = originalImage;
    image = uint8(image);
    sliceAppended = true;
end

%% ROI Selection and Characterisation of Cropped Image and ROI Ellipse.
% Determine central slice to select ROI on.
centralSliceNumber = int16(ceil(size(image,3)/2));
centralSlice = image(:,:,centralSliceNumber);

% Show central slice.
imshow(centralSlice,[0 255])

% Let the user draw an elliptical ROI covering the prostate area and all
% needles. The ellipse should neither be too small nor too wide and the 
% short semiaxis should divide the prostate into two virtually symmetric 
% parts. ROI selection is confirmed by double clicking on it.
roi = drawellipse;

% Stop further code execution until the ROI is selected by the user.
wait(roi)

% Create logical mask from the elliptical ROI.
roiMask = uint8(createMask(roi,image));

% Use the mask to mask the TRUS image.
prostateRegion = zeros(size(image,1,2,3));
for i = 1:size(image,3)
    prostateRegion(:,:,i) = roiMask.*image(:,:,i);
end

% Close the ROI selection image.
close

% Stop first timer.
%toc
% Start timer number two.
%tic

% Get the position and the dimensions of the smallest rectangle containing 
% the ROI ellipse. Use this rectangle to select a cuboid image space applying 
% it to all slices of the image stack.
boundingBox = regionprops(roiMask,'BoundingBox');

% Define cropping dimensions.
xMin = boundingBox.BoundingBox(1);
yMin = boundingBox.BoundingBox(2);
zMin = 1;
width = boundingBox.BoundingBox(3);
height = boundingBox.BoundingBox(4);
depth = size(image,3)-1;

% Crop the smallest cuboid shape containing the ellipse.
croppedImage = imcrop3(prostateRegion,[xMin yMin zMin width height depth]);
croppedMask = imcrop(roiMask,[xMin yMin width height]);

% Use edge detection to identify and store the ROI ellipse's vertices.
croppedMaskRim = edge(croppedMask);
[ellipseVertices(:,1),ellipseVertices(:,2)] = find(croppedMaskRim);

% Determine missing ellipse vertices that come to lie on the rectangular 
% border.
firstCounter = 1;
secondCounter = 1;
for i = 1:size(ellipseVertices,1)
    if ellipseVertices(i,1) == 1
        limitsFirstDim{firstCounter} = ellipseVertices(i,2);
        firstCounter = firstCounter+1;
    elseif ellipseVertices(i,2) == 1
        limitsSecondDim{secondCounter} = ellipseVertices(i,1);
        secondCounter = secondCounter+1;
    end
end
fillingFirstDim = linspace(limitsFirstDim{1}+1,limitsFirstDim{2}-1,limitsFirstDim{2}-...
    limitsFirstDim{1}-1);
fillingSecondDim = linspace(limitsSecondDim{1}+1,limitsSecondDim{2}-1,limitsSecondDim{2}-...
    limitsSecondDim{1}-1);

verticesToAdd = ones(length(fillingFirstDim)+length(fillingSecondDim),2);
verticesToAdd(1:1:length(fillingFirstDim),2) = fillingFirstDim;
verticesToAdd(length(fillingFirstDim)+1:1:length(fillingFirstDim)+length(fillingSecondDim),1) = ...
    fillingSecondDim;

% Add these additonal vertices to the ellipse vertices obtained by edge
% detection.
ellipseVertices = cat(1,ellipseVertices,verticesToAdd);

% Store the cropped image in the image variable to perform the segmentation 
% on.
image = croppedImage;

% Determine the image centre in all three dimensions.
centreX = int16(floor(size(image,1)/2));
centreY = int16(floor(size(image,2)/2));
centreZ = int16(floor(size(image,3)/2));

% Show the cropped image.
% figure
% imshow3Dfull(image)

% Compute the sagittal image view for the cross correlation with the needle model.
sagittalImage = flip(permute(image,[3 1 2]));

%% Creation of the Needle Model.
% Define the needle model's dimensions. 
needleLength = size(image,3);
needleWidth = 5;
needleDepth = 5;
halfLength = floor(needleLength/2);
halfWidth = floor(needleWidth/2);
halfDepth = floor(needleDepth/2);

% Define a longitudinal elongation to make sure that the needle still spans  
% the needle model space irrespective of its orientation after rotation.
longitudinalElongation = 6;
halfLongitudinalElongation = longitudinalElongation/2;

% Determine the box dimensions of the needle model space.
needleBoxDimension = int16(ceil(sqrt(needleDepth*needleDepth+needleLength*needleLength)));
% Make sure that the box dimensions are odd to have a well defined centre.
if mod(needleBoxDimension,2) == 0
    needleBoxDimension = needleBoxDimension+1;
end

% Determine the centre of the needle model array.
needleBoxCentre = ceil(needleBoxDimension/2);

% Create the needle box.
needleStraight = zeros(needleBoxDimension+longitudinalElongation,needleBoxDimension,needleBoxDimension);
% Add the straight needle core as a line of ones.
needleStraight(:,needleBoxCentre,needleBoxCentre) = 1;

% Define a radially weighted axial needle plane to combine with the rotated
% straight needle core.
needleAxialPlane = ones(5);
needleAxialPlane(2:4,2:4) = 2;
needleAxialPlane(3,3) = 3;

% Define shortened axial plane versions for needle points near the rim.
lowerRimAxialPlane = needleAxialPlane(3:1:5,:);
lowerRimPlusOneAxialPlane = needleAxialPlane(2:1:5,:);
upperRimAxialPlane = needleAxialPlane(1:1:3,:);
upperRimMinusOneAxialPlane = needleAxialPlane(1:1:4,:);

%% Cross Correlation of Sagittal Image Patches with the Rotated Needle Model.
% Select rotation angles in dorso-ventral direction i. e. rotation in the 
% sagittal plane.
dorsoVentralAngles = linspace(-40,10,6);

localMax = cell({0});
maxCounter = 1;
localMaxCounter = 1;
% For all rotation angles.
for z=1:length(dorsoVentralAngles)
    % Rotate the needle model in dorso-ventral direction.
    needleRotated = imrotate3(needleStraight,dorsoVentralAngles(z),[0 0 1],'nearest','crop');
    
    % Reduce the elongated rotated needle to overlay dimensions.
    needleRotatedTemp = needleRotated(1+halfLongitudinalElongation:1:size(needleStraight,1)-...
        halfLongitudinalElongation,:,:);
    needleRotated = zeros(size(needleRotatedTemp,1,2,3));
    for i = 1:size(needleRotatedTemp,1)
        slice = squeeze(needleRotatedTemp(i,:,:));
        % Identify the position of the needle core on every slice.
        [row,column] = find(slice,1);
        % Combine the needle core with the radially weighted axial needle 
        % planes to obtain a rotated weighted needle model for cross correlation.
        if ~isempty(row)
            if row == 1
                needleRotated(i,row:1:row+2,column-2:1:column+2) = lowerRimAxialPlane;
            elseif row == 2
                needleRotated(i,row-1:1:row+2,column-2:1:column+2) = lowerRimPlusOneAxialPlane;
            elseif row == size(slice,2)
                needleRotated(i,row-2:1:row,column-2:1:column+2) = upperRimAxialPlane;
            elseif row == size(slice,2)-1
                needleRotated(i,row-2:1:row+1,column-2:1:column+2) = upperRimMinusOneAxialPlane;
            else
                needleRotated(i,row-2:1:row+2,column-2:1:column+2) = needleAxialPlane;
            end
        end
    end   
    
    % Create an array to store the alignment values of the needle model 
    % and the image patch obtained in the successive cross correlation.
    alignment = zeros(3,(size(sagittalImage,2)-needleBoxDimension+1)*...
        (size(sagittalImage,3)-needleBoxDimension+1));
    counter = 1;
    patch = zeros(needleBoxDimension,needleBoxDimension,needleBoxDimension);
    for i=1:size(sagittalImage,2)-needleBoxDimension+1 
        for j=1:size(sagittalImage,3)-needleBoxDimension+1
            % Get an image patch which is translated in every loop iteration.
            patch(needleBoxCentre-halfLength:1:needleBoxCentre+halfLength,:,:) = ...
                sagittalImage(:,i:1:i+needleBoxDimension-1,j:1:j+needleBoxDimension-1);
            % Store the voxel indices of the needle core on the central slice.
            alignment(1,counter) = i+needleBoxCentre;
            alignment(2,counter) = j+needleBoxCentre;
            % Cross correlation: Multiply the image patch and the needle model and 
            % sum the voxelwise products as alignment value. 
            alignment(3,counter) = sum(sum(sum(patch.*needleRotated)));
            counter=counter+1;
        end
    end

    % Find the local maxima of the obtained alignment values for one 
    % rotational configuration. The alignment values are sorted into  
    % a 2D array at positions corrsponding to the needle centre overlay
    % position.
    identifiedLocalMaxima = findLocalMaxima2D(alignment,[],i,j,needleBoxCentre,false);

    % Assess identified local maxima in close neighbourhood further and find 
    % the local maximum of these local maxima within a Euclidean distance of 19.
    [identifiedLocalMaxima,~] = mergeIdenticalMaxima(identifiedLocalMaxima,19);

    % Compare newly identified local extrema for the current orientation 
    % to already identified local extrema in former iterations and assess 
    % for doublings.
    for n = 1:size(identifiedLocalMaxima,2)
        % If this is the very first local maximum detected, store it. 
        if localMax{1} == 0
            localMaxI{maxCounter} = identifiedLocalMaxima(1,n);
            localMaxJ{maxCounter} = identifiedLocalMaxima(2,n);
            localMax{maxCounter} = identifiedLocalMaxima(3,n);
            detectedNeedle{maxCounter} = needleRotated;
            needleDorsoVentralAngle{maxCounter} = dorsoVentralAngles(z);
            maxCounter = maxCounter+1;
        else
            % For all local maxima already saved test whether the newly 
            % identified local maximum is close to an existing local maximum. 
            % If so, replace the old local maximum in case the new local
            % maximum is higher.
            for m = 1:size(localMax,2)
                if abs(identifiedLocalMaxima(1,n)-localMaxI{m})<20 &&...
                        abs(identifiedLocalMaxima(2,n)-localMaxJ{m})<=20
                    if identifiedLocalMaxima(3,n) > localMax{m}
                        localMaxI{m} = identifiedLocalMaxima(1,n);
                        localMaxJ{m} = identifiedLocalMaxima(2,n);
                        localMax{m} = identifiedLocalMaxima(3,n);
                        detectedNeedle{m} = needleRotated;
                        needleDorsoVentralAngle{maxCounter} = dorsoVentralAngles(z);
                    end
                else 
                    localMaxCounter = localMaxCounter+1;
                end
            end
            % If the newly identified local maximum is not in close
            % proximity to another identified local maximum, store it as a
            % new local maximum.
            if localMaxCounter > size(localMax,2)
                localMaxCounter = 1;
                localMaxI{maxCounter} = identifiedLocalMaxima(1,n);
                localMaxJ{maxCounter} = identifiedLocalMaxima(2,n);
                localMax{maxCounter} = identifiedLocalMaxima(3,n);
                detectedNeedle{maxCounter} = needleRotated;
                needleDorsoVentralAngle{maxCounter} = dorsoVentralAngles(z);
                maxCounter = maxCounter+1;
            end
        end
    end
end

% After cross correlation in all orientations, store the needle core indices on 
% the central slice, the alignment values and a reference to the respective 
% needle model rotational configuration and rotation angle in one array.
finalAlignment = zeros(6,size(localMax,2));
for d = 1:size(localMax,2)
    finalAlignment(1,d) = localMaxI{d};
    finalAlignment(2,d) = localMaxJ{d};
    finalAlignment(3,d) = localMax{d};
    finalAlignment(4,d) = d;
    finalAlignment(5,d) = needleDorsoVentralAngle{d};
end

%% Definition of Potential Needle Detections though Assessment of Cross Correlation Maxima.
% Identify the local maxima of all local maxima obtained for all needle 
% orientations in 2D.
newFinalAlignment = findLocalMaxima2D(finalAlignment,needleDorsoVentralAngle,i,j,needleBoxCentre,true);

% Find the local maximum of the identified local maxima in close neighbourhood.
[newFinalAlignment,~] = mergeIdenticalMaxima(newFinalAlignment,19);

% Save finally identified local maxima as potential needle positions.
for i = 1:size(newFinalAlignment,2)
    newLocalMaxI{i} = newFinalAlignment(1,i);
    newLocalMaxJ{i} = newFinalAlignment(2,i);
    newLocalMax{i} = newFinalAlignment(3,i);
    newDetectedNeedle{i} = detectedNeedle{newFinalAlignment(4,i)};
end

% Use needle cores corresponding to the alignment local maxima to define 
% potential needle detections for segmentation.
sagittalDetections = zeros(size(sagittalImage,1,2,3));
for s = 1:size(newLocalMax,2)
    % Define needles as intense centre of the needle model and label each 
    % needle.
    needleCore = s*double(newDetectedNeedle{s}==3);
    % Sort into saggital segmentation array.
    sagittalDetections(:,newLocalMaxI{s}-needleBoxCentre:1:newLocalMaxI{s}-needleBoxCentre+needleBoxDimension-1,...
        newLocalMaxJ{s}-needleBoxCentre:1:newLocalMaxJ{s}-needleBoxCentre+needleBoxDimension-1) = ...
        needleCore(needleBoxCentre-halfLength:1:needleBoxCentre+halfLength,:,:);
end

% Compute axial view of this initial segmentation.
axialDetections = flip(permute(sagittalDetections,[2 3 1]),3);

% Plot potentially identified needles.

% Thicken the segmentation contour for better visibility. Colour 8-neighbourhood.
se = strel('square',3);
for i = 1:size(image,3)
    axialDetectionsThick(:,:,i) = imdilate(axialDetections(:,:,i),se);
end

% Convert image and potential detections into rgb format.
[rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(image,axialDetectionsThick,cmap);
rgbImage = rgbImage*255/max(max(max(max(rgbImage))));
rgbImage = uint8(rgbImage);
   
% Superimpose image and inital needle detections.
superposition = superimposeImageAndSegmentation(rgbImage,rgbNeedles,0.4);

% Plot superposition.
%figure
%imshow3Dfullcolour(superposition,[],cmap)

% Determine the number of labels.
numberOfLabels = (max(max(max(axialDetections))));

%% Segmentation of Bent and 2D-Rotated Needles at Positions Obtained by Cross Correlation.
% Initialise array to store needle segmentation in.
needleROIs = zeros(size(image,1,2,3));

% Define an area around the initial needle detection in which the segmented 
% needle points may come to lie. These parameter values were set in an
% optimisation process to make sure that needles can be segmented
% adequately.
dorsoVentralRangeCentralSlice = 41;
lateralRangeCentralSlice = 33;
dorsoVentralExtentCentralSlice = floor(dorsoVentralRangeCentralSlice/2);
lateralExtentCentralSlice = floor(lateralRangeCentralSlice/2);
dorsoVentralRange = 13;
dorsoVentralExtent = floor(dorsoVentralRange/2);
lateralRange = 13;
lateralExtent = floor(lateralRange/2);

% Positional weighting for central slice search favouring central
% positions.
positionalWeightingCentralSlice = ones(dorsoVentralRangeCentralSlice,lateralRangeCentralSlice);
positionalWeightingCentralSlice(11:1:31,6:1:28) = 1.15;

% Slightly asymmetric dorso-ventral weighting to introduce bias for a needle 
% instead of an echo while favouring central positions.
dorsoVentralPositionalWeighting = 1.05*ones(dorsoVentralRange,lateralRange);
dorsoVentralPositionalWeighting(1:1:2,:,:) = 1.0;
dorsoVentralPositionalWeighting(dorsoVentralExtent+1:1:size(dorsoVentralPositionalWeighting,1)-...
    (dorsoVentralExtent-1),:) = 1.10;

% Symmetric dorso-ventral weighting for the inner one third of slices to avoid 
% bias besides favouring central positions.
dorsoVentralPositionalWeightingEnd = 1.1*ones(dorsoVentralRange,lateralRange);
dorsoVentralPositionalWeightingEnd(1,:,:) = 1.0;
dorsoVentralPositionalWeightingEnd(dorsoVentralRange,:,:) = 1.0;
dorsoVentralPositionalWeightingEnd(dorsoVentralExtent:1:dorsoVentralExtent+2,:) = 1.15;

% Symmetric lateral positional weighting favouring central positions.
lateralPositionalWeighting = 1.1*ones(dorsoVentralRange,lateralRange);
lateralPositionalWeighting(:,1) = 1.0;
lateralPositionalWeighting(:,lateralRange) = 1.0;
lateralPositionalWeighting(:,lateralExtent:1:lateralExtent+2) = 1.15;

% Initalise array to store the points of the obtained labelled segmentation.
if sliceAppended == true
    pointsOfLabelledSegmentation = zeros(2,size(image,3)-1,numberOfLabels);
else
    pointsOfLabelledSegmentation = zeros(2,size(image,3),numberOfLabels);
end
% For all labels.
for i = 1:numberOfLabels
    % For all slices from the central slice outwards.
    for j = (0:1:centreZ)
        % If search is performed on the central slice.
        if j == 0
            % Get the position of the labelled needle detection on the
            % central slice as determined in the cross correlation process.
            currentSliceLabelPosition = axialDetections(:,:,centralSliceNumber+j) == i;
            [y,x] = find(currentSliceLabelPosition);
            
            % Get an axial image excerpt according to the defined lateral and 
            % dorso-ventral extent around the detected needle position in the 
            % central slice and its two neighbouring slices.
            intensity1 = getIntensities(image,y,dorsoVentralExtentCentralSlice,x,...
                lateralExtentCentralSlice,centralSliceNumber,j-1);
            intensity2 = getIntensities(image,y,dorsoVentralExtentCentralSlice,x,...
                lateralExtentCentralSlice,centralSliceNumber,j);
            intensity3 = getIntensities(image,y,dorsoVentralExtentCentralSlice,x,...
                lateralExtentCentralSlice,centralSliceNumber,j+1);
            % Average the image intensities over the three central slices.
            meanIntensity = 1/3*(intensity1+intensity2+intensity3);
            
            % Apply positional weighting to obtained mean intensities.
            weightedIntensities = meanIntensity.*positionalWeightingCentralSlice;
            
            % Find the position of the maximum within the weighted mean intensity array.
            [maximum,index] = max(weightedIntensities,[],'all','linear');
            [row,column] = ind2sub(size(weightedIntensities,1,2),index);
            
            % Store the position of the maximum for the successive positional 
            % comparison for voxel classification on adjacent slices.
            currentLowerRowPosition = y-dorsoVentralExtentCentralSlice-1+row;
            currentUpperRowPosition = currentLowerRowPosition;
            currentLowerColumnPosition = getCurrentColumnPosition(x,column,2*lateralExtentCentralSlice+1);
            currentUpperColumnPosition = currentLowerColumnPosition;
            
            % Classify maximum voxel as needle voxel and store in segmentation.
            needleROIs(currentLowerRowPosition,currentLowerColumnPosition,centralSliceNumber+j) = i;
            % Store point of labelled segmentation.
            pointsOfLabelledSegmentation(1,centralSliceNumber+j,i) = currentLowerRowPosition;
            pointsOfLabelledSegmentation(2,centralSliceNumber+j,i) = currentLowerColumnPosition;
        else
            % Get the image intensities for the axial image excerpt from the
            % previous slice.
            intensity = getIntensities(image,currentLowerRowPosition,dorsoVentralExtent,...
                currentLowerColumnPosition,lateralExtent,centralSliceNumber,-j);
            
            % Create a directional weighting array taking the preceding 
            % positional change in this direction into account.
            dorsoVentralDirectionalWeighting = ones(dorsoVentralRange,lateralRange);
            lateralDirectionalWeighting = ones(dorsoVentralRange,lateralRange);
            % If this is not the first step away from the central
            % slice, i. e. if a positional change in this direction has already occurred.
            if j > 1
                if lowerRowDirection == -dorsoVentralExtent
                    dorsoVentralDirectionalWeighting(1:1:2,:) = 1.1;
                elseif lowerRowDirection == dorsoVentralExtent
                    dorsoVentralDirectionalWeighting(2*dorsoVentralExtent:1:2*dorsoVentralExtent+1,:) = 1.1;
                else
                    dorsoVentralDirectionalWeighting(lowerRowDirection+...
                        dorsoVentralExtent:1:lowerRowDirection+dorsoVentralExtent+2,:) = 1.1;
                end
                if lowerColumnDirection == -lateralExtent
                    lateralDirectionalWeighting(:,1:1:2) = 1.1;
                elseif lowerColumnDirection == lateralExtent
                    lateralDirectionalWeighting(:,2*lateralExtent:1:2*lateralExtent+1) = 1.1;
                else
                    lateralDirectionalWeighting(:,lowerColumnDirection+...
                        lateralExtent:1:lowerColumnDirection+lateralExtent+2) = 1.1;
                end
            end
            
            % Weighting of the image intensity excerpt with the positional and
            % directional weighting factors.
            weightedIntensities = intensity.*dorsoVentralPositionalWeighting.*...
                lateralPositionalWeighting.*dorsoVentralDirectionalWeighting.*lateralDirectionalWeighting;
           
            % Update positions.
            oldLowerRowPosition = currentLowerRowPosition;
            oldLowerColumnPosition = currentLowerColumnPosition;
            
            % Find position of the maximum within the weighted image excerpt.
            [maximum,index] = max(weightedIntensities,[],'all','linear');
            [row,column] = ind2sub(size(weightedIntensities,1,2),index);

            % Store positions.
            currentLowerRowPosition = oldLowerRowPosition-dorsoVentralExtent-1+row;
            currentLowerColumnPosition = getCurrentColumnPosition(oldLowerColumnPosition,column,lateralRange);
            
            % Write into segmentation array.
            needleROIs(currentLowerRowPosition,currentLowerColumnPosition,centralSliceNumber-j) = i;
            
            % Store point in segmentation point storage array.
            pointsOfLabelledSegmentation(1,centralSliceNumber-j,i) = currentLowerRowPosition;
            pointsOfLabelledSegmentation(2,centralSliceNumber-j,i) = currentLowerColumnPosition;
            
            % Determine axial segmentation direction.
            lowerRowDirection = currentLowerRowPosition-oldLowerRowPosition;
            lowerColumnDirection = currentLowerColumnPosition-oldLowerColumnPosition;
            
            % Terminate the loop if the last slice is an empty appended image slice.
            if j== centreZ && sliceAppended == true
                break
            end
            
            % Get the image intensities for the axial image excerpt from the
            % subsequent slice.
            intensity = getIntensities(image,currentUpperRowPosition,dorsoVentralExtent,...
                currentUpperColumnPosition,lateralExtent,centralSliceNumber,j);
            
            % Create a directional weighting array.
            dorsoVentralDirectionalWeighting = ones(dorsoVentralRange,lateralRange);
            lateralDirectionalWeighting = ones(dorsoVentralRange,lateralRange);
            
            % If this is not the first step away from the central
            % slice.
            if j > 1
                if upperRowDirection == -dorsoVentralExtent
                    dorsoVentralDirectionalWeighting(1:1:2,:) = 1.1;
                elseif upperRowDirection == dorsoVentralExtent
                    dorsoVentralDirectionalWeighting(2*dorsoVentralExtent:1:2*dorsoVentralExtent+1,:) = 1.1;
                else
                    dorsoVentralDirectionalWeighting(upperRowDirection+...
                        dorsoVentralExtent:1:upperRowDirection+dorsoVentralExtent+2,:) = 1.1;
                end
                if upperColumnDirection == -lateralExtent
                    lateralDirectionalWeighting(:,1:1:2) = 1.1;
                elseif upperColumnDirection == lateralExtent
                    lateralDirectionalWeighting(:,2*lateralExtent:1:2*lateralExtent+1) = 1.1;
                else
                    lateralDirectionalWeighting(:,upperColumnDirection+...
                        lateralExtent:1:upperColumnDirection+lateralExtent+2) = 1.1;
                end
            end
            
            % Weighting of the image excerpt intensities with the positional and
            % directional weighting factors.
            if j + centralSliceNumber < 2/3*size(image,3)
                weightedIntensities = intensity.*dorsoVentralPositionalWeighting.*...
                    lateralPositionalWeighting.*dorsoVentralDirectionalWeighting.*lateralDirectionalWeighting;
            else
                % Use symmetric dorso-ventral positional weighting for the last
                % one third of slices.
                weightedIntensities = intensity.*dorsoVentralPositionalWeightingEnd.*...
                    lateralPositionalWeighting.*dorsoVentralDirectionalWeighting.*lateralDirectionalWeighting;
            end
            
            % Update positions.
            oldUpperRowPosition = currentUpperRowPosition;
            oldUpperColumnPosition = currentUpperColumnPosition;
            
            % Find the position of the maximum within the image excerpt.
            [maximum,index] = max(weightedIntensities,[],'all','linear');
            [row,column] = ind2sub(size(weightedIntensities,1,2),index);

            % Store the maximum position.
            currentUpperRowPosition = oldUpperRowPosition-dorsoVentralExtent-1+row;
            currentUpperColumnPosition = getCurrentColumnPosition(oldUpperColumnPosition,column,lateralRange);
            
            % Write into segmentation array.
            needleROIs(currentUpperRowPosition,currentUpperColumnPosition,centralSliceNumber+j) = i;
            
            % Store point.
            pointsOfLabelledSegmentation(1,centralSliceNumber+j,i) = currentUpperRowPosition;
            pointsOfLabelledSegmentation(2,centralSliceNumber+j,i) = currentUpperColumnPosition;
            
            % Determine axial segmentation direction.
            upperRowDirection = currentUpperRowPosition-oldUpperRowPosition;
            upperColumnDirection = currentUpperColumnPosition-oldUpperColumnPosition;
        end
    end
end

% Plot the segmentation.

% Thicken the segmentation contour for better visibility. Colour 8-neighbourhood.
se = strel('square',3);
for i = 1:size(image,3)
    needleROIsThick(:,:,i) = imdilate(needleROIs(:,:,i),se);
end

% Convert image and segmentation into rgb format.
[rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(image,needleROIsThick,cmap);
rgbImage = rgbImage*255/max(max(max(max(rgbImage))));
rgbImage = uint8(rgbImage);
   
% Superimpose image and segmentation.
superposition = superimposeImageAndSegmentation(rgbImage,rgbNeedles,0.4);

% Plot superposition.
%figure
%imshow3Dfullcolour(superposition,[],cmap)

% Calculate the sum of intensity values within all voxels belonging to the
% same label.
pointIntensities = zeros(1,numberOfLabels);
% For every label.
for i = 1:numberOfLabels
    % For all points of a label.
    for j = 1:size(pointsOfLabelledSegmentation,2)
        pointIntensities(1,i) = pointIntensities(1,i) + image(pointsOfLabelledSegmentation(1,j,i),...
            pointsOfLabelledSegmentation(2,j,i),j);
    end
end

% Store all labels.
allLabels = 1:numberOfLabels;

% For the central half of all slices, store the segmentation point indices for 
% all labels and the respective summed intensity value as in the former
% alignment arrays to, again, merge identical maxima, this time removing
% doubled segmentations in close proximity.
% For the central half of all slices.
for j = floor(1/4*size(image,3)):floor(3/4*size(image,3))
    totalIntensity = zeros(1,length(allLabels));
    counter = 1;
    % For all labels.
    for i = allLabels
        totalIntensity(1,counter) = pointsOfLabelledSegmentation(1,j,i);
        totalIntensity(2,counter) = pointsOfLabelledSegmentation(2,j,i);
        totalIntensity(3,counter) = pointIntensities(1,i);
        counter = counter+1;
    end
    % Merge doubled segmentations in close proximity.
    [keptMaxima,keptIndices] = mergeIdenticalMaxima(totalIntensity,11);
    keptIndices = sort(keptIndices);
    
    % Update labels discarding labels that have been excluded in the merge process.
    counter = 1;
    oldAllLabels = allLabels;
    clear allLabels;
    for k = 1:length(keptIndices)
        allLabels(counter) = oldAllLabels(1,keptIndices(1,k));
        counter = counter+1;
    end 
end

% Get valid summed intensities.
totalIntensity = keptMaxima(3,:);

% Correct the number of labels.
numberOfLabels = length(allLabels);

% Update point storage arrays as well.
reliabilityOfSegmentation = zeros(13,numberOfLabels);
oldPointsOfLabelledSegmentation = pointsOfLabelledSegmentation;
clear pointsOfLabelledSegmentation
pointsOfLabelledSegmentation = zeros(2,size(oldPointsOfLabelledSegmentation,2),numberOfLabels);
counter = 1;
for i = 1:length(allLabels)
    pointsOfLabelledSegmentation(:,:,counter) = oldPointsOfLabelledSegmentation(:,:,allLabels(i));
    counter = counter +1;
end

%% Assessment of the Reliabilty of Needle Detections.
% Set reliability criteria counter.
criterium = 1;

% ASSESS MULTIPLE DETECTIONS OF THE SAME NEEDLE.
% Use the mergeIdenticalMaxima function to identify detections that come close
% to each other in the outer slices that were not examined in the preceding
% merge process.

counter = 1;
% Select outer slices to assess needles for remaining proximity. Do not consider 
% the innermost slices as the needles terminate earlier.
outerSlices = floor(3/4*size(image,3))+1:1:size(pointsOfLabelledSegmentation,2);

% For the selected slices.
for i = outerSlices
    % Write points and corresponding total intensity into an array to feed
    % into the merge function.
    alignmentReliability = zeros(3,numberOfLabels);
    % For all labels.
    for j = 1:numberOfLabels
        alignmentReliability(1,j) = pointsOfLabelledSegmentation(1,i,j);
        alignmentReliability(2,j) = pointsOfLabelledSegmentation(2,i,j);
        alignmentReliability(3,:) = totalIntensity;
    end
    % Determine segmentations that come to lie within a Euclidean distance
    % of 13 on these outer slices.
    [~,keptMultidetections{counter}] = mergeIdenticalMaxima(alignmentReliability,13);
    counter = counter+1;
end

% Initialise a scoring array containing two for all labels.
multipleDetectionScore = 2*ones(1,numberOfLabels);

% Initialise a cell array to store missing labels.
missingLabelNumbers = {};

% Get all possible labels.
possibleLabels = 1:numberOfLabels;

counter = 1;
% For the multidetections in each slice.
for i = 1:size(keptMultidetections,2)
    % Only if labels are missing.
    if size(keptMultidetections{i},2) < numberOfLabels
        % Find missing label numbers comparing all possible label numbers
        % to those kept in the merge process.
        missingLabelsMask = ~ismember(possibleLabels,keptMultidetections{i});
        missingLabelNumbers{counter} = missingLabelsMask.*possibleLabels; 
    end
    counter = counter+1;
end
missingLabelNumbers = cell2mat(missingLabelNumbers);
missingLabelNumbers = unique(missingLabelNumbers);
missingLabelNumbers(:,~any(missingLabelNumbers,1)) = [];

% Set the multiple detection score zero for missing labels to decline their total
% reliability score.
for i = missingLabelNumbers
    multipleDetectionScore(1,i) = 0;
end

% Write multiple detection score into reliability of segemtation array.
reliabilityOfSegmentation(criterium,:) = multipleDetectionScore;
criterium = criterium+1;

% ASSESS SEGMENTATION INTENSITIES.
% Assess the total intensity within the segmented needles to identify
% spurious detections with very low intensity.

% Calculate the mean label intensity for all labelled segmentations.
meanLabelIntensity = mean(totalIntensity);
% Caluclate the standard deviation of all total intensity values.
stdLabelIntensity = std(totalIntensity);

% Score the intensity of each segmented needle in relation to the mean
% intensity of all segmented needles and write this score into the
% reliability of segmentation array.
reliabilityOfSegmentation(criterium,:) = uint8(totalIntensity>(meanLabelIntensity-0.0*stdLabelIntensity))...
    +uint8(totalIntensity>(meanLabelIntensity-1.2*stdLabelIntensity));
criterium = criterium+1;

% ASSESS TRANSAXIAL CRITERIA.
% Analyse transaxial characteristics of the segmented needles to score the
% credibility of detections.

% Choose slices to consider for transaxial analysis.
numberOfSlices = 11;
outerSlice = floor(numberOfSlices/2);

% Initialise arrays for local grey value range and transaxial mean grey
% value.
localRange = zeros(numberOfSlices,numberOfLabels);
localMean = zeros(numberOfSlices,numberOfLabels);

% Set the dimensions of the square transaxial area to analyse.
halfSquare = 8;

% Select angles for Radon transform anisotropy analysis.
angles = linspace(0,170,18);

% Initialise array for Radon transform results.
radonMaxima = zeros(numberOfSlices,numberOfLabels,length(angles));

sliceCounter = 1;
centralSliceX = zeros(1,numberOfLabels);
% For the selected number of central slices.
for h = linspace(-outerSlice,outerSlice,numberOfSlices)
    % For all labels.
    for i = 1:numberOfLabels
        % Get the needle segmentation point on the current slice.
        pointX = pointsOfLabelledSegmentation(2,needleBoxCentre+h,i);
        pointY = pointsOfLabelledSegmentation(1,needleBoxCentre+h,i);
        % Store pointX on central slice for later weighting of range
        % criterium.
        if h == 0
        centralSliceX(1,i) = pointX;
        end
        % Only if there is a point unequal zero.
        if pointX*pointY > 0
            % Get transaxial image excerpt.
            if pointY-halfSquare <= 0
                neighbourhood = image(1:1:pointY+halfSquare,pointX-halfSquare:1:pointX+...
                    halfSquare,needleBoxCentre+h);
            elseif pointY+halfSquare > size(image,1)
                neighbourhood = image(pointY-halfSquare:1:size(image,1),pointX...
                    -halfSquare:1:pointX+halfSquare,needleBoxCentre+h);
            elseif pointX-halfSquare <= 0
                neighbourhood = image(pointY-halfSquare:1:pointY+halfSquare,...
                    1:1:pointX+halfSquare,needleBoxCentre+h);
            elseif pointX+halfSquare > size(image,2)
                neighbourhood = image(pointY-halfSquare:1:pointY+halfSquare,...
                    pointX-halfSquare:1:size(image,2),needleBoxCentre+h);
            else
                neighbourhood = image(pointY-halfSquare:1:pointY+halfSquare,...
                    pointX-halfSquare:1:pointX+halfSquare,needleBoxCentre+h);
            end
            neighbourhood = squeeze(neighbourhood);
            
            % Set zeros NaN to exclude these values from further calculations.
            neighbourhoodNaN = neighbourhood;
            neighbourhoodNaN(neighbourhood==0) = NaN;
            
            % Determine minimum, maximum and mean intensity of the
            % transaxial image excerpt.
            localMin = min(min(neighbourhoodNaN));
            localMax = max(max(neighbourhoodNaN));
            localMean(sliceCounter,i) = mean(neighbourhood,'all');
            
            % Calculate the local grey value range.
            localRange(sliceCounter,i) = localMax-localMin;
            
            % Assess the local anisotropy using the Radon transform.
            for j = 1:length(angles)
                radonTransform = radon(neighbourhood,angles(j));
                % Find and store the maximum value obtained by the Radon
                % transformation.
                radonMaxima(sliceCounter,i,j) = max(radonTransform);
            end   
        end
    end
    sliceCounter = sliceCounter+1;
end

% Determine the mean local grey value range for every labelled segmentation.
meanLabelLocalRange = mean(localRange,1);
% Calculate the mean local grey value range for all labels.
meanLocalRange = mean(meanLabelLocalRange);
% Calculate the corresponding standard deviation.
stdLocalRange = std(meanLabelLocalRange);
% Only award points to segmentations with sufficiently high grey value range.
reliabilityOfSegmentation(criterium,:) = (meanLabelLocalRange-meanLocalRange) > -0.6*stdLocalRange;
criterium = criterium+1;

% Initialise array to store the anisotropy result in.
anisotropyLabels = zeros(numberOfSlices,numberOfLabels);
% For every slice.
for h = 1:size(radonMaxima,1)
    % For every needle detection.
    for i = 1:size(radonMaxima,2)
        % Comparing all maxima determined in the radon transform for each label and
        % slice, find the lowest and highest maximum among the maxima for different
        % angles.
        projectionMin = min(radonMaxima(h,i,:));
        projectionMax = max(radonMaxima(h,i,:));
        % Calculate the anisotropy for every slice and every label as the 
        % difference between projection maximum and minimum.
        anisotropyLabels(h,i) = projectionMax-projectionMin;
    end
end

% Calculate the mean anisotropy for every label.
meanLabelAnisotropy = mean(anisotropyLabels);
% Calculate the mean anisotropy for all labels.
meanAnisotropy = mean(meanLabelAnisotropy);
% Determine the corresponding standard deviation.
stdAnisotropy = std(meanLabelAnisotropy);

% Award points to segmentations with sufficiently high anisotropy in the
% transaxial image excerpts.
reliabilityOfSegmentation(criterium,:) = (meanLabelAnisotropy-meanAnisotropy) > -0.7*stdAnisotropy;
criterium = criterium+1;

% ASSESS EXTREME NEEDLE PATHS.
% Analyse every segmentation's needle range in rows and columns
% and punish extreme deviations from the mean range.

% Initialise arrays to store row and column index ranges in.
xRange = zeros(1,numberOfLabels);
yRange = zeros(1,numberOfLabels);

% Define weighting factors for the range criterium to account for TRUS
% characteristics where objects of the same size are smeared over a larger
% range of rows and columns in axial view the further away they are from the
% TRUS probe.
weightingRange = ones(1,numberOfLabels);
weightingRange(centralSliceX<2/3*size(image,1)) = 0.75;
weightingRange(centralSliceX<1/3*size(image,1)) = 0.5;

stopSlice = size(pointsOfLabelledSegmentation,2);
% Determine ranges of rows and columns spanned by the segmentations.
for i = 1:numberOfLabels
    xRange(1,i) = max(pointsOfLabelledSegmentation(1,1:1:stopSlice,i))-...
        min(pointsOfLabelledSegmentation(1,1:1:stopSlice,i));
    yRange(1,i) = max(pointsOfLabelledSegmentation(2,1:1:stopSlice,i))-...
        min(pointsOfLabelledSegmentation(2,1:1:stopSlice,i));
end

% Apply weighting factors to the determined x- and y-ranges.
xRange = xRange.*weightingRange;
yRange = yRange.*weightingRange;

% Calculate the mean x-range.
meanXRange = mean(xRange);
% Calculate the corresponding standard deviation.
stdXRange = std(xRange);
% Award points to segmentations that have a moderate x-range.
reliabilityOfSegmentation(criterium,:) = xRange>meanXRange-2.5*stdXRange & xRange<meanXRange+2.5*stdXRange;
criterium = criterium+1;

% Determine the mean y-range.
meanYRange = mean(yRange);
% Caluclate the corresponding standard deviation.
stdYRange = std(yRange);
% Award points to segmentations which have a moderate y-range.
reliabilityOfSegmentation(criterium,:) = yRange>meanYRange-2.5*stdYRange & yRange<meanYRange+2.5*stdYRange;
criterium = criterium+1;

% ASSESS LONGITUDINAL CONTINUITY.
% Analyse the grey value evolvement along the segmentated needles.

% Calculate the change in the mean local grey values from slice to slice
% determined previously in the transaxial analysis for every label.
greyValueDifferences = abs(diff(localMean,1));
% Calculate the mean grey value change from slice to slice for every label.
meanGreyValueDifferences = mean(greyValueDifferences);
% Determine the overall mean grey value change.
meanGreyValueDifference = mean(meanGreyValueDifferences);
% Calculate the corresponding standard deviation.
stdGreyValueDifference = std(meanGreyValueDifferences);

% Find extremely high grey value differences for every needle from slice to
% slice.
largeGreyValueDifferences = greyValueDifferences > meanGreyValueDifference+4.5*stdGreyValueDifference;
% Calculate the number of these huge grey value fluctuations per needle.
greyValueFluctuation = sum(largeGreyValueDifferences,1);

% Do only award points to segmentations with a maximum of 2 huge grey value
% changes.
reliabilityOfSegmentation(criterium,:) = greyValueFluctuation <= 2;
criterium = criterium+1; 

% ASSESS POSITIONAL CRITERIA.
% Analyse the segmented needle's position on the central slice in relation 
% to the image centre, the ROI rim and areas above or below the centre that 
% are likely to contain the highly intense urethra catheter or hyperechoic 
% tissue between prostate and rectum.

% Define a point at slice centre.
sliceCentre = [double(centreX) double(centreY)];

% Initialise arrays for distances to the centre, the rim and the ROI
% ellipse vertices.
distanceToCentre = zeros(1,numberOfLabels);
distanceToRim = zeros(1,numberOfLabels);
distancesToVertices = zeros(1,size(ellipseVertices,1));

% Define the range of columns around the centre in which the urethra
% catheter may likely come to lie.
urethraBoxHalfWidth = int16(floor(size(sagittalImage,2)/7));

% Define a range of columns around the centre in which hyperechoic tissue
% between prostate and rectum is likely to come to lie within the elliptical
% ROI.
prostateBoxHalfWidth = int16(floor(size(sagittalImage,2)/4));

% Initialise score arrays for the urethra catheter and hyperechoic tissue
% between prostate and rectum criteria.
isLikelyUrethra = ones(1,numberOfLabels);
isLikelyBehindProstate = 2*ones(1,numberOfLabels);

% Define where the lower axial image half starts.
lowerImageHalf = size(image,1)-centreX;

% For all labels.
for i = 1:numberOfLabels
    % Get the point indices of the needle on the central slice.
    needle = [pointsOfLabelledSegmentation(1,centralSliceNumber,i),...
        pointsOfLabelledSegmentation(2,centralSliceNumber,i)];
    % Compute the distance from the needle centre to the slice centre.
    distanceToCentre(i) = norm(needle-sliceCentre);
    
    % For all ellipse vertices.
    for j = 1:size(ellipseVertices,1)
        % Find the Euclidean distances from the needle point to every ellipse vertex.
        distancesToVertices(j) = norm(needle-[ellipseVertices(j,1) ellipseVertices(j,2)]);
    end
    % Find the minimum distance to one of the rim points.
    [minimumDistance,~] = min(distancesToVertices);
    distanceToRim(i) = minimumDistance;
    
    % Find needles that are likely false urethra catheter detections
    % because of their positioning in an area where the urehtra typically
    % comes to lie. Punish them by point withdrawal.
    if needle(1) < centreX+0.3*lowerImageHalf && ((needle(2) > centreY-urethraBoxHalfWidth) &&...
            (needle(2) < centreY+urethraBoxHalfWidth))
        isLikelyUrethra(i) = 0;
    end
    
    % Find needles that lie in a box located in the lower area between prostate 
    % and rectum excluding no columns.
    if needle(1) > centreX+0.65*lowerImageHalf
        isLikelyBehindProstate(i) = isLikelyBehindProstate(i)-1;
    end      
    % Find needles within a box between prostate and rectum restricted to
    % more central columns.
    if needle(1) > centreX+0.75*lowerImageHalf && ((needle(2) > centreY-prostateBoxHalfWidth) &&...
            (needle(2) < centreY+prostateBoxHalfWidth))
        isLikelyBehindProstate(i) = isLikelyBehindProstate(i)-1;
    end
end

% Calculate the mean distance of the segmented needles to the image centre.
meanDistanceToCentre = mean(distanceToCentre);
% Calculate the corresponding standard deviation.
stdDistanceToCentre = std(distanceToCentre);

% Award points to detections that are not too close to the image centre in
% a two level graduated system.
reliabilityOfSegmentation(criterium,:) = ((distanceToCentre-meanDistanceToCentre) > -0.7*stdDistanceToCentre) +...
    ((distanceToCentre-meanDistanceToCentre) > -1.2*stdDistanceToCentre);
criterium = criterium+1;

% Caluclate the mean distance of detections to the ellipse vertices.
meanDistanceToRim = mean(distanceToRim);
% Determine the corresponding standard deviation.
stdDistanceToRim = std(distanceToRim);

% Award points to detections that are not too close to the rim in this relative
% distance measure.
reliabilityOfSegmentation(criterium,:) = (distanceToRim-meanDistanceToRim) > -1.0*stdDistanceToRim;
criterium = criterium+1;

% Write the score for potential urethra catheter detections into the reliability
% array.
reliabilityOfSegmentation(criterium,:) = isLikelyUrethra;
criterium = criterium+1;

% Fill in the negative score for lower detections between prostate and rectum.
reliabilityOfSegmentation(criterium,:) = isLikelyBehindProstate;
criterium = criterium+1;

% Initialise array to store the absolute distances of segmentation points
% to the ellipse vertices.
absoluteRimLocalisation = {};
% For all labels.
for i = 1:numberOfLabels
    % For a range of slices.
    for j = floor(1/4*size(image,3)):1:size(pointsOfLabelledSegmentation,2)
        % Get the needle coordinates.
        needle = [pointsOfLabelledSegmentation(1,j,i),pointsOfLabelledSegmentation(2,j,i)];
        % For all ellipse vertices.
        for k = 1:size(ellipseVertices,1)
            % Find the Euclidean distances to the rim points.
            distancesToVertices(k) = norm(needle-[ellipseVertices(k,1) ellipseVertices(k,2)]);
        end
        % Find the minimum absolute distance to a rim point.
        [minimumDistance,~] = min(distancesToVertices);
        % If this minimum absolute distance is less than three, store the
        % label.
        if minimumDistance < 3
            absoluteRimLocalisation{counter} = i;
            counter = counter+1;
        end
    end
end
absoluteRimLocalisation = cell2mat(absoluteRimLocalisation);
absoluteRimLocalisation = unique(absoluteRimLocalisation);

% Initialise the scoring array for absolute distance to rim criterium with
% two.
absoluteProximityToRim = 2*ones(1,numberOfLabels);

% Set score to zero if the absolute distance to the rim is small.
for i = absoluteRimLocalisation 
    absoluteProximityToRim(i) = 0;
end

% Write score for the absolute distance to the rim into the reliability of
% segmentation array.
reliabilityOfSegmentation(criterium,:) = absoluteProximityToRim;
criterium = criterium+1;

% If a needle fails to fulfill the grey value range criterium and the
% isLikelyBehindProstate criterium, subtract another point classifying this
% result combination as an indicator of an erroneous detection.
for i = 1:numberOfLabels
    if reliabilityOfSegmentation(5,i) == 0 && reliabilityOfSegmentation(11,i) < 2
        reliabilityOfSegmentation(criterium,i) = -1;
    end
end

%% Evaluation of the Reliability of Segmentation.
% Set a threshold of 15 for the sum of the reliability score points above 
% which a segmentation is considered to truly represent a needle. The
% maximum achievable reliability count is 17.
reliabilityThreshold = 15;

% Calculate the reliability count.
reliabilityCount = sum(reliabilityOfSegmentation,1);

% Remove the excluded segmentations from the segmentation point storage
% array and total intensity storage array.
counter = 1;
for i = 1:size(reliabilityCount,2)
    if reliabilityCount(1,i) < reliabilityThreshold
        pointsOfLabelledSegmentation(:,:,i) = 0;
        totalIntensity(i) = 0;
    end
end
pointsOfLabelledSegmentation(:,:,~any(pointsOfLabelledSegmentation,[1 2])) = [];
totalIntensity(:,~any(totalIntensity,1)) = [];

% Update the number of labels.
numberOfLabels = size(pointsOfLabelledSegmentation,3);

% Initialise an array to store the maximum row positions of the segmentations in.
maximumRowPositions = zeros(1,numberOfLabels);
% For all labels.
for i = 1:numberOfLabels
    % Get the maximum row position in all slices for every label.
    maximumRowPositions(1,i) = max(pointsOfLabelledSegmentation(1,:,i));
end

% Sort these row positions in descending order.
[maximumRowPositionsSorted,indices] = sort(maximumRowPositions,'descend');

% To assess the data for extreme row positions only consider the lower half
% of row positions in the axial image view.
lowerRowPositions = maximumRowPositionsSorted(1:1:floor(numberOfLabels/2));
indices = indices(1:1:floor(numberOfLabels/2));

% For these selected row positions determine the mean position.
meanLowerRowPositions = mean(lowerRowPositions);
% Calculate the corresponding standard deviation.
stdLowerRowPositions = std(lowerRowPositions);

% Exclude detections with extraordinarily low positions in the axial image
% that are likely to represent the hyperechoic tissue between prostate and
% rectum instead of needles.
mostLikelyBehindProstate = (lowerRowPositions-meanLowerRowPositions) > 1.6*stdLowerRowPositions;

% Get the labels corresponding to the indices to exclude from the storage
% arrays.
indicesToExclude = mostLikelyBehindProstate.*indices;
indicesToExclude(:,~any(indicesToExclude,1)) = [];

% Set positions to exclude to zero.
for i = indicesToExclude
    pointsOfLabelledSegmentation(:,:,i) = 0;
    totalIntensity(1,i) = 0;
end

% Update arrays.
pointsOfLabelledSegmentation(:,:,~any(pointsOfLabelledSegmentation,[1 2])) = [];
totalIntensity(:,~any(totalIntensity,1)) = [];

% Determine the final number of labels.
numberOfLabels = size(pointsOfLabelledSegmentation,3);

% Write into a new segmentation array.
newNeedleROIs = zeros(size(image,1,2,3));
% For every label.
for i = 1:length(totalIntensity)
    % For all slices.
    for j = 1:size(pointsOfLabelledSegmentation,2)
        newNeedleROIs(pointsOfLabelledSegmentation(1,j,i),pointsOfLabelledSegmentation(2,j,i),j) = i;
    end
end

% Thicken the segmentation contour for better visibility. Colour 8-neighbourhood.
se = strel('square',3);
for i = 1:size(image,3)
    newNeedleROIs(:,:,i) = imdilate(newNeedleROIs(:,:,i),se);
end

% Convert image and segmentation into rgb format.
[rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(image,newNeedleROIs,cmap);
rgbImage = rgbImage*255/max(max(max(max(rgbImage))));
rgbImage = uint8(rgbImage);
   
% Superimpose image and segmentation.
superposition = superimposeImageAndSegmentation(rgbImage,rgbNeedles,0.4);

% Plot superposition.
%figure
%imshow3Dfullcolour(superposition,[],cmap)

%% Needle Tip Detection.
% Initialise an array to store all the relevant segmentation characteristics 
% affecting the decision about the needle tip definition.
segmentationCharacteristics = zeros(8,size(pointsOfLabelledSegmentation,2),numberOfLabels);

% Fill in all indices for the points included in the segmentation.
segmentationCharacteristics(1:1:2,:,:) = pointsOfLabelledSegmentation;

% Initialise an array for the local transaxial grey value ranges.
greyValueRanges = zeros(1,size(pointsOfLabelledSegmentation,2));

% Initialise an array to store the grey values having the same 2D position
% as the segmented voxel on the current slice, but lying on the next slice.
greyValueSamePositionNextSlice = zeros(1,size(pointsOfLabelledSegmentation,2)-1);

% For every label.
for i = 1:numberOfLabels
    % For every slice.
    for j = 1:size(pointsOfLabelledSegmentation,2)
        pointX = pointsOfLabelledSegmentation(2,j,i);
        pointY = pointsOfLabelledSegmentation(1,j,i);
        % Store the grey value of every segmentation point.
        segmentationCharacteristics(3,j,i) = image(pointY,pointX,j);
        
        % If this is not the first slice, determine the grey value at the
        % current 2D position but on the next slice.
        if j > 1
            greyValueSamePositionNextSlice(1,j-1) = image(pointY,pointX,j-1);
        end
        
        % Choose a transaxial neighbourhood around the segmented point.
        if pointX*pointY > 0
            if pointY-halfSquare <= 0
                neighbourhood = image(1:1:pointY+halfSquare,pointX-halfSquare:1:pointX+halfSquare,j);
            elseif pointY+halfSquare > size(image,1)
                neighbourhood = image(pointY-halfSquare:1:size(image,1),pointX-halfSquare:1:pointX+halfSquare,j);
            elseif pointX-halfSquare <= 0
                neighbourhood = image(pointY-halfSquare:1:pointY+halfSquare,1:1:pointX+halfSquare,j);
            elseif pointX+halfSquare > size(image,2)
                neighbourhood = image(pointY-halfSquare:1:pointY+halfSquare,pointX-halfSquare:1:size(image,2),j);
            else
                neighbourhood = image(pointY-halfSquare:1:pointY+halfSquare,pointX-halfSquare:1:pointX+halfSquare,j);
            end
            neighbourhood = squeeze(neighbourhood);
            neighbourhoodNaN = neighbourhood;
            % Set zero entries to NaN to avoid distortion of the grey value
            % range determination.
            neighbourhoodNaN(neighbourhood==0) = NaN;
            
            % Determine the minimum and the maximum grey value and the
            % corresponding grey value range.
            localMin = min(min(neighbourhoodNaN));
            localMax = max(max(neighbourhoodNaN));
            greyValueRanges(1,j) = localMax-localMin;
        end
    end
    
    % Determine the change in x-Position of the segmentation points from slice to slice.
    differencesX = diff(segmentationCharacteristics(1,:,i));
    % Determine the change in intensity of the segmentation points from slice to slice.
    differencesIntensity = diff(segmentationCharacteristics(3,:,i));
    % Determine the change in the grey value range around the segmentation
    % points from slice to slice.
    differencesGreyValueRange = diff(greyValueRanges,1,2);
    
    % Fill all these into the segmentation characteristics array.
    segmentationCharacteristics(4,2:1:size(segmentationCharacteristics,2),i) = differencesX;
    segmentationCharacteristics(5,2:1:size(segmentationCharacteristics,2),i) = differencesIntensity; 
    segmentationCharacteristics(6,2:1:size(segmentationCharacteristics,2),i) = differencesGreyValueRange;
    % Fill in the grey value changes of the segmentation point and the
    % point at the same 2D position on the next slice.
    segmentationCharacteristics(7,2:1:size(segmentationCharacteristics,2),i) = ...
        segmentationCharacteristics(3,2:1:size(pointsOfLabelledSegmentation,2),i)-greyValueSamePositionNextSlice;
end

% Store an unnormalised versioin of the segmentation characteristics array. 
segmentationCharacteristicsUnnormalised = segmentationCharacteristics;

% Use the inner half of all slices to find the maxima for normalisation.
slicesToNorm = 2;

% Maximum-normalisation of the four last segmentation characteristics.
outerSliceLimit = round(size(segmentationCharacteristics,2)/slicesToNorm);
for i = 1:numberOfLabels
    norm4 = max(abs(segmentationCharacteristics(4,1:1:outerSliceLimit,i)));
    norm5 = max(abs(segmentationCharacteristics(5,1:1:outerSliceLimit,i)));
    norm6 = max(abs(segmentationCharacteristics(6,1:1:outerSliceLimit,i)));
    norm7 = max(abs(segmentationCharacteristics(7,1:1:outerSliceLimit,i)));
    segmentationCharacteristics(4,:,i) = segmentationCharacteristics(4,:,i)/norm4;
    segmentationCharacteristics(5,:,i) = segmentationCharacteristics(5,:,i)/norm5;
    segmentationCharacteristics(6,:,i) = segmentationCharacteristics(6,:,i)/norm6;
    segmentationCharacteristics(7,:,i) = segmentationCharacteristics(7,:,i)/norm7;
end

% Set all infite values to NaN.
segmentationCharacteristics(segmentationCharacteristics == inf) = NaN;
segmentationCharacteristics(segmentationCharacteristics == -inf) = NaN;

% Use the absolute values for the x-positional change characteristic in the
% normalised data array.
segmentationCharacteristics(4,:,:) = abs(segmentationCharacteristics(4,:,:));

% Calculate the sum of the four normalised criteria and store it in the
% last line of the segmentation characteristics array.
numberOfSummedCriteria = 4;
segmentationCharacteristics(8,:,:) = round(sum(segmentationCharacteristics([4,5,6,7],:,:),1,'omitnan'));

% Determine the mean positional change in x direction for each dataset to 
% sort all datasets into groups which are assigned different limits for
% kinks and intensity to terminate the needle segmentation at the needle tips.
meanXPositionalChange = mean(abs(segmentationCharacteristicsUnnormalised(4,:,:)),'all');

% Set group-specific needle kink and intensity thresholds.
if meanXPositionalChange < 0.85
    kinkLimit = 3;
    intensityLimit = 0.8;
elseif meanXPositionalChange > 1.23
    kinkLimit = 6;
    intensityLimit = 0.6;
else
    kinkLimit = 4;
    intensityLimit = 0.75;
end

% For all labels.
for i = 1:numberOfLabels
    % Use the highest sum of all four normalised criteria to decide upon the needle
    % tip position.
    [maximum,determinedSlice] = max(segmentationCharacteristics(8,1:1:outerSliceLimit,i));
    % Stop later if the determined slice number is very big and there is another maximum
    % with the same value further in the patient.
    [sameMaximum,indexFurtherIn] = max(segmentationCharacteristics(8,determinedSlice:-1:1,i));
    if determinedSlice > round(size(segmentationCharacteristics,2)/(slicesToNorm+2)) && sameMaximum == maximum
        determinedSlice = determinedSlice-indexFurtherIn+1;
    end
    % Stop earlier if there is a huge intensity drop within the
    % segentation.
    bigIntensityLoss = segmentationCharacteristicsUnnormalised(3,determinedSlice:1:outerSliceLimit,i)...
        <= intensityLimit*segmentationCharacteristicsUnnormalised(3,determinedSlice+1:1:outerSliceLimit+1,i);
    if sum(bigIntensityLoss) > 0
        indices = find(bigIntensityLoss);
        % Do only stop earlier if this big intensity loss is not just a
        % small gap.
        if determinedSlice>2 && segmentationCharacteristicsUnnormalised(3,determinedSlice+max(indices)-2,i)...
                < 0.80*segmentationCharacteristicsUnnormalised(3,determinedSlice+max(indices),i)
            determinedSlice = determinedSlice+max(indices);
            % There is no gap test possible for determined slices close to
            % the stack end.
        elseif determinedSlice <= 2
            determinedSlice = determinedSlice+max(indices);
        end
    end

    % Terminate the needle later if its positional change is small and the
    % corresponding grey value change is more dominant at the second highest 
    % sum metric value further in the patient.
    if abs(segmentationCharacteristicsUnnormalised(4,determinedSlice,i)) <= 2
        % Find the second highest sum metric and use it if the difference to 
        % the initially indentified maximum sum metric is just one.
        [secondMaximum,secondIndex] = max(segmentationCharacteristics(8,determinedSlice-1:-1:1,i));
        if secondMaximum+1 >= maximum && segmentationCharacteristicsUnnormalised(5,determinedSlice-secondIndex,i)...
                > segmentationCharacteristicsUnnormalised(5,determinedSlice,i)
            determinedSlice = determinedSlice-secondIndex;
        end
    end
    
    % Stop earlier if the positional change along the contour is too big.
    bigPositionalChanges = abs(segmentationCharacteristicsUnnormalised...
        (4,determinedSlice:1:outerSliceLimit,i)) >= kinkLimit;
    if sum(bigPositionalChanges) > 0
        indices = find(bigPositionalChanges);
        % Stop at the innermost big positional change.
        determinedSlice = determinedSlice+max(indices)-1;
    end
    
    % Stop later at the second highest sum metric if the needle is terminated 
    % very early and  its grey value recovers thereby accepting gaps of two
    % slices.
    if determinedSlice > round(size(segmentationCharacteristics,2)/(slicesToNorm+1))...
            && (segmentationCharacteristicsUnnormalised(3,determinedSlice-3,i)...
            >= 0.8*segmentationCharacteristicsUnnormalised(3,determinedSlice,i) || ...
            segmentationCharacteristicsUnnormalised(3,determinedSlice-2,i) >= ...
            0.8*segmentationCharacteristicsUnnormalised(3,determinedSlice,i))
        [secondMaximum,secondIndex] = max(segmentationCharacteristics(8,determinedSlice-1:-1:1,i));
        determinedSlice = determinedSlice-secondIndex;
    end  
    segmentationCharacteristics(:,determinedSlice-1:-1:1,i) = NaN;
end
 
% Update the segmentation points storage array taking into account the
% needle tip definition.
pointsOfLabelledSegmentation = segmentationCharacteristics(1:1:2,:,:);

% Write the complete segmentation into a new segmentation array.
newNeedleROIs = zeros(size(image,1,2,3));
% For every label.
for i = 1:numberOfLabels
    % For all slices.
    for j = 1:size(pointsOfLabelledSegmentation,2)
        if isnan(pointsOfLabelledSegmentation(1,j,i)) == false
            newNeedleROIs(pointsOfLabelledSegmentation(1,j,i),pointsOfLabelledSegmentation(2,j,i),j) = i;
        end
    end
end

% Thicken the segmentation contour for better visibility. Colour 8-neighbourhood.
se = strel('square',3);
for i = 1:size(image,3)
    newNeedleROIs(:,:,i) = imdilate(newNeedleROIs(:,:,i),se);
end

% Produce maximally perceptually distinguishable colourmap containing as
% many colours as there are labels.
%cmap = distinguishable_colors(numberOfLabels);

% Remove the initially appended slice.
image = image(:,:,1:1:size(originalImage,3));
newNeedleROIs = newNeedleROIs(:,:,1:1:size(originalImage,3));

% Convert image and segmentation into rgb format.
[rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(image,newNeedleROIs,cmap);
rgbImage = rgbImage*255/max(max(max(max(rgbImage))));
rgbImage = uint8(rgbImage);
   
% Superimpose image and segmentation.
superposition = superimposeImageAndSegmentation(rgbImage,rgbNeedles,0.4);

% Plot superposition.
%figure
%imshow3Dfullcolour(superposition,[],cmap)

%% Restoration of the Initial Image Dimensions.
% Use the voxel indices as initially determined for cropping of the
% original image to write a segmentation array covering the original image
% dimensions.
xMin = floor(xMin);
yMin = floor(yMin);
needleSegmentationFullSize = zeros(size(originalImage,1,2,3));
needleSegmentationFullSize(yMin:1:yMin+height,xMin:1:xMin+width,:) = newNeedleROIs;

% Correct pointsOfLabelledSegmentation for the new full image dimensions.
pointsOfLabelledSegmentationOriginal(1,:,:) = pointsOfLabelledSegmentation(1,:,:)+yMin;
pointsOfLabelledSegmentationOriginal(2,:,:) = pointsOfLabelledSegmentation(2,:,:)+xMin;

% Thicken the segmentation contour for better visibility. Colour 8-neighbourhood.
% se = strel('square',3);
% for i = 1:size(originalImage,3)
%     needleSegmentationFullSize(:,:,i) = imdilate(needleSegmentationFullSize(:,:,i),se);
% end
% 
% Convert image and segmentation into rgb format.
% [rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(originalImage,needleSegmentationFullSize,cmap);
%    
% Superimpose image and segmentation.
% superposition = superimposeImageAndSegmentation(rgbImage,rgbNeedles,0.4);
% 
% Plot superposition.
% figure
% imshow3Dfullcolour(superposition,[],cmap)

%% Application of the Physical Needle Model.

% Pixel size on a slice in mm.
pixelSizeOnSlice = 0.1412; 
% Needle outer radius in mm.
needleOuterRadius = 1.473/2; 

% Define the origin of the TRUS image acquisition for all slices in the DICOM 
% patient-based coordinate system.
originTRUS = [(size(originalImage,1)-1)*pixelSizeOnSlice,(round(size(originalImage,2)/2)-1)*pixelSizeOnSlice];

% Initialise an array to store the segmentation points corrected with the
% physical needle model in. This correction is applied in the axial plane
% for every slice.
pointsOfPhysicalSegmentation = NaN(size(pointsOfLabelledSegmentationOriginal,1,2,3));
% For every label.
for i = 1:numberOfLabels
    % For every slice.
    for j = 1:size(pointsOfLabelledSegmentation,2)
        % If the needle was found on that slice.
        if isnan(pointsOfLabelledSegmentationOriginal(1,j,i)) == false
            % Get the segmentation point and calculate its distance to the
            % origin of the DICOM patient-based coordinate system.
            point = [(pointsOfLabelledSegmentationOriginal(1,j,i)-1)*pixelSizeOnSlice,...
                (pointsOfLabelledSegmentationOriginal(2,j,i)-1)*pixelSizeOnSlice];
            
            % Calculate the length of the hypotenuse in mm.
            hypotenuseLength = sqrt((point(1)-originTRUS(1))^2+(point(2)-originTRUS(2))^2);
            
            % Calculate the total distance of the needle centre from
            % the TRUS origin.
            totalLength = hypotenuseLength+needleOuterRadius;
            
            % Determine the voxel indices for the corrected needle
            % positions.
            needleCentrePositionX = round((originTRUS(1)+totalLength/hypotenuseLength*...
                (point(1)-originTRUS(1)))/pixelSizeOnSlice);
            needleCentrePositionY = round((originTRUS(2)+totalLength/hypotenuseLength*...
                (point(2)-originTRUS(2)))/pixelSizeOnSlice);
            
            % Store the transaxially corrected segmentation points.
            pointsOfPhysicalSegmentation(1,j,i) = needleCentrePositionX;
            pointsOfPhysicalSegmentation(2,j,i) = needleCentrePositionY;
        end
    end
end

% Next, correct the tip position from the TRUS-visible tip position to the
% Vitesse-defined needle end which ignores the dead space at the needle
% tip. Therefore, 4 mm of dead space need to be subtracted from the visibly
% detected needle end.

% Convert the transaxially corrected segmentation points into DICOM
% coordinates.
deadEnd = 4;
pointsOfPhysicalSegmentationCoordinates = {};
% For all labels.
for i = 1:numberOfLabels
    counter = 1;
    % For all slices.
    for j = 1:size(pointsOfPhysicalSegmentation,2)
        if isnan(pointsOfPhysicalSegmentation(1,j,i)) == false
            pointsOfPhysicalSegmentationCoordinates{i}(1,counter) = pointsOfPhysicalSegmentation(1,j,i);
            pointsOfPhysicalSegmentationCoordinates{i}(2,counter) = pointsOfPhysicalSegmentation(2,j,i);
            pointsOfPhysicalSegmentationCoordinates{i}(3,counter) = j;
            counter = counter+1;
        end
    end
    pointsOfPhysicalSegmentationCoordinates{i} = voxels2coordinates(info,pointsOfPhysicalSegmentationCoordinates{i});
    
    % Get the position of the visible needle tip.
    visibleTip = pointsOfPhysicalSegmentationCoordinates{i}(:,1);
    
    % Find the first needle point behind the 4mm shortened needle end. 
    firstIndex = find((pointsOfPhysicalSegmentationCoordinates{i}(3,:) < (visibleTip(3,1)-deadEnd)),1);
    
    % Interpolate new tip position between the visible tip and the next
    % point available on the shortened needle.
    interpolationVector = visibleTip-pointsOfPhysicalSegmentationCoordinates{i}(:,firstIndex);
    interpolationFactor = (interpolationVector(3,1)-deadEnd)/interpolationVector(3,1);
    
    % Shorten the needle end by 4mm to obtain a needle end using the
    % Vitesse needle end definition.
    tipVitesse = pointsOfPhysicalSegmentationCoordinates{i}(:,firstIndex)+interpolationFactor*interpolationVector;
   
    % For all slices of the final segmenation.
    for k = 1:size(pointsOfPhysicalSegmentationCoordinates{i},2)-(firstIndex-2)
        if k == 1
            % Store the new needle tip.
            finalNeedlePointsCoordinates{i}(:,k) = tipVitesse;
        else
            % Insert the remaining needle segmentation points behind the
            % newly defined endpoint.
            finalNeedlePointsCoordinates{i}(:,k) = pointsOfPhysicalSegmentationCoordinates{i}...
                (:,firstIndex+(k-(firstIndex-1)));
        end
    end
end

% Convert the finally identified segmentation points back to voxels to plot the
% contours.
% For all labels.
for i = 1:numberOfLabels
    finalNeedlePointsVoxels{i} = coordinates2voxels(info,finalNeedlePointsCoordinates{i});
end

% Write the segmentation incorporating the physical needle model into a new segmentation array.
physicalNeedleROIs = zeros(size(originalImage,1,2,3));
% For every label.
for i = 1:numberOfLabels
    % For all slices.
    for j = 1:size(finalNeedlePointsVoxels{i},2)
        physicalNeedleROIs(finalNeedlePointsVoxels{i}(1,j),finalNeedlePointsVoxels{i}(2,j),...
            finalNeedlePointsVoxels{i}(3,j)) = i;
    end
end

% Thicken the segmentation contour for better visibility. Colour 8-neighbourhood.
se = strel('square',3);
for i = 1:size(image,3)
    physicalNeedleROIs(:,:,i) = imdilate(physicalNeedleROIs(:,:,i),se);
end

% Convert image and segmentation into rgb format.
[rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(originalImage,physicalNeedleROIs,cmap);
   
% Superimpose image and segmentation.
superposition = superimposeImageAndSegmentation(rgbImage,rgbNeedles,0.4);

% Plot superposition.
%figure
%imshow3Dfullcolour(superposition,[],cmap)

% Display the same figure in ROI size.
if sliceAppended == true
    physicalNeedleROIsCropped = imcrop3(physicalNeedleROIs,[xMin yMin zMin width height depth-1]);
else
    physicalNeedleROIsCropped = imcrop3(physicalNeedleROIs,[xMin yMin zMin width height depth]);
end

% Convert image and segmentation into rgb format.
[rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(image,physicalNeedleROIsCropped,cmap);
rgbImage = uint8(rgbImage);
   
% Superimpose image and segmentation.
superposition = superimposeImageAndSegmentation(rgbImage,rgbNeedles,0.4);

% Plot superposition.
figure
imshow3Dfullcolour(superposition,[],cmap)

%% Spline Interpolation of Needle Segmentations
% splineSegmentation = zeros(size(image,1,2,3));
% % For every label.
% for i = allLabels
%     x = pointsOfLabelledSegmentation(1,:,i);
%     y = pointsOfLabelledSegmentation(2,:,i);
%     z = 1:1:size(image,3);
%     % use every 4th point only to test.
%     x = x(1:3:size(image,3));
%     y = y(1:3:size(image,3));
%     z = z(1:3:size(image,3));
%     [pp,p]=csaps(z,[x;y]);
%     znew = 1:1:size(image,3);
%     value=fnval(pp,znew);
%     value = round(value);
%     for j = 1:size(value,2)
%         splineSegmentation(value(1,j),value(2,j),j) = i;
%     end
% end
%% Time Measurement.
% Get time elapsed.
%toc

%% Local Functions.

function [imageAllSlices,infoAllSlices] = loadTRUS(path)
% LOADTRUS - Load Transrectal Ultrasound (TRUS) DICOM data
% The function loads, reads and stores TRUS DICOM data and the associated 
% DICOM information from a folder. TRUS data need to be stored labelled as 
% Modality 'MR'.
%
% Inputs:
%    path - path of the folder containing the image data
%
% Outputs:
%    imageAllSlices - images of all acquired TRUS slices in one 4D dataset
%                     (8 bit)
%    infoAllSlices - DICOM information for all acquired slices in one cell
%                    array
%
% Other m-files required: none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% April 2020

%------------- CODE --------------

    % Get general metadata from ultrasound DICOM series.
    info = dicominfo([path '/MR001.dcm']);

    % Get all ultrasound files (labeled as MR).
    files = dir([path '/MR*.dcm']);

    % Get general image parameters from DICOM info.
    nRows = info.Rows;
    nCols = info.Columns;
    nPlanes = info.SamplesPerPixel;
    nFrames = size(files,1); 

    % Read and save image and DICOM info.
    imageAllSlices = repmat(uint8(0), [nRows, nCols, nPlanes, nFrames]);
    infoAllSlices = cell(1,length(files));
    
    for i = 1:nFrames
        % Get current filename.
        name = [path '/' vertcat(files(i).name)];
        
        % Get ultrasound DICOM metadata.
        infoAllSlices{i} = dicominfo(name);
        
        % Get image data.
        imageAllSlices(:,:,:,i) = dicomread(name);   
    end
end

function [newAlignment] = findLocalMaxima2D(alignment,needleZAngle,i,j,centre,needleStorage)
% FINDLOCALMAXIMA2D - Find local maxima in 2D space
% The function finds local maxima in 2D space. The minimum separation of
% extrema along each of the two directions and the maximum number of 
% identified extrema in each row or column can be selected in lines 48 and 
% 49. Plotting of the 2D distribution and its identified local maxima can 
% be enabled in the code.
%
% Inputs:
%    alignment - array with variable number of rows containing cross correlation  
%                parameters and results: index i, index j, function value, 
%                reference to needle rotational configuration (facultative)
%    i - first index of the image patch origin 
%    j - second index of the image patch origin
%    needleStorage - storage of array reference index to trace correct
%                    needle rotational configuration (boolean)
% Outputs:
%    newAlignment - alignment array containing alignment data of identified
%                   local maxima
%
% Other m-files required: 
% none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% May 2020

%------------- CODE --------------
 
    maxDistribution = zeros(i+centre,j+centre);
    if needleStorage == true
        storageCell = zeros(i+centre,j+centre);
    end
    % Create 2D ij-space (image patch centre indices) containing the
    % cross correlation results.
    for a = 1:size(alignment,2)
        maxDistribution(alignment(1,a),alignment(2,a)) = alignment(3,a);
        % Store reference index for needle configuration.
        if needleStorage == true
            storageCell(alignment(1,a),alignment(2,a)) = alignment(4,a);
        end
    end

    % Plot distribution of cross correlation values in ij-space.
    %figure
    %imshow(maxDistribution,[min(min(min(maxDistribution))) max(max(max(maxDistribution)))])

    % Determine local extrema in i- and j- direction.
    extremaI = islocalmax(maxDistribution,2,'MinSeparation',11,'MaxNumExtrema',6);
    extremaJ = islocalmax(maxDistribution,1,'MinSeparation',11,'MaxNumExtrema',6);
    % Get common local extrema.
    extrema = extremaI.*extremaJ.*maxDistribution;

    % Plot identified local extrema.
    %figure
    %imshow(extrema)

    % Extract local extrema and their origin and reference parameters as 
    % new matrix of cross correlation information.
    lengthExtrema = length(nonzeros(extrema));
    newAlignment = zeros(size(alignment,1),lengthExtrema);
    extremaCounter = 1;
    for b = centre:i+centre
        for c = centre:j+centre
            if extrema(b,c)~=0
                newAlignment(1,extremaCounter) = b;
                newAlignment(2,extremaCounter) = c;
                newAlignment(3,extremaCounter) = extrema(b,c);
                if needleStorage == true
                    newAlignment(4,extremaCounter) = storageCell(b,c);
                    newAlignment(5,extremaCounter) = needleZAngle{storageCell(b,c)};
                end
                extremaCounter = extremaCounter+1;
            end
        end
    end
end

function [newAlignment,determinedIndices] = mergeIdenticalMaxima(alignment,distance)
% MERGEIDENTICALMAXIMA - Find local maximum of closely adjacent local maxima
% The function identifies closely adjacent local maxima and returns the 
% alignment array containing only the local maxima of these local maxima 
% clusters and non-clustering local maxima.
%
% Inputs:
%    alignment - array with variable number of rows containing cross correlation  
%                parameters and results: index i, index j, function value, 
%                reference to needle rotational configuration (facultative)
%    distance -  threshold value for Euclidean distance below which local
%                maxima are classified as being representations of the same
%                local maximum
% Outputs:
%    newAlignment - alignment array with merged local maxima
%    determinedIndices - original indices that were kept during the merge
%                        process
%
% Other m-files required: 
% none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% May 2020

%------------- CODE --------------

    % Identify closely adjacent entries in alignment array based on Euclidean
    % distance.
    pairCounter = 1; 
    noPairs = linspace(1,size(alignment,2),size(alignment,2));
    for m = 1:size(alignment,2)
        for n = m+1:size(alignment,2)
            calculatedNorm = norm([alignment(1,m) alignment(2,m)]-[alignment(1,n) alignment(2,n)]);
            % Treat and save two identified local maxima as a pair.
            if calculatedNorm < distance
                pairs{pairCounter} = [m,n];
                % Remove paired local maxima from index array to keep track of
                % unpaired alignment entries.
                noPairs(noPairs == m) = 0;
                noPairs(noPairs == n) = 0;
                pairCounter = pairCounter+1;
            end
        end
    end
    if exist('pairs','var') == 1
        % Determine clusters of local maxima sorting connected pairs.
        maxCounter = 1;
        % For all local maxima pairs.
        for i = 1:length(pairs)
            % If there is already a pair sorted into a cluster.
            if exist('sameMax','var') == 1
                paired = false;
                % Identify whether one of the pair components is already
                % sorted into a cluster, if so, add pair numbers to cluster.
                for j = 1:size(sameMax,2)
                    sameElements = ismember(pairs{i},sameMax{j});
                    if sum(sameElements) > 0
                        sameMax{j}(size(sameMax{j},2)+1) = pairs{i}(1);
                        sameMax{j}(size(sameMax{j},2)+1) = pairs{i}(2);
                        % Register that at least one pair was found.
                        paired = true;
                    end
                end
                % If no matching cluster was identified. Define as new cluster
                % centre.
                if paired == false
                    sameMax{maxCounter} = pairs{i};
                    maxCounter = maxCounter+1;
                end
            % If this is the very first pair to be sorted into the clusters, save
            % as cluster centre.
            else 
                sameMax{maxCounter} = pairs{i};
                maxCounter = maxCounter+1;
            end
        end
        % Remove number doublings in clusters.
        determinedMaxima = zeros(1,size(sameMax,2));
        for i = 1:size(sameMax,2)
            sameMax{i} = unique(sameMax{i});
            toCompare = zeros(2,size(sameMax{i},2));
            % Find the local maximum within each cluster of local maxima.
            for j = 1:size(sameMax{i},2)
                toCompare(1,j) = sameMax{i}(j);
                toCompare(2,j) = alignment(3,sameMax{i}(j));
                [~,index] = max(toCompare(2,:));
                determinedMaxima(i) = toCompare(1,index);
            end
        end
        % Store indices of all unpaired entries in alignment array.
        noPairs(:,~any(noPairs,1)) = [];
        % Combine unpaired indices with indices of determined local maxima.
        determinedIndices = cat(2,noPairs,determinedMaxima);

        % Write new alignment array with merged local maxima.
        newAlignment = zeros(size(alignment,1),size(determinedIndices,2));
        for i = 1:size(determinedIndices,2)
            newAlignment(:,i) = alignment(:,determinedIndices(i));
        end
    else
        newAlignment = alignment;
        determinedIndices = 1:size(alignment,2);
    end
end

function [rgbImage,rgbNeedles] = imageAndLabelledNeedles2rgb(greyValueImage,labelledNeedles,cmap)
% IMAGEANDLABELLEDNEEDLES2RGB - convert image and segmentation data to rgb data
% This function converts a 3D grey value image with its corresponding
% labelled 3D segmentation into rgb data.
%
% Inputs:
%    greyValueImage - 3D grey value image to be converted into rgb
%    labelledNeedles - 3D labelled ROI data to be converted into rgb
%    cmap - colourmap to be used
%
% Outputs:
%    rgbImage - image converted to rgb
%    rgbNeedles - ROI data converted to rgb
%
% Other m-files required: 
% label2rgb3d
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% June 2020

%------------- CODE --------------

    % Convert labelled grey value segmentation data into rgb.
    rgbNeedles = label2rgb3d(labelledNeedles,cmap,[0,0,0]);
    rgbNeedles = 255*uint8(rgbNeedles);

    % Convert grey value image to rgb.
    rgbImage = cat(4,greyValueImage,greyValueImage,greyValueImage);
end

function rgb3d=label2rgb3d(varargin)
[label,map,zerocolor,order,fcnflag] = parse_inputs(varargin{:});
% label = 3d image with labels
% map = specified color map
% Source:
% payel ghosh (2020). label2rgb3D (https://www.mathworks.com/matlabcentral/fileexchange/8355-label2rgb3d), 
% MATLAB Central File Exchange. Retrieved June 16, 2020; adapted by LFAnt, 
% 07.12.2019.
%==================================
    numregion = double(max(label(:)));

    % If MAP is a function, evaluate it. Make sure that the evaluated function
    % returns a valid colormap.
    if fcnflag == 1
        cmap = feval(map, numregion);
        if ~isreal(cmap) || any(cmap(:) > 1) || any(cmap(:) < 0) || ...
            ~isequal(size(cmap,2),3) || size(cmap,1) < 1
            eid = sprintf('Images:%s:functionReturnsInvalidColormap',mfilename);
            msg = 'function handle MAP must return a n x 3 colormap array';
            error(eid,'%s',msg);
        end
    else
        cmap = map;
    end

    % If ORDER is set to 'shuffle', save original state. The SHUFFLE keyword
    % uses the same state every time it is called. After shuffling is completed,
    % go back to original state.
    if isequal(order,'shuffle')
        S = rand('state');
        rand('state', 0);
        index = randperm(numregion);
        cmap = cmap(index,:,:);
        rand('state', S);
    end

    % Issue a warning if the zerocolor (boundary color) matches the color of one
    % of the regions.
    for i=1:numregion
        if isequal(zerocolor,cmap(i,:))
            wid= sprintf('Images:%s:zerocolorSameAsRegionColor',mfilename);
            msg= sprintf('Region number %d has the same color as the ZEROCOLOR.',i);
            warning(wid,'%s',msg);
        end
    end
    cmap = [zerocolor;cmap];

    % if label is of type double, need to pass 'label + 1' into IND2RGB.
    % IND2RGB does not like double arrays containing zero values.
    if isa(label, 'double')
        rgb3d = ind2rgb3d(label + 1, cmap);
    else
        rgb3d = ind2rgb3d(label, cmap);
    end
end

%======================================================================
%======================================================================
function [L, Map, Zerocolor, Order, Fcnflag] = parse_inputs(varargin)
% L label 3d matrix: matrix containing non-negative values.
% Map colormap: name of standard colormap, user-defined map, function
% handle.
% Zerocolor RGB triple or Colorspec
% Order keyword if specified: 'shuffle' or 'noshuffle'
% Fcnflag flag to indicating that Map is a function

    valid_order = {'shuffle', 'noshuffle'};
    narginchk(1,4);

    % set defaults
    L = varargin{1};
    Map = 'jet';
    Zerocolor = [1 1 1];
    Order = 'noshuffle';
    Fcnflag = 0;

    % parse inputs
    if nargin > 1
        Map = varargin{2};
    end
    if nargin > 2
        Zerocolor = varargin{3};
    end
    if nargin > 3
        Order = varargin{4};
    end

    % error checking for L
    validateattributes(L,{'numeric' 'logical'}, ...
    {'real' 'nonsparse' 'finite' 'nonnegative' 'integer'}, ...
    mfilename,'L',1);

    % error checking for Map
    [fcn, fcnchk_msg] = fcnchk(Map);
    if isempty(fcnchk_msg)
        Map = fcn;
        Fcnflag = 1;
    else
        if isnumeric(Map)
            if ~isreal(Map) || any(Map(:) > 1) || any(Map(:) < 0) || ...
                ~isequal(size(Map,2), 3) || size(Map,1) < 1
                eid = sprintf('Images:%s:invalidColormap',mfilename);
                msg = 'Invalid entry for MAP.';
                error(eid,'%s',msg);
            end
        else
            eid = sprintf('Images:%s:invalidFunctionforMAP',mfilename);
            error(eid,'%s',fcnchk_msg);
        end
    end

    % error checking for Zerocolor
    if ~ischar(Zerocolor)
        % check if Zerocolor is a RGB triple
        if ~isreal(Zerocolor) || ~isequal(size(Zerocolor),[1 3]) || ...
            any(Zerocolor> 1) || any(Zerocolor < 0)
            eid = sprintf('Images:%s:invalidZerocolor',mfilename);
            msg = 'Invalid RGB triple entry for ZEROCOLOR.';
            error(eid,'%s',msg);
        end
    else
        assert(isscalar(Zerocolor) && ismember(Zerocolor, 'krgybmcw'), 'Only one letter color names allowed for now.');
        Zerocolor = bitget(find('krgybmcw'==Zerocolor)-1,1:3);
    end

    % error checking for Order
    idx = strmatch(lower(Order), valid_order);
    eid = sprintf('Images:%s:invalidEntryForOrder',mfilename);
    if isempty(idx)
        msg = 'Valid entries for ORDER are ''shuffle'' or ''noshuffle''.';
        error(eid,'%s',msg);
    elseif length(idx) > 1
        msg = sprintf('Ambiguous string for ORDER: %s.', Order);
        error(eid,'%s',msg);
    else
        Order = valid_order{idx};
    end
end

%================================================================
%=================================================================
function [rout,g,b] = ind2rgb3d(a,cm)
%IND2RGB Convert indexed image to RGB image.
% RGB = IND2RGB(X,MAP) converts the 3d matrix X and corresponding
% colormap MAP to RGB (truecolor) format.
%
% Class Support
% -------------
% X can be of class uint8, uint16, or double. RGB is an
% M-by-N-by-3 array of class double.
%
% See also IND2GRAY, RGB2IND (in the Image Processing Toolbox).

    if ~isa(a, 'double')
        a = double(a)+1; % Switch to one based indexing
    end

    error(nargchk(2,2,nargin));

    % Make sure A is in the range from 1 to size(cm,1)
    a = max(1,min(a,size(cm,1)));

    % Extract r,g,b components
    r = zeros(size(a)); r(:) = cm(a,1);
    g = zeros(size(a)); g(:) = cm(a,2);
    b = zeros(size(a)); b(:) = cm(a,3);

    if nargout==3
        rout = r;
    else
        rout = zeros([size(r),3]);
        rout(:,:,:,1) = r;
        rout(:,:,:,2) = g;
        rout(:,:,:,3) = b;
    end 
end

function [superposition] = superimposeImageAndSegmentation(image,segmentation,alpha)
% SUPERIMPOSEIMAGEANDSEGMENTATION - superimpose image and segmentation
% This function superimposes the image and the segmentation returning the
% weighted combined dataset. Datatypes of both images have to be the same.
%
% Inputs:
%    image - image to be superimposed
%    segmentation - segmentation to be superimposed
%    alpha - weighting parameter for the image
%
% Outputs:
%    superposition - weighted superimposed dataset
%
% Other m-files required: none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% April 2020

%------------- CODE --------------
    superposition = alpha*image+(1-alpha)*segmentation;
end

function [intensity] = getIntensities(image,currentRowPosition,dorsoVentralExtent,...
    currentColumnPosition,lateralExtent,centralSliceNumber,j)
% GETINTENSITIES - Get image intensities in an image patch.
% This function extracts an image patch which is defined by its centre
% position and a range of pixels.
%
% Inputs:
%    image - axial image
%    currentRowPosition - the labelled needle's row position
%    dorsoVentralExtent - extent of dorso-ventral bending/rotation to be captured
%    currentColumnPosition - the labelled needle's column position
%    lateralExtent - lateral extent of pixels to be considered in maximum
%                       analysis
%    centralSliceNumber - number of central slice
%    j - directional distance of slice from central slice
% Outputs:
%    intensity - 2D image patch 
%
% Other m-files required: 
% none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% June 2020

%------------- CODE --------------

    % Get the image patch if dimensional requirements are fulfilled.
    if (currentRowPosition-dorsoVentralExtent) > 0 && (currentColumnPosition-lateralExtent) > 0 && ...
            currentRowPosition+dorsoVentralExtent <= size(image,1) && currentColumnPosition+lateralExtent <= size(image,2)
        intensity = image(currentRowPosition-dorsoVentralExtent:1:currentRowPosition+dorsoVentralExtent,...
            currentColumnPosition-lateralExtent:1:currentColumnPosition+lateralExtent,centralSliceNumber+j);
    else 
        % Calculate corrections for all four image limits.
        correctionLowerRow = 0;
        correctionLowerColumn = 0;
        correctionUpperRow = 0;
        correctionUpperColumn = 0;
        startRow = 1;
        startColumn = 1;
        stopRow = 2*dorsoVentralExtent+1;
        stopColumn = 2*lateralExtent+1;

        if currentRowPosition-dorsoVentralExtent <= 0
            correctionLowerRow = abs(currentRowPosition-dorsoVentralExtent)+1;
            startRow = 1+correctionLowerRow;
        end
        if currentColumnPosition-lateralExtent <= 0
            correctionLowerColumn = abs(currentColumnPosition-lateralExtent)+1;
            startColumn = 1+correctionLowerColumn;
        end
        if currentRowPosition+dorsoVentralExtent > size(image,1)
            correctionUpperRow = size(image,1)-(currentRowPosition+dorsoVentralExtent);
            stopRow = 2*dorsoVentralExtent+1+correctionUpperRow;
        end
        if currentColumnPosition+lateralExtent > size(image,2)
            correctionUpperColumn = size(image,2)-(currentColumnPosition+lateralExtent);
            stopColumn = 2*lateralExtent+1+correctionUpperColumn;
        end
        % Get image patch for these special cases.
        intensity = zeros(2*dorsoVentralExtent+1,2*lateralExtent+1);
        intensity(startRow:1:stopRow,startColumn:1:stopColumn) =...
            image(currentRowPosition-dorsoVentralExtent+correctionLowerRow:1:...
            currentRowPosition+dorsoVentralExtent+correctionUpperRow,currentColumnPosition...
            -lateralExtent+correctionLowerColumn:1:currentColumnPosition+lateralExtent+...
            correctionUpperColumn,centralSliceNumber+j);
    end
end

function [currentColumnPosition] = getCurrentColumnPosition(oldColumnPosition,column,xRotationalRange)
% GETCURRENTCOLUMNPOSITION - Determine the labelled needle's column
% position on the respective slice.
% The function finds the needle's column position for segmentation on the 
% respective slice based on the column position of the same needle in the
% previous slice located closer to the image stack centre and the column
% containing the current slice's maximum.
%
% Inputs:
%    oldColumnPosition - the labelled needle's column position in the
%                        adjacent slice closer to image centre
%    column - column position of the intensity maximim on the current slice
%    xRotationalRange - x-range in pixels to be considered in maximum
%                       analysis
% Outputs:
%    currentColumnPosition - labelled needle's column position on the
%                            current slice
%
% Other m-files required: 
% none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% June 2020

%------------- CODE --------------

    % Create an array containing all summands that can be added to the
    % needle position to calculate the new column position from the one on 
    % the previous slice.
    numbers = 0:1:floor(xRotationalRange/2);
    negativeNumbers = -flip(numbers(2:end));
    summands = cat(2,negativeNumbers,numbers);

    % Update the column position using the column position of the current
    % slice's maximum.
    currentColumnPosition = oldColumnPosition+summands(column);
end

function [voxels] = coordinates2voxels(info,coordinates)
% COORDINATES2VOXELS - Convert 3D point coordinates from the DICOM patient-
% based coordinate system into voxels of the 3D image matrix as defined by
% DICOM parameters ImageOrientationPatient, ImagePositionPatient and 
% PixelSpacing. 
% The mapping is based on explanations by the Neuroimaging for Python 
% community (NIPY). These can be found under:
% https://nipy.org/nibabel/dicom/dicom_orientation.html.
%
% Inputs:
%    info - dicominfo of the RTSTRUCT-assigned DICOM images
%    coordinates - 3D DICOM point coordinates to be converted into voxels
%                  as a 3 x n array
%
% Outputs:
%    voxels - array with 3 x n matlab voxel indices
%
% Other m-files required: none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% May 2020

%------------- CODE --------------

    % Get x,y and z point coordinates in the patient-based coordinate system.
    px = coordinates(1,:);
    py = coordinates(2,:);
    pz = coordinates(3,:);
    
    % Array for point voxels to be calculated.
    voxels = zeros(size(coordinates,1,2));
    
    % For all points.
    for i = 1:length(px)
        % Get single point.
        pointCoordinates = ([px(i); py(i); pz(i); 1]);
        % Get total number of slices of the 3D image stack.
        numberOfSlices = size(info,2);
        % Conversion matrix as described in the reference.
        conversionMatrix = ([info{1}.ImageOrientationPatient(4)*info{1}.PixelSpacing(1)...
            info{1}.ImageOrientationPatient(1)*info{1}.PixelSpacing(2) ((info{1}.ImagePositionPatient(1)...
            -info{numberOfSlices}.ImagePositionPatient(1))/(1-numberOfSlices)) info{1}.ImagePositionPatient(1);
                            info{1}.ImageOrientationPatient(5)*info{1}.PixelSpacing(1)...
                            info{1}.ImageOrientationPatient(2)*info{1}.PixelSpacing(2)...
                            ((info{1}.ImagePositionPatient(2)-info{numberOfSlices}.ImagePositionPatient(2))...
                            /(1-numberOfSlices)) info{1}.ImagePositionPatient(2);
                            info{1}.ImageOrientationPatient(6)*info{1}.PixelSpacing(1)...
                            info{1}.ImageOrientationPatient(3)*info{1}.PixelSpacing(2)...
                            ((info{1}.ImagePositionPatient(3)-info{numberOfSlices}.ImagePositionPatient(3))...
                            /(1-numberOfSlices)) info{1}.ImagePositionPatient(3);
                            0 0 0 1]);
                        
        % Solve and round to voxel indices.
        voxelOutput = int16(round(conversionMatrix\pointCoordinates));
       
        % Adapt to matlab indexing as compared to Python indexing.
        voxelOutput(1) = voxelOutput(1)+1;
        voxelOutput(2) = voxelOutput(2)+1;
        voxelOutput(3) = voxelOutput(3)+1;
        
        % Store in voxel array.
        voxels(:,i) = voxelOutput(1:3);
    end
end

function [coordinates] = voxels2coordinates(info,voxels)
% VOXELS2COORDINATES - Convert 3D voxel indices into the DICOM patient-
% based coordinate system.
% The mapping is based on explanations by the Neuroimaging for Python 
% community (NIPY). These can be found under:
% https://nipy.org/nibabel/dicom/dicom_orientation.html.
%
% Inputs:
%    info - dicominfo of the RTSTRUCT-assigned DICOM images
%    voxels - 3D voxel indices as a 3 x n array to be converted into DICOM 
%    point coordinates
%
% Outputs:
%    coordinates - array with 3 x n DICOM coordinates
%
% Other m-files required: none
%
% Author: Juliane Peter
% email: juliane-peter@web.de
% June 2020

%------------- CODE --------------

    % Get x,y and z voxel indices.
    px = voxels(1,:);
    py = voxels(2,:);
    pz = voxels(3,:);
    
    % Array for point voxels to be calculated.
    coordinates = zeros(size(voxels,1,2));
    
    % For all points.
    for i = 1:length(px)
        % Get single point in matlab indices.
        voxelIndices = ([px(i); py(i); pz(i); 2]);
        
        % Adapt to Python indexing.
        voxelIndices = voxelIndices-1;
        
        % Get total number of slices of 3D image stack.
        numberOfSlices = size(info,2);
        % Conversion matrix as described in reference.
        conversionMatrix = ([info{1}.ImageOrientationPatient(4)*info{1}.PixelSpacing(1)...
            info{1}.ImageOrientationPatient(1)*info{1}.PixelSpacing(2) ((info{1}.ImagePositionPatient(1)...
            -info{numberOfSlices}.ImagePositionPatient(1))/(1-numberOfSlices)) info{1}.ImagePositionPatient(1);
                            info{1}.ImageOrientationPatient(5)*info{1}.PixelSpacing(1)...
                            info{1}.ImageOrientationPatient(2)*info{1}.PixelSpacing(2)...
                            ((info{1}.ImagePositionPatient(2)-info{numberOfSlices}.ImagePositionPatient(2))/...
                            (1-numberOfSlices)) info{1}.ImagePositionPatient(2);
                            info{1}.ImageOrientationPatient(6)*info{1}.PixelSpacing(1)...
                            info{1}.ImageOrientationPatient(3)*info{1}.PixelSpacing(2)...
                            ((info{1}.ImagePositionPatient(3)-info{numberOfSlices}.ImagePositionPatient(3))/...
                            (1-numberOfSlices)) info{1}.ImagePositionPatient(3);
                            0 0 0 1]);
        % Solve.
        coordinateOutput = conversionMatrix*voxelIndices;
        
        % Store in coordinate array.
        coordinates(:,i) = coordinateOutput(1:3);
    end
end

function  imshow3Dfull( Img, disprange )
%IMSHOW3DFULL displays 3D grayscale or RGB images from three perpendicular
%views (i.e. axial, sagittal, and coronal) in slice by slice fashion with
%mouse based slice browsing and window and level adjustment control.
%
% Usage:
% imshow3Dfull ( Image )
% imshow3Dfull ( Image , [] )
% imshow3Dfull ( Image , [LOW HIGH] )
%   
%    Image:      3D image MxNxKxC (K slices of MxN images) C is either 1
%                (for grayscale images) or 3 (for RGB images)
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the widest available range)
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment). Window and level adjustment control works only for
% grayscale images.
%
% Use 'A', 'S', and 'C' buttons to switch between axial, sagittal and
% coronal views, respectivelly.
% 
% "Auto W/L" button adjust the window and level automatically 
%
% While "Fine Tune" check box is checked the window/level adjustment gets
% 16 times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse based window and level adjustment is set
% based on the user defined display intensity range; the wider the range
% the more sensitivity to mouse drag.
% 
% 
%   Example
%   --------
%       % Display an image (MRI example)
%       load mri 
%       Image = squeeze(D); 
%       figure, 
%       imshow3Dfull(Image) 
%
%       % Display the image, adjust the display range
%       figure,
%       imshow3Dfull(Image,[20 100]);
%
%   See also IMSHOW.

%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% - Revision: 2.0.0   Date: 2014/08/05
% - Revision: 2.5.0   Date: 2016/09/22
% - Revision: 2.5.1   Date: 2018/10/29
% 

sno = size(Img);  % image size
sno_a = sno(3);  % number of axial slices
S_a = round(sno_a/2);
sno_s = sno(2);  % number of sagittal slices
S_s = round(sno_s/2);
sno_c = sno(1);  % number of coronal slices
S_c = round(sno_c/2);
S = S_a;
sno = sno_a;

global InitialCoord;

MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end 

ImgAx = Img;
if verLessThan('matlab', '8')
    ImgSg = flipdim(permute(Img, [3 1 2 4]),1);   % Sagittal view image
    ImgCr = flipdim(permute(Img, [3 2 1 4]),1);   % Coronal view image
else
    ImgSg = flip(permute(Img, [3 1 2 4]),1);   % Sagittal view image
    ImgCr = flip(permute(Img, [3 2 1 4]),1);   % Coronal view image
end

View = 'A';

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
VwFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
elseif numel(disprange) == 0
    [Rmin Rmax] = WL2R(Win, LevV);
else
    LevV = (double(disprange(2)) + double(disprange(1))) / 2;
    Win = double(disprange(2)) - double(disprange(1));
    WLAdjCoe = (Win + 1)/1024;
    [Rmin Rmax] = WL2R(Win, LevV);
end

clf
hdl_im = axes('position',[0,0.2,1,0.8]);
imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])

FigPos = get(gcf,'Position');
S_Pos = [50 50 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [50 70 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [20 20 60 20];
Wval_Pos = [75 20 60 20];
Ltxt_Pos = [140 20 45 20];
Lval_Pos = [180 20 60 20];
BtnStPnt = uint16(FigPos(3)-210)+1;
if BtnStPnt < 360
    BtnStPnt = 360;
end
Btn_Pos = [BtnStPnt 20 80 20];
ChBx_Pos = [BtnStPnt+90 20 100 20];
Vwtxt_Pos = [255 20 35 20];
VAxBtn_Pos = [290 20 15 20];
VSgBtn_Pos = [310 20 15 20];
VCrBtn_Pos = [330 20 15 20];

if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end    
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);
Vwtxthand = uicontrol('Style', 'text','Position', Vwtxt_Pos,'String','View: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
VAxBtnhand = uicontrol('Style', 'pushbutton','Position', VAxBtn_Pos,'String','A', 'FontSize', BtnSz, 'Callback' , @AxialView);
VSgBtnhand = uicontrol('Style', 'pushbutton','Position', VSgBtn_Pos,'String','S', 'FontSize', BtnSz, 'Callback' , @SagittalView);
VCrBtnhand = uicontrol('Style', 'pushbutton','Position', VCrBtn_Pos,'String','C', 'FontSize', BtnSz, 'Callback' , @CoronalView);

set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-210)+1;
        if BtnStPnt < 360
            BtnStPnt = 360;
        end
        Btn_Pos = [BtnStPnt 20 80 20];
        ChBx_Pos = [BtnStPnt+90 20 100 20];
        if sno > 1
            set(shand,'Position', S_Pos);
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
        set(Vwtxthand,'Position', Vwtxt_Pos);
        set(VAxBtnhand,'Position', VAxBtn_Pos);
        set(VSgBtnhand,'Position', VSgBtn_Pos);
        set(VCrBtnhand,'Position', VCrBtn_Pos);
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

% -=< Axial view callback function >=-
    function AxialView(object,eventdata)
        if View == 'S'
            S_s = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'A';
        
        Img = ImgAx;
        S = S_a;
        sno = sno_a;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
    end

% -=< Sagittal view callback function >=-
    function SagittalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'S';

        Img = ImgSg;
        S = S_s;
        sno = sno_s;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])
        daspect([17.7054,1,1]);

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);

    end

% -=< Coronal view callback function >=-
    function CoronalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'S'
            S_s = S;
        end            
        View = 'C';
        
        Img = ImgCr;
        S = S_c;
        sno = sno_c;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax])
        daspect([17.7054,1,1]);

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
    end

end
% -=< Maysam Shahedi (mshahedi@gmail.com), September 22, 2016>=-

function  imshow3Dfullcolour( Img, disprange, cmap )
%IMSHOW3DFULL displays 3D grayscale or RGB images from three perpendicular
%views (i.e. axial, sagittal, and coronal) in slice by slice fashion with
%mouse based slice browsing and window and level adjustment control.
%
% Usage:
% imshow3Dfull ( Image )
% imshow3Dfull ( Image , [] )
% imshow3Dfull ( Image , [LOW HIGH] )
%   
%    Image:      3D image MxNxKxC (K slices of MxN images) C is either 1
%                (for grayscale images) or 3 (for RGB images)
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the widest available range)
%    cmap:       Colourmap
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment). Window and level adjustment control works only for
% grayscale images.
%
% Use 'A', 'S', and 'C' buttons to switch between axial, sagittal and
% coronal views, respectivelly.
% 
% "Auto W/L" button adjust the window and level automatically 
%
% While "Fine Tune" check box is checked the window/level adjustment gets
% 16 times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse based window and level adjustment is set
% based on the user defined display intensity range; the wider the range
% the more sensitivity to mouse drag.
% 
% 
%   Example
%   --------
%       % Display an image (MRI example)
%       load mri 
%       Image = squeeze(D); 
%       figure, 
%       imshow3Dfull(Image) 
%
%       % Display the image, adjust the display range
%       figure,
%       imshow3Dfull(Image,[20 100]);
%
%   See also IMSHOW.

%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% - Revision: 2.0.0   Date: 2014/08/05
% - Revision: 2.5.0   Date: 2016/09/22
% - Revision: 2.5.1   Date: 2018/10/29
% - Modified by Juliane Peter, April 2020

sno = size(Img);  % image size
sno_a = sno(3);  % number of axial slices
S_a = round(sno_a/2);
sno_s = sno(2);  % number of sagittal slices
S_s = round(sno_s/2);
sno_c = sno(1);  % number of coronal slices
S_c = round(sno_c/2);
S = S_a;
sno = sno_a;

global InitialCoord;

MinV = 0;
MaxV = max(Img(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Img,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint32')
    MaxV = uint32(Inf);
    MinV = uint32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'uint64')
    MaxV = uint64(Inf);
    MinV = uint64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int8')
    MaxV = int8(Inf);
    MinV = int8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int16')
    MaxV = int16(Inf);
    MinV = int16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int32')
    MaxV = int32(Inf);
    MinV = int32(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'int64')
    MaxV = int64(Inf);
    MinV = int64(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Img,'logical')
    MaxV = 0;
    MinV = 1;
    LevV =0.5;
    Win = 1;
    WLAdjCoe = 0.1;
end 

ImgAx = Img;
if verLessThan('matlab', '8')
    ImgSg = flipdim(permute(Img, [3 1 2 4]),1);   % Sagittal view image
    ImgCr = flipdim(permute(Img, [3 2 1 4]),1);   % Coronal view image
else
    ImgSg = flip(permute(Img, [3 1 2 4]),1);   % Sagittal view image
    ImgCr = flip(permute(Img, [3 2 1 4]),1);   % Coronal view image
end

View = 'A';

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
VwFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

if (nargin < 2)
    [Rmin Rmax] = WL2R(Win, LevV);
elseif numel(disprange) == 0
    [Rmin Rmax] = WL2R(Win, LevV);
else
    LevV = (double(disprange(2)) + double(disprange(1))) / 2;
    Win = double(disprange(2)) - double(disprange(1));
    WLAdjCoe = (Win + 1)/1024;
    [Rmin Rmax] = WL2R(Win, LevV);
end

clf
hdl_im = axes('position',[0,0.2,1,0.8]);
imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax],'Colormap',cmap)

FigPos = get(gcf,'Position');
S_Pos = [50 50 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [50 70 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [20 20 60 20];
Wval_Pos = [75 20 60 20];
Ltxt_Pos = [140 20 45 20];
Lval_Pos = [180 20 60 20];
BtnStPnt = uint16(FigPos(3)-210)+1;
if BtnStPnt < 360
    BtnStPnt = 360;
end
Btn_Pos = [BtnStPnt 20 80 20];
ChBx_Pos = [BtnStPnt+90 20 100 20];
Vwtxt_Pos = [255 20 35 20];
VAxBtn_Pos = [290 20 15 20];
VSgBtn_Pos = [310 20 15 20];
VCrBtn_Pos = [330 20 15 20];

if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end    
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);
Vwtxthand = uicontrol('Style', 'text','Position', Vwtxt_Pos,'String','View: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
VAxBtnhand = uicontrol('Style', 'pushbutton','Position', VAxBtn_Pos,'String','A', 'FontSize', BtnSz, 'Callback' , @AxialView);
VSgBtnhand = uicontrol('Style', 'pushbutton','Position', VSgBtn_Pos,'String','S', 'FontSize', BtnSz, 'Callback' , @SagittalView);
VCrBtnhand = uicontrol('Style', 'pushbutton','Position', VCrBtn_Pos,'String','C', 'FontSize', BtnSz, 'Callback' , @CoronalView);

set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-210)+1;
        if BtnStPnt < 360
            BtnStPnt = 360;
        end
        Btn_Pos = [BtnStPnt 20 80 20];
        ChBx_Pos = [BtnStPnt+90 20 100 20];
        if sno > 1
            set(shand,'Position', S_Pos);
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
        set(Vwtxthand,'Position', Vwtxt_Pos);
        set(VAxBtnhand,'Position', VAxBtn_Pos);
        set(VSgBtnhand,'Position', VSgBtn_Pos);
        set(VCrBtnhand,'Position', VCrBtn_Pos);
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Img)
        S = round(get(hObj,'Value'));
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata)
        UPDN = eventdata.VerticalScrollCount;
        S = S - UPDN;
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
    end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        set(gcf, 'WindowButtonMotionFcn', '')
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Img(:))-min(Img(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Img(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

% -=< Axial view callback function >=-
    function AxialView(object,eventdata)
        if View == 'S'
            S_s = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'A';
        
        Img = ImgAx;
        S = S_a;
        sno = sno_a;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax],'Colormap',cmap)

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
    end

% -=< Sagittal view callback function >=-
    function SagittalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'S';

        Img = ImgSg;
        S = S_s;
        sno = sno_s;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax],'Colormap',cmap)
        daspect([17.7054,1,1]);

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);

    end

% -=< Coronal view callback function >=-
    function CoronalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'S'
            S_s = S;
        end            
        View = 'C';
        
        Img = ImgCr;
        S = S_c;
        sno = sno_c;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        imshow(squeeze(Img(:,:,S,:)), [Rmin Rmax],'Colormap',cmap)
        daspect([17.7054,1,1]);

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Img});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        set(get(gca,'children'),'cdata',squeeze(Img(:,:,S,:)))
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
    end

end
% -=< Maysam Shahedi (mshahedi@gmail.com), September 22, 2016>=-
