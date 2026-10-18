function BSE_Pore_Analysis
% Generic BSE pore segmentation example; requires Image Processing Toolbox.
% Not the original CNN/DTD research code. Validate masks before research use.
clc;
% ---------------- User settings ----------------
poreMode = 'bright'; % 'bright': metal-filled pores; 'dark': ordinary BSE pores
cropRows = 890;      % Inf keeps the full image; exclude scale bars manually
pixelSizeUm = NaN;  % Set calibrated micrometres/pixel for physical sizes
manualThreshold = []; % [] = Otsu; otherwise a value from 0 to 1
minAreaPixels = 4;  % Removes small components; changes measured porosity
sigma = 1;         % Gaussian smoothing; 0 disables it
% ------------------------------------------------
inputDir = uigetdir(pwd,'Select BSE image folder');
if isequal(inputDir,0), return; end
outputDir = fullfile(inputDir,'Pore_Results');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
files = dir(inputDir);
files = files(~[files.isdir]);
keep = false(size(files));
for k=1:numel(files)
    [~,~,ext] = fileparts(files(k).name);
    keep(k) = any(strcmpi(ext,{'.tif','.tiff','.png','.jpg','.jpeg','.bmp'}));
end
files = files(keep);
if isempty(files), error('No supported image files found.'); end
summary = table();
for k=1:numel(files)
    try
        [raw,map] = imread(fullfile(inputDir,files(k).name));
        if ~isempty(map), raw = ind2rgb(raw,map); end
        if size(raw,3)==3, raw = rgb2gray(raw); end
        I = im2double(raw);
        if any(~isfinite(I(:))) || any(I(:)<0) || any(I(:)>1)
            error('Image must have finite intensities in [0,1] after im2double.');
        end
        I = I(1:min(size(I,1),cropRows),:);
        if max(I(:))==min(I(:)), error('Constant-intensity image.'); end
        J = I;
        if sigma>0, J = imgaussfilt(J,sigma); end
        if isempty(manualThreshold), threshold = graythresh(J);
        else
            validateattributes(manualThreshold,{'numeric'},{'scalar','>=',0,'<=',1});
            threshold = manualThreshold;
        end
        switch lower(poreMode)
            case 'bright', mask = J>threshold;
            case 'dark', mask = J<threshold;
            otherwise, error('poreMode must be bright or dark.');
        end
        mask = bwareaopen(mask,minAreaPixels,8);
        stats = regionprops(mask,'Area','EquivDiameter','Centroid');
        areaPx = reshape([stats.Area],[],1);
        diameterPx = reshape([stats.EquivDiameter],[],1);
        if isempty(stats), centers=zeros(0,2);
        else, centers=vertcat(stats.Centroid); end
        pores = table((1:numel(stats))',areaPx,diameterPx,centers(:,1),centers(:,2), ...
            'VariableNames',{'PoreID','Area_px','EquivalentDiameter_px','CentroidX_px','CentroidY_px'});
        if isfinite(pixelSizeUm) && pixelSizeUm>0
            pores.Area_um2 = areaPx*pixelSizeUm^2;
            pores.EquivalentDiameter_um = diameterPx*pixelSizeUm;
        end
        % Blend red over segmented regions, leaving the rest unchanged.
        overlay = repmat(I,1,1,3);
        red = overlay(:,:,1); green=overlay(:,:,2); blue=overlay(:,:,3);
        red(mask)=0.55+0.45*red(mask);
        green(mask)=0.45*green(mask); blue(mask)=0.45*blue(mask);
        overlay=cat(3,red,green,blue);
        [~,base,ext] = fileparts(files(k).name);
        base=[base '_' ext(2:end)]; % avoid collisions across input formats
        imwrite(mask,fullfile(outputDir,[base '_mask.png']));
        imwrite(overlay,fullfile(outputDir,[base '_overlay.png']));
        writetable(pores,fullfile(outputDir,[base '_pores.csv']));
        porePercent=100*nnz(mask)/numel(mask);
        row=table({files(k).name},threshold,numel(mask),nnz(mask),porePercent,numel(stats), ...
            'VariableNames',{'Image','Threshold','TotalPixels','PorePixels','AreaPorosity_percent','PoreCount'});
        summary=[summary;row]; %#ok<AGROW>
        fprintf('%s: %.3f %%\n',files(k).name,porePercent);
    catch ME
        warning('Failed to process %s: %s',files(k).name,ME.message);
    end
end
if ~isempty(summary)
    writetable(summary,fullfile(outputDir,'Porosity_Summary.csv'));
    fprintf('Completed %d of %d images. Output: %s\n',height(summary),numel(files),outputDir);
else
    warning('No images were successfully processed.');
end
% Area fraction is a 2D descriptive metric; bright phases may include solids.
% Connected component diameters describe segmented 2D regions, not pore throats.
end
