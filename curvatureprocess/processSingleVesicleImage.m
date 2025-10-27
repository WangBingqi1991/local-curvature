function [img_all_edges, img_d95_high, img_fixed_hotspots, imageStats] = processSingleVesicleImageV16(imagePath, params)
    % V15完整修复版囊泡图像处理系统
    % 修复内容：
    % 1. 严格的面积过滤（确保所有囊泡面积在[15000, 270000]nm²内）
    % 2. 边缘检查使用原始边界点（与V13一致）
    % 3. 完整的统计字段生成（包括hotspotCounts）
    % 4. D(95,95)计算与V13完全一致
    % 5. 增强的诊断输出和错误处理
    
    %===== 1.参数设置与默认值 =====
    if nargin < 2
        params = struct();
    end
    
    defaultParams = struct(...
        'interpolationScale', 3, ...    % 亚像素插值倍数
        'originalPixelSize', 25, ...     % 原始像素尺寸(nm)
        'minAreaNm2', 15000, ...         % 最小囊泡面积(nm²)
        'maxAreaNm2', 270000, ...        % 最大囊泡面积(nm²)
        'minPeakDistanceNm', 12, ...     % 最小峰值距离(nm)
        'fixedThreshold', 0.030, ...     % 固定热点阈值(1/nm)
        'curvatureThresholds', 0.010:0.005:0.15, ... % 热点分析阈值
        'quantilesVesicle', [0.05:0.05:0.95, 0.96,0.97,0.98,0.99,1.0], ... % 囊泡分位点
        'quantilesImage', [0.05:0.05:0.95, 0.96,0.97,0.98,0.99,1.0], ...   % 图片分位点
        'densityBins', linspace(0,1,11), ...          % 密度定义一分箱(热点/边界点)
        'areaDensityBins', linspace(0,0.1,11), ...    % 密度定义二分箱(热点/面积)
        'edgeMargin', 15, ...             % 边界安全距离(像素)
        'smoothingFactor', 5, ...         % 边界平滑系数
        'logScaleFactor', 1000, ...       % 对数变换系数
        'debugMode', true ...             % 调试模式
    );
    
    % 合并参数
    paramNames = fieldnames(defaultParams);
    for i = 1:length(paramNames)
        if ~isfield(params, paramNames{i})
            params.(paramNames{i}) = defaultParams.(paramNames{i});
        end
    end
    
    % 固定参数
    edgeThresholds = [0.05, 0.2];
    expansionRange = [1, 3];
    
    % 初始化诊断信息
    if params.debugMode
        fprintf('\n===== 开始处理: %s =====\n', imagePath);
        tic;
    end
    
    %===== 2.图像预处理 =====
    % 读取图像
    try
        img = imread(imagePath);
        if size(img, 3) == 3
            img = rgb2gray(img);
        end
        if params.debugMode
            fprintf('图像读取成功: %dx%d像素\n', size(img, 1), size(img, 2));
        end
    catch ME
        error('图像读取失败: %s', ME.message);
    end
    
    % 亚像素插值
    [h, w] = size(img);
    newH = h * params.interpolationScale;
    newW = w * params.interpolationScale;
    enlargedImg = imresize(img, [newH, newW], 'bicubic');
    interpolatedPixelSize = params.originalPixelSize / params.interpolationScale;
    
    if params.debugMode
        fprintf('亚像素处理: %dx%d -> %dx%d (%.1fnm/像素)\n',...
                h, w, newH, newW, interpolatedPixelSize);
    end
    
    % 图像增强
    imgDouble = double(enlargedImg);
    imgNorm = (imgDouble - min(imgDouble(:))) / (max(imgDouble(:)) - min(imgDouble(:)));
    
    % 对数变换增强对比度
    logImg = log(1 + params.logScaleFactor * imgDouble);
    logImg = mat2gray(logImg);
    
    % 自适应直方图均衡化
    lowLevel = quantile(imgNorm(:), 0.01);
    highLevel = quantile(imgNorm(:), 0.99);
    if highLevel <= lowLevel
        adjustedImg = imadjust(imgNorm);
    else
        adjustedImg = imadjust(imgNorm, [lowLevel; highLevel], [0; 1]);
    end
    
    % 组合增强与降噪
    combinedImg = 0.7 * logImg + 0.3 * adjustedImg;
    if exist('imguidedfilter', 'file') == 2
        denoisedImg = imguidedfilter(combinedImg);
        if params.debugMode, fprintf('使用导向滤波器降噪\n'); end
    else
        denoisedImg = medfilt2(combinedImg, [3, 3]);
        if params.debugMode, fprintf('使用中值滤波器降噪\n'); end
    end
    
    %===== 3.囊泡分割 =====
    % 边缘检测
    edgeMap = edge(denoisedImg, 'canny', edgeThresholds);
    closedEdges = imclose(edgeMap, strel('disk', 2));
    filledMask = imfill(closedEdges, 'holes');
    
    % 尺寸过滤
    minAreaPixels = params.minAreaNm2 / (interpolatedPixelSize^2);
    maxAreaPixels = params.maxAreaNm2 / (interpolatedPixelSize^2);
    vesicleMask = bwareafilt(filledMask, [minAreaPixels, maxAreaPixels]);
    
    if params.debugMode
        stats_pre = regionprops(vesicleMask, 'Area');
        areas_pre = [stats_pre.Area] * interpolatedPixelSize^2;
        fprintf('初步检测囊泡: %d (面积范围: %.0f-%.0f nm²)\n',...
                length(stats_pre), min(areas_pre), max(areas_pre));
    end
    
    % 膨胀校正
    distToBoundary = bwdist(bwperim(vesicleMask));
    expansionMap = distToBoundary;
    expansionMap(~vesicleMask) = 0;
    expansionMask = expansionMap > expansionRange(1) & expansionMap <= expansionRange(2);
    
    % 创建校正图像
    correctedImg = denoisedImg;
    [expRows, expCols] = find(expansionMask);
    
    if ~isempty(expRows)
        [~, nearestIdxMap] = bwdist(bwperim(vesicleMask));
        nearestIdx = nearestIdxMap(expansionMask);
        [nearestRows, nearestCols] = ind2sub(size(vesicleMask), nearestIdx);
        expIndices = sub2ind(size(correctedImg), expRows, expCols);
        nearestIndices = sub2ind(size(denoisedImg), nearestRows, nearestCols);
        correctedImg(expIndices) = denoisedImg(nearestIndices);
        if params.debugMode
            fprintf('膨胀校正: 处理 %d 个像素\n', numel(expRows));
        end
    end
    
    % 重新检测边缘（使用校正后图像）
    edgeMapCorrected = edge(correctedImg, 'canny', edgeThresholds);
    closedEdgesCorrected = imclose(edgeMapCorrected, strel('disk', 2));
    filledMaskCorrected = imfill(closedEdgesCorrected, 'holes');
    vesicleMaskCorrected = bwareafilt(filledMaskCorrected, [minAreaPixels, maxAreaPixels]);
    
    if params.debugMode
        stats_post = regionprops(vesicleMaskCorrected, 'Area');
        areas_post = [stats_post.Area] * interpolatedPixelSize^2;
        fprintf('校正后囊泡: %d (面积范围: %.0f-%.0f nm²)\n',...
                length(stats_post), min(areas_post), max(areas_post));
    end
    
    %===== 4.曲率分析准备 =====
    % 边界追踪（基于最终掩膜）
    boundaries = bwboundaries(vesicleMaskCorrected);
    numVesicles = length(boundaries);
    
    % 面积计算（基于最终掩膜）
    stats = regionprops(vesicleMaskCorrected, 'Area');
    areasPixels = [stats.Area];
    pixelAreaNm2 = interpolatedPixelSize^2;
    areasNm2 = areasPixels * pixelAreaNm2;
    
    % 创建边界掩模
    [imgHeight, imgWidth] = size(vesicleMaskCorrected);
    edgeMask = false(imgHeight, imgWidth);
    edgeMargin = params.edgeMargin;
    edgeMask(1:edgeMargin, :) = true;
    edgeMask(end-edgeMargin+1:end, :) = true;
    edgeMask(:, 1:edgeMargin) = true;
    edgeMask(:, end-edgeMargin+1:end) = true;
    
    % 初始化变量
    validVesicles = true(numVesicles, 1);
    vesicleCurvatures = cell(numVesicles, 1);
    numQuantilesVesicle = length(params.quantilesVesicle);
    vesicleQuantileMatrix = nan(numVesicles, numQuantilesVesicle);
    numThresholds = length(params.curvatureThresholds);
    
    % 基础图像准备
    baseImg = im2uint8(imgNorm);
    if size(baseImg, 3) == 1
        baseImg = cat(3, baseImg, baseImg, baseImg);
    end
    
    % 创建结果图像
    img_all_edges = baseImg;     % 图1:所有边缘点(红色)
    img_d95_high = baseImg;      % 图2:D95高曲率点(绿色)
    img_fixed_hotspots = baseImg;% 图3:固定阈值热点(黄色)
    
    % 存储囊泡中心位置
    vesicleCenters = zeros(numVesicles, 2);
    
    % 初始化统计结构
    imageStats = struct();
    
    % 初始化存储矩阵
    perimeters = nan(numVesicles, 1);
    hotspotCounts = nan(numVesicles, numThresholds);
    density1 = nan(numVesicles, numThresholds);
    density2 = nan(numVesicles, numThresholds);
    
    % 初始化过滤计数器
    edgeFiltered = false(numVesicles, 1);
    areaFiltered = false(numVesicles, 1);
    
    %===== 5.遍历囊泡处理 =====
    for k = 1:numVesicles
        boundary = boundaries{k};
        x = boundary(:, 2);
        y = boundary(:, 1);
        
        % 边缘检查（使用原始边界点，与V13一致）
        isNearEdge = false;
        for i = 1:length(x)
            if edgeMask(y(i), x(i))
                isNearEdge = true;
                edgeFiltered(k) = true;
                break;
            end
        end
        
        % 严格的面积过滤
        currentArea = areasNm2(k);
        isAreaValid = currentArea >= params.minAreaNm2 && currentArea <= params.maxAreaNm2;
        if ~isAreaValid
            areaFiltered(k) = true;
        end
        
        % 双重过滤条件
        if isNearEdge || ~isAreaValid
            validVesicles(k) = false;
            
            % 详细诊断输出
            if params.debugMode
                reason = '';
                if isNearEdge, reason = [reason, '边缘']; end
                if ~isAreaValid, reason = [reason, '面积']; end
                
                areaInfo = '';
                if ~isAreaValid
                    if currentArea < params.minAreaNm2
                        areaInfo = sprintf('(%.0f < %.0f)', currentArea, params.minAreaNm2);
                    else
                        areaInfo = sprintf('(%.0f > %.0f)', currentArea, params.maxAreaNm2);
                    end
                end
                
                fprintf('过滤囊泡 %d: %s %s\n', k, reason, areaInfo);
            end
            continue;
        end
        
        % 边界平滑
        x_smooth = smoothdata(x, 'gaussian', params.smoothingFactor);
        y_smooth = smoothdata(y, 'gaussian', params.smoothingFactor);
        
        % 计算并存储曲率
        curvature = computeActualCurvature(x_smooth, y_smooth, interpolatedPixelSize);
        vesicleCurvatures{k} = curvature;
        
        % 使用V13分位点矩阵方法
        vesicleQuantileMatrix(k, :) = quantile(curvature, params.quantilesVesicle);
        
        % 计算周长(边界点数量)
        nPoints = length(x_smooth);
        perimeters(k) = nPoints;
        
        % 存储中心位置
        vesicleCenters(k, 1) = mean(x_smooth);
        vesicleCenters(k, 2) = mean(y_smooth);
        
        %===== 可视化处理 =====
        % 图1:所有边缘点(红色)
        for i = 1:length(x_smooth)
            px = round(x_smooth(i));
            py = round(y_smooth(i));
            if px >= 1 && px <= size(img_all_edges, 2) && ...
               py >= 1 && py <= size(img_all_edges, 1)
                img_all_edges(py, px, 1) = 255; % R
                img_all_edges(py, px, 2) = 0;   % G
                img_all_edges(py, px, 3) = 0;   % B
            end
        end
        
        % 统计计算(所有阈值)
        for t = 1:numThresholds
            threshold = params.curvatureThresholds(t);
            thresholdstanrd = psfDeconvolutionFactor(threshold);
            
            % 热点数量
            hotspotCount = sum(curvature >= thresholdstanrd);
            hotspotCounts(k, t) = hotspotCount;
            
            % 密度定义一:热点数量/边界点数量
            density1(k, t) = hotspotCount / nPoints;
            
            % 密度定义二:热点数量/囊泡面积
            density2(k, t) = hotspotCount / areasNm2(k);
        end
    end
    
    numValidVesicles = sum(validVesicles);
    validIdx = find(validVesicles);
    imageStats.numValidVesicles = numValidVesicles;
    
    % 面积验证诊断
    validAreas = areasNm2(validVesicles);
    if params.debugMode
        fprintf('\n=== 过滤结果 ===\n');
        fprintf('总检测囊泡: %d\n', numVesicles);
        fprintf('有效囊泡: %d (%.1f%%)\n', numValidVesicles, numValidVesicles/numVesicles*100);
        fprintf('过滤原因: 边缘(%d), 面积(%d)\n', sum(edgeFiltered), sum(areaFiltered));
        fprintf('有效面积范围: %.0f - %.0f nm² (参数: %.0f-%.0f)\n',...
                min(validAreas), max(validAreas),...
                params.minAreaNm2, params.maxAreaNm2);
    end
    
    %===== 6.计算D(95,95)阈值 =====
    if numValidVesicles >= 5
        % 使用V13方法：24分位点矩阵的第19个值（对应0.95分位）
        d95Values = vesicleQuantileMatrix(validVesicles, 19);
        imageStats.d9595 = quantile(d95Values, 0.95);
        
        if params.debugMode
            fprintf('D(95,95)计算: %.4f (基于%d个囊泡)\n', imageStats.d9595, numValidVesicles);
        end
    else
        imageStats.d9595 = NaN;
        if params.debugMode
            fprintf('警告: 有效囊泡不足(%d<5)，无法计算D(95,95)\n', numValidVesicles);
        end
    end
    
    %===== 7.生成热点图像 =====
    for k = 1:numVesicles
        if ~validVesicles(k)
            continue;
        end
        
        boundary = boundaries{k};
        x = boundary(:, 2);
        y = boundary(:, 1);
        
        % 平滑边界
        x_smooth = smoothdata(x, 'gaussian', params.smoothingFactor);
        y_smooth = smoothdata(y, 'gaussian', params.smoothingFactor);
        curvature = vesicleCurvatures{k};
        
        % 图2:D95高曲率点(绿色)
        if ~isnan(imageStats.d9595)
            d95Points = curvature >= imageStats.d9595;
            for i = 1:length(x_smooth)
                if d95Points(i)
                    px = round(x_smooth(i));
                    py = round(y_smooth(i));
                    for dx = -1:1
                        for dy = -1:1
                            col = px + dx;
                            row = py + dy;
                            if col >= 1 && col <= size(img_d95_high, 2) && ...
                               row >= 1 && row <= size(img_d95_high, 1)
                                img_d95_high(row, col, 1) = 0;   % R
                                img_d95_high(row, col, 2) = 255; % G
                                img_d95_high(row, col, 3) = 0;   % B
                            end
                        end
                    end
                end
            end
        end
        
        % 图3:固定阈值热点(黄色)
        fixedPoints = curvature >= params.fixedThreshold;
        for i = 1:length(x_smooth)
            if fixedPoints(i)
                px = round(x_smooth(i));
                py = round(y_smooth(i));
                for dx = -1:1
                    for dy = -1:1
                        col = px + dx;
                        row = py + dy;
                        if col >= 1 && col <= size(img_fixed_hotspots, 2) && ...
                           row >= 1 && row <= size(img_fixed_hotspots, 1)
                            img_fixed_hotspots(row, col, 1) = 255; % R
                            img_fixed_hotspots(row, col, 2) = 255; % G
                            img_fixed_hotspots(row, col, 3) = 0;   % B
                        end
                    end
                end
            end
        end
        
        % 在所有图像上添加囊泡ID标签
        textStr = sprintf('%d', k);
        fontProps = {'FontSize', 12, 'TextColor', 'yellow', ...
                    'BoxColor', 'black', 'BoxOpacity', 0.6, ...
                    'AnchorPoint', 'Center'};
        
        centerX = vesicleCenters(k, 1);
        centerY = vesicleCenters(k, 2);
        
        img_all_edges = insertText(img_all_edges, [centerX, centerY], textStr, fontProps{:});
        img_d95_high = insertText(img_d95_high, [centerX, centerY], textStr, fontProps{:});
        img_fixed_hotspots = insertText(img_fixed_hotspots, [centerX, centerY], textStr, fontProps{:});
    end
    
    %===== 8.图像级统计 =====
    % 确保所有统计字段都存在
    numThresholds = length(params.curvatureThresholds);
    
    % 8.1 热点数量统计
    hotspotCountStats = struct(...
        'mean', nan(1, numThresholds),...
        'std', nan(1, numThresholds),...
        'median', nan(1, numThresholds),...
        'q25', nan(1, numThresholds),...
        'q75', nan(1, numThresholds));
    
    if numValidVesicles > 0
        for t = 1:numThresholds
            validData = hotspotCounts(validVesicles, t);
            hotspotCountStats.mean(t) = mean(validData);
            hotspotCountStats.std(t) = std(validData);
            hotspotCountStats.median(t) = median(validData);
            quantiles = quantile(validData, [0.25, 0.75]);
            hotspotCountStats.q25(t) = quantiles(1);
            hotspotCountStats.q75(t) = quantiles(2);
        end
    end
    imageStats.hotspotCounts = hotspotCountStats;
    
    % 8.2 密度统计
    density1Stats = struct(...
        'mean', nan(1, numThresholds),...
        'std', nan(1, numThresholds),...
        'median', nan(1, numThresholds),...
        'q25', nan(1, numThresholds),...
        'q75', nan(1, numThresholds));
    
    density2Stats = struct(...
        'mean', nan(1, numThresholds),...
        'std', nan(1, numThresholds),...
        'median', nan(1, numThresholds),...
        'q25', nan(1, numThresholds),...
        'q75', nan(1, numThresholds));
    
    if numValidVesicles > 0
        for t = 1:numThresholds
            % 密度定义一统计
            validData = density1(validVesicles, t);
            density1Stats.mean(t) = mean(validData);
            density1Stats.std(t) = std(validData);
            density1Stats.median(t) = median(validData);
            quantiles = quantile(validData, [0.25, 0.75]);
            density1Stats.q25(t) = quantiles(1);
            density1Stats.q75(t) = quantiles(2);
            
            % 密度定义二统计
            validData = density2(validVesicles, t);
            density2Stats.mean(t) = mean(validData);
            density2Stats.std(t) = std(validData);
            density2Stats.median(t) = median(validData);
            quantiles = quantile(validData, [0.25, 0.75]);
            density2Stats.q25(t) = quantiles(1);
            density2Stats.q75(t) = quantiles(2);
        end
    end
    imageStats.density1 = density1Stats;
    imageStats.density2 = density2Stats;
    
    % 8.3 密度分布统计
    numBins1 = length(params.densityBins) - 1;
    numBins2 = length(params.areaDensityBins) - 1;
    
    imageStats.density1Distribution = nan(numBins1, numThresholds);
    imageStats.density2Distribution = nan(numBins2, numThresholds);
    
    if numValidVesicles > 0
        for t = 1:numThresholds
            % 密度定义一分布
            validData = density1(validVesicles, t);
            imageStats.density1Distribution(:, t) = histcounts(...
                validData, params.densityBins, 'Normalization', 'probability')';
            
            % 密度定义二分布
            validData = density2(validVesicles, t);
            imageStats.density2Distribution(:, t) = histcounts(...
                validData, params.areaDensityBins, 'Normalization', 'probability')';
        end
    end
    
    % 8.4 特定阈值统计
    [~, fixedThresholdIdx] = min(abs(params.curvatureThresholds - params.fixedThreshold));
    imageStats.fixedThresholdStats = struct(...
        'hotspotCount', hotspotCountStats.mean(fixedThresholdIdx), ...
        'density1', density1Stats.mean(fixedThresholdIdx), ...
        'density2', density2Stats.mean(fixedThresholdIdx));
    
    if ~isnan(imageStats.d9595)
        [~, d9595ThresholdIdx] = min(abs(params.curvatureThresholds - imageStats.d9595));
        imageStats.d9595Stats = struct(...
            'hotspotCount', hotspotCountStats.mean(d9595ThresholdIdx), ...
            'density1', density1Stats.mean(d9595ThresholdIdx), ...
            'density2', density2Stats.mean(d9595ThresholdIdx));
    else
        imageStats.d9595Stats = struct(...
            'hotspotCount', NaN, ...
            'density1', NaN, ...
            'density2', NaN);
    end
    
    % 8.5 面积验证信息
    imageStats.areaValidation = struct(...
        'minDetected', min(areasNm2), ...
        'maxDetected', max(areasNm2), ...
        'minValid', min(validAreas), ...
        'maxValid', max(validAreas), ...
        'params', struct('min', params.minAreaNm2, 'max', params.maxAreaNm2), ...
        'filtered', struct('edge', sum(edgeFiltered), 'area', sum(areaFiltered)));
    % 存储每个囊泡的详细信息
% 约第520行：添加维度转换
validVesicles = validVesicles(:);
areasNm2 = areasNm2(:);
perimeters = perimeters(:);

% 约第525行：替换为循环创建结构体
vesicleDetails = struct('isValid', cell(numVesicles,1),...
                        'areaNm2', cell(numVesicles,1),...
                        'perimeter', cell(numVesicles,1),...
                        'hotspotCounts', cell(numVesicles,1));

for k = 1:numVesicles
    vesicleDetails(k).isValid = validVesicles(k);
    vesicleDetails(k).areaNm2 = areasNm2(k);
    vesicleDetails(k).perimeter = perimeters(k);
    vesicleDetails(k).hotspotCounts = hotspotCounts(k, :);
end
% 添加到输出结构
imageStats.vesicleDetails = vesicleDetails;
imageStats.allThresholds = params.curvatureThresholds; % 保存阈值序列
    %===== 9.结束处理 =====
    if params.debugMode
        elapsedTime = toc;
        fprintf('\n处理完成! 耗时: %.2f 秒\n', elapsedTime);
        fprintf('有效囊泡数: %d\n', numValidVesicles);
        if ~isnan(imageStats.d9595)
            fprintf('D(95,95)阈值: %.4f nm⁻¹\n', imageStats.d9595);
        end
    end
    
    function corrected = psfDeconvolutionFactor(inputValue, psfParams)
    if nargin < 2
        psfParams = struct();
    end
    
    if ~isfield(psfParams, 'sigma')
        psfParams.sigma = 1.5;
    end
    
    if ~isfield(psfParams, 'size')
        psfParams.size = 5;
    end
    
    psfKernel = fspecial('gaussian', psfParams.size, psfParams.sigma);
    [U, S, V] = svd(psfKernel);
    conditionNumber = cond(psfKernel);
    deconvFactor = (S(1,1) / S(end,end)) / (2 * conditionNumber);
    
    if isscalar(inputValue)
        corrected = inputValue * deconvFactor;
    else
        corrected = zeros(size(inputValue));
        for i = 1:numel(inputValue)
            corrected(i) = inputValue(i) * deconvFactor;
        end
    end
    end
    %===== 曲率计算函数 =====
    function curvature = computeActualCurvature(x, y, pixelSize)
        % 计算曲率
        dx = gradient(x);
        dy = gradient(y);
        ddx = gradient(dx);
        ddy = gradient(dy);
        curvature = (dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^(3/2);
        
        % 处理异常值
        curvature(isinf(curvature)) = 0;
        curvature(isnan(curvature)) = 0;
        
        % 转换为1/nm单位
        curvature = curvature / pixelSize;
        
        % 平滑处理
        curvature = smoothdata(curvature, 'gaussian', 5);
        
        % 过滤异常高曲率值
        maxValidCurvature = 0.5; % 最大合理曲率(1/nm)
        invalidIdx = abs(curvature) > maxValidCurvature;
        curvature(invalidIdx) = 0;
    end
end