function exportImageStatsToExcel(imageStatsArray, outputPath)
    summaryTable = table();
    
    thresholds = imageStatsArray{1}.allThresholds;
    numThresholds = length(thresholds);
    
    densityColNames = arrayfun(@(t) sprintf('Density_%.3f', t), thresholds, 'UniformOutput', false);
    
    for imgIdx = 1:length(imageStatsArray)
        stats = imageStatsArray{imgIdx};
        
        row = table();
        row.ImageIndex = imgIdx;
        row.ImageName = {stats.imageName};
        row.NumValidVesicles = stats.numValidVesicles;
        row.D9595 = stats.d9595;
        
        if isfield(stats, 'density1') && isfield(stats.density1, 'mean')
            densityValues = stats.density1.mean;
            for t = 1:numThresholds
                row.(densityColNames{t}) = densityValues(t);
            end
        else
            for t = 1:numThresholds
                row.(densityColNames{t}) = NaN;
            end
        end
        
        summaryTable = [summaryTable; row];
        
        [matDir, ~] = fileparts(outputPath);
        
        if ~exist(matDir, 'dir')
            mkdir(matDir);
        end
        
        vesicleDetails = stats.vesicleDetails;
        matPath = fullfile(matDir, sprintf('%s_vesicles.mat', stats.imageName));
        save(matPath, 'vesicleDetails');
    end
    
    [excelDir, ~] = fileparts(outputPath);
    if ~exist(excelDir, 'dir')
        mkdir(excelDir);
    end
    
    writetable(summaryTable, outputPath, 'Sheet', 'ImageLevelStats');
    
    try
        thresholds = thresholds(:);
        thresholdInfo = table(thresholds, 'VariableNames', {'ThresholdSequence'});
        writetable(thresholdInfo, outputPath, 'Sheet', 'ThresholdInfo');
    catch
        writecell({'ThresholdSequence'}, outputPath, 'Sheet', 'ThresholdInfo', 'Range', 'A1');
        writematrix(thresholds, outputPath, 'Sheet', 'ThresholdInfo', 'Range', 'A2');
    end
    
    fprintf('Image-level statistics exported to: %s\n', outputPath);
    fprintf('Vesicle-level information saved as MAT files\n');
end