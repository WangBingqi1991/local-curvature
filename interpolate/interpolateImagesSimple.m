function interpolateImagesSimple(folderPath, scaleFactor)
    [parentDir, folderName] = fileparts(folderPath);
    outputFolderName = "Standard" + folderName;
    outputFolder = fullfile(parentDir, outputFolderName);
    
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
        fprintf('Output folder created: %s\n', outputFolder);
    end
    
    imageFiles = dir(fullfile(folderPath, '*.tif'));
    
    for i = 1:numel(imageFiles)
        imgPath = fullfile(folderPath, imageFiles(i).name);
        img = imread(imgPath);
        
        if size(img, 3) == 3
            img = img(:, :, 1);
        end
        
        [h, w] = size(img);
        newH = floor(h * scaleFactor);
        newW = floor(w * scaleFactor);
        
        interpolated = imresize(img, [newH, newW], 'bicubic');
        
        outputPath = fullfile(outputFolder, imageFiles(i).name);
        imwrite(interpolated, outputPath);
    end
    
    fprintf('Interpolation completed! Results saved to: %s\n', outputFolder);
end