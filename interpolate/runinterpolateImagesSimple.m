clearvars; close all; clc;
tic; 
inputFolder = "your picture.tif";
scaleFactor = 5; 
interpolateImagesSimple(inputFolder, scaleFactor);
