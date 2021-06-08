%    This file is part of cellMat.
%     Copyright (C) 2016 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
%
%     This library is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published
%     by the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This software is provided "as is",
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details
%     <http://www.gnu.org/licenses/>.


function [replaceStrIn,replaceStrOut,deleteLineKeyword] = parseElastixFindShiftParam(options)
% This function parse the parameters needed to run elastix



% Initialize the keywords used to delete useless lines in the parameter file
deleteKeywordCount = 1; deleteLineKeyword{deleteKeywordCount} = [];
% Initialize the strings that have to be replaced in the parameter file
strCount = 1; replaceStrIn{strCount} = []; replaceStrOut{strCount} = [];
% Parse required parameters
replaceStrIn{strCount} = 'pxType'; replaceStrOut{strCount} = 'float'; strCount = strCount+1;
replaceStrIn{strCount} = 'imDim'; replaceStrOut{strCount} = num2str(options.imDim); strCount = strCount+1;
replaceStrIn{strCount} = 'numResolutionVal'; replaceStrOut{strCount} = num2str(options.multiscale.num); strCount = strCount+1;
replaceStrIn{strCount} = 'imPyramidScheduleVal'; replaceStrOut{strCount} = numMtx2Str(options.multiscale.pyramidSchedule); strCount = strCount+1; % [x,y,z] for each scale
switch options.metric.type
    case 'normXCorr'
        replaceStrIn{strCount} = 'metricType'; replaceStrOut{strCount} = 'AdvancedNormalizedCorrelation'; strCount = strCount+1;
        deleteLineKeyword{deleteKeywordCount} = 'NumberOfHistogramBins'; deleteKeywordCount = deleteKeywordCount+1;
    case 'mutualInfo'
        replaceStrIn{strCount} = 'metricType'; replaceStrOut{strCount} = 'AdvancedMattesMutualInformation'; strCount = strCount+1;
        replaceStrIn{strCount} = 'numHistBinVal'; replaceStrOut{strCount} = numMtx2Str(options.metric.mutualInfo.numHistBins); strCount = strCount+1; % it could be a different value for each scale (typically 16 or 32)
    case 'meanSquares'
        replaceStrIn{strCount} = 'metricType'; replaceStrOut{strCount} = 'AdvancedMeanSquares'; strCount = strCount+1;
        deleteLineKeyword{deleteKeywordCount} = 'NumberOfHistogramBins'; deleteKeywordCount = deleteKeywordCount+1;
end
switch options.optim.method
    case 'qnLBFGS'
        replaceStrIn{strCount} = 'optimType'; replaceStrOut{strCount} = 'QuasiNewtonLBFGS'; strCount = strCount+1;
        deleteLineKeyword{deleteKeywordCount} = 'AutomaticParameterEstimation'; deleteKeywordCount = deleteKeywordCount+1;
        deleteLineKeyword{deleteKeywordCount} = 'UseAdaptiveStepSizes'; deleteKeywordCount = deleteKeywordCount+1;
        deleteLineKeyword{deleteKeywordCount} = 'ASGDParameterEstimationMethod'; deleteKeywordCount = deleteKeywordCount+1;
        replaceStrIn{strCount} = 'newSamplesFlag';replaceStrOut{strCount} = 'false'; strCount = strCount+1;
    case 'adapStochGradDesc'
        replaceStrIn{strCount} = 'optimType'; replaceStrOut{strCount} = 'AdaptiveStochasticGradientDescent'; strCount = strCount+1;
        if not(options.optim.asgd.speed)
            deleteLineKeyword{deleteKeywordCount} = 'ASGDParameterEstimationMethod'; deleteKeywordCount = deleteKeywordCount+1;
        end
        replaceStrIn{strCount} = 'newSamplesFlag';replaceStrOut{strCount} = 'true'; strCount = strCount+1;
end

%Replace the transform type
replaceStrIn{strCount} = 'txType'; replaceStrOut{strCount} = options.txType; strCount = strCount+1;


replaceStrIn{strCount} = 'iterNum'; replaceStrOut{strCount} = numMtx2Str(options.optim.iterNum); strCount = strCount+1; %  for each scale
if isempty(options.optim.evalSampleNum)
    replaceStrIn{strCount} = 'imageSamplerVal'; replaceStrOut{strCount} = 'Full'; strCount = strCount+1;
    deleteLineKeyword{deleteKeywordCount} = 'NumberOfSpatialSamples'; deleteKeywordCount = deleteKeywordCount+1;
    deleteLineKeyword{deleteKeywordCount} = 'NewSamplesEveryIteration'; deleteKeywordCount = deleteKeywordCount+1;
else
     replaceStrIn{strCount} = 'imageSamplerVal'; replaceStrOut{strCount} = 'RandomCoordinate'; strCount = strCount+1;
    replaceStrIn{strCount} = 'NumberOfSpatialSamples'; replaceStrOut{strCount} = numMtx2Str(options.optim.evalSampleNum); strCount = strCount+1; % for each scale
end
replaceStrIn{strCount} = 'outputImType';  replaceStrOut{strCount} = options.outImType; strCount = strCount+1;

