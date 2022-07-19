function [bestModel, bestModelFit, bestModelOutput] = bestModelFinder(models,data)
    bestModelFit = 0;
    bestModelOutput=[];
    for i = 1:length(models)
        %for each estimated model
        contFit = 0;
        modelOutput=[];
        for j = 1:length(data)        
            [y,fit] = compare(data{j}, models{i});
            contFit = contFit + fit;
            y1 = cell2mat(get(y).OutputData);
            modelOutput=[modelOutput; y1];
        end
        contFit = contFit/length(data);
        if contFit > bestModelFit
            bestModelFit = contFit;
            bestModelOutput=modelOutput;
            bestModel=models{i};
        end
%         fprintf('Model %d mean fit: %.2f\n', i, contFit);
    end
%     fprintf('\nBest model fit %.2f \n', bestModelFit);
%     bestModel
end

