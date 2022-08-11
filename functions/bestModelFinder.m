% For each model in models, which is a structure of IDTF objects, the function
% compares the output signal of the iddata contained in the data structure,
% with the ouput obtained by feeding the model with the input of the
% iddata, a cumulative fit score is recorded. The model which gives the
% best fit score is returned
function [bestModel, bestModelFit, bestModelOutput] = bestModelFinder(models,data)
    bestModelFit = -1000;
    bestModelOutput=[];
    for i = 1:length(models)
        %for each estimated model
        contFit = 0;
        modelOutput=[];
        for j = 1:length(data)        
            [y,fit] = compare(data{j}, models{i});
            contFit = contFit + fit;
            y1 = cell2mat(get(y).OutputData)';
            modelOutput(j,:)=y1;
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

