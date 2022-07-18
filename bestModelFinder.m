function [bestModelIdx, bestModelFit] = bestModelFinder(models,data)
    bestModelIdx = 1;
    bestModelFit = 0;
    for i = 1:length(models)
        %for each estimated model
        contFit = 0;
        for j = 1:length(data)        
            [y,fit] = compare(data{j}, models{i});
            contFit = contFit + fit;
        end
        contFit = contFit/length(data);
        if contFit > bestModelFit
            bestModelFit = contFit;
            bestModelIdx = i;
        end
        fprintf('Model %d mean fit: %.2f\n', i, contFit);
    end
    fprintf('\nBest model %d with %.2f fit\n', bestModelIdx, bestModelFit);

end

