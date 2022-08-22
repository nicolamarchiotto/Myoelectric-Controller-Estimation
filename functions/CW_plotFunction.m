function CW_plotFunction(plotC, plotW, allTrimmedTorque, allTrimmedPos, C_bestModelOutput, W_bestModelOutput, saveImages, path, architecture, betaPath)

if plotC
    for j=1:1:ceil(size(allTrimmedTorque,1)/10)
        figure(j)

        sgtitle('Torque of best estimated C');

        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(size(allTrimmedTorque,1)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)];
            hold on;
            title(titleStr)

            plot(1:size(allTrimmedTorque,2),allTrimmedTorque(testIdx,:));    
            plot(1:size(C_bestModelOutput,2),C_bestModelOutput(testIdx,:));
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
        if saveImages
            savePath = path + "C_estimation\\" + betaPath + string(architecture) + "_" + j + ".png";
            print(gcf, savePath ,'-dpng','-r300');
        end
    end
    h=j;
else
    h=0;
end

if plotW
    for j=1:1:ceil(size(allTrimmedTorque,1)/10)
        figure(j+h)

        sgtitle('Position of best estimated W');

        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(size(allTrimmedTorque,1)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)];
            hold on;
            title(titleStr)

            plot(1:size(allTrimmedPos,2),allTrimmedPos(testIdx,:));    
            plot(1:size(W_bestModelOutput,2),W_bestModelOutput(testIdx,:));
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
        if saveImages
            savePath = path + "W_estimation\\" + betaPath + string(architecture) + "_" + j + ".png";
            print(gcf, savePath ,'-dpng','-r300');
        end
    end
end

