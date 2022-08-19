function G_plotFunction(allTrimmedTorque, allTrimmedPos, G_bestModelOutput, saveImages, path, architecture)

for j=1:1:ceil(size(allTrimmedTorque,1)/10)
    figure(j)
    sgtitle('Position best estimated G');
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
        plot(1:size(G_bestModelOutput,2),G_bestModelOutput(testIdx,:));
        legend('Location','southoutside')
        legend('testing', 'estimated model output')
    end
    if saveImages
        savePath = path + string(architecture) + "_" + j + ".png";
        print(gcf, savePath ,'-dpng','-r300');
    end
end


