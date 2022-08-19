function G_plotFunction(allTrimmedTorque, allTrimmedPos, allTrimmedVel, G_bestModelOutput, saveImages, path, architecture, G_estCase)

for j=1:1:ceil(size(allTrimmedTorque,1)/10)
    figure(j)
  
    if G_estCase == GEstCasesEnum.TORQUE_VEL_ESTIMATION
          sgtitle('Velocity best estimated G');
    else
        sgtitle('Position best estimated G');
    end
    for i=1:1:10
        testIdx=(j-1)*10+i;
        if(size(allTrimmedTorque,1)<testIdx)
            continue
        end
        subplot(2,5,i)
        titleStr=['Test ',num2str(testIdx)];
        hold on;
        title(titleStr)
        
        if G_estCase == GEstCasesEnum.TORQUE_VEL_ESTIMATION
            plot(1:size(allTrimmedVel,2),allTrimmedVel(testIdx,:));
        else
            plot(1:size(allTrimmedPos,2),allTrimmedPos(testIdx,:));    
        end
            
        plot(1:size(G_bestModelOutput,2),G_bestModelOutput(testIdx,:));
        legend('Location','southoutside')
        legend('testing', 'estimated model output')
    end
    if saveImages
        savePath = path + string(architecture) + "_" + j + ".png";
        print(gcf, savePath ,'-dpng','-r300');
    end
end


