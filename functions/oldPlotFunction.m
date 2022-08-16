function oldPlotFunction(plot_C,plot_W,plot_G,allTrimmedTorque,allTrimmedPos,C_bestModelOutput,W_bestModelOutput,G_bestModelOutput)
% Controller C
if plot_C
    for j=1:1:ceil(length(allTrimmedTorque)/10)
        figure(j)
        sgtitle('Torque of best estimated C'); 
        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(length(allTrimmedTorque)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)]; 
            hold on;
            title(titleStr) 
            plot(0:length(allTrimmedTorque{testIdx})-1,allTrimmedTorque{testIdx});
            plot(0:length(C_bestModelOutput{testIdx})-1,C_bestModelOutput{testIdx});
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
    end
    k=j;
else
     k=0;
end
% Whole model W

if plot_W
     for j=1:1:ceil(length(allTrimmedTorque)/10)
        figure(j+k)
        sgtitle('Position best estimated W'); 
        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(length(allTrimmedPos)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)]; 
            hold on;
            title(titleStr) 
            plot(0:length(allTrimmedPos{testIdx})-1,allTrimmedPos{testIdx});
            plot(0:length(W_bestModelOutput{testIdx})-1,W_bestModelOutput{testIdx});
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
     end
    h=j;
else
    h=0;
end
% Motor G

if plot_G
     for j=1:1:ceil(length(allTrimmedTorque)/10)
        figure(j+k+h)
        sgtitle('Position best estimated G'); 
        for i=1:1:10
            testIdx=(j-1)*10+i;
            if(length(allTrimmedPos)<testIdx)
                continue
            end
            subplot(2,5,i)
            titleStr=['Test ',num2str(testIdx)]; 
            hold on;
            title(titleStr) 
            plot(0:length(allTrimmedPos{testIdx})-1,allTrimmedPos{testIdx});
            plot(0:length(G_bestModelOutput{testIdx})-1,G_bestModelOutput{testIdx});
            legend('Location','southoutside')
            legend('testing', 'estimated model output')
        end
     end
end
end