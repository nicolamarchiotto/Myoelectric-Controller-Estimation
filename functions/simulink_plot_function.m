function simulink_plot_function(saveImages, path, architecture, betaPath, out, allTrimmedPosRef, allTrimmedPos)

%
%
% C estimation from position error and command torque
%
%

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
grid on;
plot(out.ref.Time(1:end), out.ref.data(1:end),'LineWidth', 2);

switch architecture
    %GRAV COMP
    case ArchitectureEnum.COMP_NONE
        plot(out.C_est_grav.Time(1:end), out.C_est_grav.data(1:end),'LineWidth', 2);
        
        % FORCE
    case ArchitectureEnum.FORCE_PLAIN_P
        plot(out.C_est_force.Time(1:end), out.C_est_force.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FORCE_MULTICH8
        plot(out.C_est_force.Time(1:end), out.C_est_force.data(1:end),'LineWidth', 2);
        
        %POS_V
    case ArchitectureEnum.POS_V_PLAIN_P
        plot(out.C_est_pos.Time(1:end), out.C_est_pos.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.POS_V_MULTICH8
        plot(out.C_est_pos.Time(1:end), out.C_est_pos.data(1:end),'LineWidth', 2);
        
        % IMP
    case ArchitectureEnum.FIX_IMP_PLAIN_P
        plot(out.C_est_imp.Time(1:end), out.C_est_imp.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FIX_IMP_MULTICH8
        plot(out.C_est_imp.Time(1:end), out.C_est_imp.data(1:end),'LineWidth', 2);
        
        % ADM
    case ArchitectureEnum.ADM_PLAIN_P
        plot(out.C_est_adm.Time(1:end), out.C_est_adm.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.ADM_MULTICH8
        plot(out.C_est_adm.Time(1:end), out.C_est_adm.data(1:end),'LineWidth', 2);
        
        % FORCE INT
    case ArchitectureEnum.FORCE_INT_PLAIN_P
        plot(out.C_est_force_int.Time(1:end), out.C_est_force_int.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FORCE_INT_MULTICH8
        plot(out.C_est_force_int.Time(1:end), out.C_est_force_int.data(1:end),'LineWidth', 2);
end

for i = 1:size(allTrimmedPosRef,1)
    if allTrimmedPosRef(1,1) == allTrimmedPosRef(i,1)
        plot(out.C_est_force_int.Time(1:10:end), allTrimmedPos(i, 1:length(out.C_est_force_int.Time(1:10:end)/10)) - allTrimmedPos(i, 1),'LineWidth', 0.15, 'Color', [0.3, 0.75, 0.25]);
    end
end

xlabel('time [s]');
ylabel('position [rad]');
titleStr = sprintf('%s\n%s', string(architecture), 'C estimated from position error and command torque');
title(titleStr)

legend('reference','response','exps output', 'Location', 'northeast');

if saveImages
    savePath = path + "C_est\\" + betaPath + string(architecture) + ".png";
    print(gcf, savePath ,'-dpng','-r300');
end

%
%
% Extracted C from estimated W from reference position and output position
%
%

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
grid on;
plot(out.ref.Time(1:end), out.ref.data(1:end),'LineWidth', 2);

switch architecture
    %GRAV COMP
    case ArchitectureEnum.COMP_NONE
        plot(out.C_fromW_est_grav.Time(1:end), out.C_fromW_est_grav.data(1:end),'LineWidth', 2);
        
        % FORCE
    case ArchitectureEnum.FORCE_PLAIN_P
        plot(out.C_fromW_est_force.Time(1:end), out.C_fromW_est_force.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FORCE_MULTICH8
        plot(out.C_fromW_est_force.Time(1:end), out.C_fromW_est_force.data(1:end),'LineWidth', 2);
        
        %POS_V
    case ArchitectureEnum.POS_V_PLAIN_P
        plot(out.C_fromW_est_pos.Time(1:end), out.C_fromW_est_pos.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.POS_V_MULTICH8
        plot(out.C_fromW_est_pos.Time(1:end), out.C_fromW_est_pos.data(1:end),'LineWidth', 2);
        
        % IMP
    case ArchitectureEnum.FIX_IMP_PLAIN_P
        plot(out.C_fromW_est_imp.Time(1:end), out.C_fromW_est_imp.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FIX_IMP_MULTICH8
        plot(out.C_fromW_est_imp.Time(1:end), out.C_fromW_est_imp.data(1:end),'LineWidth', 2);
        
        % ADM
    case ArchitectureEnum.ADM_PLAIN_P
        plot(out.C_fromW_est_adm.Time(1:end), out.C_fromW_est_adm.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.ADM_MULTICH8
        plot(out.C_fromW_est_adm.Time(1:end), out.C_fromW_est_adm.data(1:end),'LineWidth', 2);
        
        % FORCE INT
    case ArchitectureEnum.FORCE_INT_PLAIN_P
        plot(out.C_fromW_est_force_int.Time(1:end), out.C_fromW_est_force_int.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FORCE_INT_MULTICH8
        plot(out.C_fromW_est_force_int.Time(1:end), out.C_fromW_est_force_int.data(1:end),'LineWidth', 2);
end

for i = 1:size(allTrimmedPosRef,1)
    if allTrimmedPosRef(1,1) == allTrimmedPosRef(i,1)
        plot(out.C_est_force_int.Time(1:10:end), allTrimmedPos(i, 1:length(out.C_est_force_int.Time(1:10:end)/10)) - allTrimmedPos(i, 1),'LineWidth', 0.15, 'Color', [0.3, 0.75, 0.25]);
    end
end

xlabel('time [s]');
ylabel('position [rad]');
titleStr = sprintf('%s\n%s', string(architecture), 'Extracted C from estimated W');
title(titleStr)
legend('reference','response','exps output', 'Location', 'northeast');

if saveImages
    savePath = path + "C_fromW_est\\" + betaPath + string(architecture) + ".png";
    print(gcf, savePath ,'-dpng','-r300');
end

%
%
% Extracted C from estimated W from reference position and output position,
% analyzed frequency response and estimate a PD controller from it
%
%

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
grid on;
plot(out.ref.Time(1:end), out.ref.data(1:end),'LineWidth', 2);

switch architecture
    %GRAV COMP
    case ArchitectureEnum.COMP_NONE
        plot(out.C_freq_est_grav.Time(1:end), out.C_freq_est_grav.data(1:end),'LineWidth', 2);
        
        % FORCE
    case ArchitectureEnum.FORCE_PLAIN_P
        plot(out.C_freq_est_force.Time(1:end), out.C_freq_est_force.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FORCE_MULTICH8
        plot(out.C_freq_est_force.Time(1:end), out.C_freq_est_force.data(1:end),'LineWidth', 2);
        
        %POS_V
    case ArchitectureEnum.POS_V_PLAIN_P
        plot(out.C_freq_est_pos.Time(1:end), out.C_freq_est_pos.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.POS_V_MULTICH8
        plot(out.C_freq_est_pos.Time(1:end), out.C_freq_est_pos.data(1:end),'LineWidth', 2);
        
        % IMP
    case ArchitectureEnum.FIX_IMP_PLAIN_P
        plot(out.C_freq_est_imp.Time(1:end), out.C_freq_est_imp.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FIX_IMP_MULTICH8
        plot(out.C_freq_est_imp.Time(1:end), out.C_freq_est_imp.data(1:end),'LineWidth', 2);
        
        % ADM
    case ArchitectureEnum.ADM_PLAIN_P
        plot(out.C_freq_est_adm.Time(1:end), out.C_freq_est_adm.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.ADM_MULTICH8
        plot(out.C_freq_est_adm.Time(1:end), out.C_freq_est_adm.data(1:end),'LineWidth', 2);
        
        % FORCE INT
    case ArchitectureEnum.FORCE_INT_PLAIN_P
        plot(out.C_freq_est_force_int.Time(1:end), out.C_freq_est_force_int.data(1:end),'LineWidth', 2);
    case ArchitectureEnum.FORCE_INT_MULTICH8
        plot(out.C_freq_est_force_int.Time(1:end), out.C_freq_est_force_int.data(1:end),'LineWidth', 2);
        
end

for i = 1:size(allTrimmedPosRef,1)
    if allTrimmedPosRef(1,1) == allTrimmedPosRef(i,1)
        plot(out.C_est_force_int.Time(1:10:end), allTrimmedPos(i, 1:length(out.C_est_force_int.Time(1:10:end)/10)) - allTrimmedPos(i, 1),'LineWidth', 0.15, 'Color', [0.3, 0.75, 0.25]);
    end
end

xlabel('time [s]');
ylabel('position [rad]');
titleStr = sprintf('%s\n%s', string(architecture), 'Extracted C from estimated W, then fitted with PD via frequency analysis');
title(titleStr)
legend('reference','response','exps output', 'Location', 'northeast');

if saveImages
    savePath = path + "C_freq_est\\" + betaPath + string(architecture) + ".png";
    print(gcf, savePath ,'-dpng','-r300');
end

end


