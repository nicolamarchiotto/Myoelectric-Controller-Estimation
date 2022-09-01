function simulink_plot_function(saveImages, path, architecture, betaPath, out)

%
%
% C estimation from position error and command torque
%
%

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
grid on;
plot(out.ref.Time(1:end), out.ref.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_est_grav.Time(1:end), out.C_est_grav.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_est_pos.Time(1:end), out.C_est_pos.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_est_force.Time(1:end), out.C_est_force.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_est_force_int.Time(1:end), out.C_est_force_int.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_est_adm.Time(1:end), out.C_est_adm.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_est_imp.Time(1:end), out.C_est_imp.data(1:end),'LineWidth', 1.2);
xlabel('time [s]');
ylabel('position [rad]');
titleStr = sprintf('%s\n%s', string(architecture), 'C estimated from position error and command torque');
title(titleStr)

legend('reference','gravity control', 'position control', 'force control', 'integral control', 'admittance control', 'impedance control', 'Location', 'southwest');

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
plot(out.ref.Time(1:end), out.ref.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_fromW_est_grav.Time(1:end), out.C_fromW_est_grav.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_fromW_est_pos.Time(1:end), out.C_fromW_est_pos.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_fromW_est_force.Time(1:end), out.C_fromW_est_force.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_fromW_est_force_int.Time(1:end), out.C_fromW_est_force_int.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_fromW_est_adm.Time(1:end), out.C_fromW_est_adm.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_fromW_est_imp.Time(1:end), out.C_fromW_est_imp.data(1:end),'LineWidth', 1.2);
xlabel('time [s]');
ylabel('position [rad]');
titleStr = sprintf('%s\n%s', string(architecture), 'Extracted C from estimated W');
title(titleStr)
legend('reference','gravity control', 'position control', 'force control', 'integral control', 'admittance control', 'impedance control', 'Location', 'southwest');

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
plot(out.ref.Time(1:end), out.ref.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_freq_est_grav.Time(1:end), out.C_freq_est_grav.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_freq_est_pos.Time(1:end), out.C_freq_est_pos.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_freq_est_force.Time(1:end), out.C_freq_est_force.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_freq_est_force_int.Time(1:end), out.C_freq_est_force_int.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_freq_est_adm.Time(1:end), out.C_freq_est_adm.data(1:end),'LineWidth', 1.2);
hold on, plot(out.C_freq_est_imp.Time(1:end), out.C_freq_est_imp.data(1:end),'LineWidth', 1.2);
xlabel('time [s]');
ylabel('position [rad]');
titleStr = sprintf('%s\n%s', string(architecture), 'Extracted C from estimated W, then fitted with PD via frequency analysis');
title(titleStr)
legend('reference','gravity control', 'position control', 'force control', 'integral control', 'admittance control', 'impedance control', 'Location', 'southwest');

if saveImages
    savePath = path + "C_freq_est\\" + betaPath + string(architecture) + ".png";
    print(gcf, savePath ,'-dpng','-r300');
end

%
%
% GRAV_COMP case, estimated W as a second order system and then extract a
% PD controller from it
%
%

if architecture == ArchitectureEnum.COMP_NONE
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on;
    grid on;
    plot(out.ref.Time(1:end), out.ref.data(1:end),'LineWidth', 2);
    hold on, plot(out.C_fromW2_est_grav.Time(1:end), out.C_fromW2_est_grav.data(1:end),'LineWidth', 1.2);
    hold on, plot(out.C_fromW2_est_pos.Time(1:end), out.C_fromW2_est_pos.data(1:end),'LineWidth', 1.2);
    hold on, plot(out.C_fromW2_est_force.Time(1:end), out.C_fromW2_est_force.data(1:end),'LineWidth', 1.2);
    hold on, plot(out.C_fromW2_est_force_int.Time(1:end), out.C_fromW2_est_force_int.data(1:end),'LineWidth', 1.2);
    hold on, plot(out.C_fromW2_est_adm.Time(1:end), out.C_fromW2_est_adm.data(1:end),'LineWidth', 1.2);
    hold on, plot(out.C_fromW2_est_imp.Time(1:end), out.C_fromW2_est_imp.data(1:end),'LineWidth', 1.2);
    xlabel('time [s]');
    ylabel('position [rad]');
    titleStr = sprintf('%s\n%s', string(architecture), 'Extracted C from estimated W as a second order system');
    title(titleStr)
    legend('reference','gravity control', 'position control', 'force control', 'integral control', 'admittance control', 'impedance control', 'Location', 'southwest');

    if saveImages
        savePath = path + "C_fromW_est\\" + betaPath + string(architecture) + "2.png";
        print(gcf, savePath ,'-dpng','-r300');
    end
end
end