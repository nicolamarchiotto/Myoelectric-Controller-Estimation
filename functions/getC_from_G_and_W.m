%% Function that given the plant G and the feedback loop W=CG/1+CG, returns C
% W=CG/1+CG
% W+WCG=CG
% CG-WCG=W
% C=W/(G-WG)
function C = getC_from_G_and_W(G, W, architecture, beta_pos, beta_force, beta_force_int, beta_adm, beta_imp, Kp_pos, Kd_pos, posPole, J_adm, D_adm, K_imp)
    s = tf('s');
    
    % implementation of position controller, posPole is a high frequency
    % pole to make the system fisible
    P = (Kd_pos*s + Kp_pos); % / (s+posPole);
    switch architecture
        case ArchitectureEnum.COMP_NONE
            C = zpk(minreal(W / minreal(G-W*G, 1e-3), 1e-2));
            
        case ArchitectureEnum.FORCE_PLAIN_P 
            C = zpk(minreal(W / minreal( G*beta_force - W*G*beta_force, 1e-3), 1e-2));
        case ArchitectureEnum.FORCE_MULTICH8
            C = zpk(minreal(W / minreal( G*beta_force - W*G*beta_force, 1e-3), 1e-2));
        
        case ArchitectureEnum.POS_V_PLAIN_P           
            L = minreal(P*G / (1 + P*G));
            C = zpk(minreal(W*s / minreal(L*beta_pos - W*L*beta_pos)));
        case ArchitectureEnum.POS_V_MULTICH8            
            L = minreal(P*G / (1 + P*G), 1e-3);
            C = zpk(minreal(W*s / minreal(L*beta_pos - W*L*beta_pos, 1e-3), 1e-2));
        
        case ArchitectureEnum.FIX_IMP_PLAIN_P 
            L = minreal((K_imp*G) / (1 + K_imp*G), 1e-3);
            C = zpk(minreal(W*s / minreal(L*beta_imp - W*L*beta_imp, 1e-3), 1e-2));
        case ArchitectureEnum.FIX_IMP_MULTICH8
            L = minreal((K_imp*G) / (1 + K_imp*G), 1e-3);
            C = zpk(minreal(W*s / minreal(L*beta_imp - W*L*beta_imp, 1e-3), 1e-2));
        
        case ArchitectureEnum.ADM_PLAIN_P 
            L = minreal(P*G / (1 + P*G), 1e-3);
            Amd = 1 / (s*(J_adm*s + D_adm));
            C = zpk(minreal(W / minreal(L*beta_adm*Amd - W*L*beta_adm*Amd, 1e-3), 1e-2));
        case ArchitectureEnum.ADM_MULTICH8
            L = minreal(P*G / (1 + P*G), 1e-3);
            Amd = 1 / (s*(J_adm*s + D_adm));
            C = zpk(minreal(W / minreal(L*beta_adm*Amd - W*L*beta_adm*Amd, 1e-3), 1e-2));
        
        case ArchitectureEnum.FORCE_INT_PLAIN_P 
            C = zpk(minreal(W*s / minreal( G*beta_force_int - W*G*beta_force_int, 1e-3), 1e-2));
        case ArchitectureEnum.FORCE_INT_MULTICH8
            C = zpk(minreal(W*s / minreal( G*beta_force_int - W*G*beta_force_int, 1e-3), 1e-2));
    
    end
end
