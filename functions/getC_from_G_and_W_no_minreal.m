%% Function that given the plant G and the feedback loop W=CG/1+CG, returns C
% W=CG/1+CG
% W+WCG=CG
% CG-WCG=W
% C=W/(G-WG)
function C = getC_from_G_and_W_no_minreal(G, W, architecture, beta_pos, beta_force, beta_force_int, beta_adm, beta_imp, Kp_pos, Kd_pos, posPole, J_adm, D_adm, K_imp)
    s = tf('s');
    
    % implementation of position controller, posPole is a high frequency
    % pole to make the system fisible
    P = (Kd_pos*s + Kp_pos); % / (s+posPole);
    switch architecture
        case ArchitectureEnum.COMP_NONE
            C = zpk(W / (G-W*G));
            
        case ArchitectureEnum.FORCE_PLAIN_P 
            C = zpk((W / ( G*beta_force - W*G*beta_force)));
        case ArchitectureEnum.FORCE_MULTICH8
            C = zpk((W / ( G*beta_force - W*G*beta_force)));
        
        case ArchitectureEnum.POS_V_PLAIN_P           
            L = (P*G / (1 + P*G));
            C = zpk((W*s / (L*beta_pos - W*L*beta_pos)));
        case ArchitectureEnum.POS_V_MULTICH8            
            L = (P*G / (1 + P*G));
            C = zpk((W*s / (L*beta_pos - W*L*beta_pos)));
        
        case ArchitectureEnum.FIX_IMP_PLAIN_P 
            L = ((K_imp*G) / (1 + K_imp*G));
            C = zpk((W*s / (L*beta_imp - W*L*beta_imp)));
        case ArchitectureEnum.FIX_IMP_MULTICH8
            L = ((K_imp*G) / (1 + K_imp*G));
            C = zpk((W*s / (L*beta_imp - W*L*beta_imp)));
        
        case ArchitectureEnum.ADM_PLAIN_P 
            L = (P*G / (1 + P*G));
            Amd = 1 / (s*(J_adm*s + D_adm));
            C = zpk((W / (L*beta_adm*Amd - W*L*beta_adm*Amd)));
        case ArchitectureEnum.ADM_MULTICH8
            L = (P*G / (1 + P*G));
            Amd = 1 / (s*(J_adm*s + D_adm));
            C = zpk((W / (L*beta_adm*Amd - W*L*beta_adm*Amd)));
        
        case ArchitectureEnum.FORCE_INT_PLAIN_P 
            C = zpk((W*s / ( G*beta_force_int - W*G*beta_force_int)));
        case ArchitectureEnum.FORCE_INT_MULTICH8
            C = zpk((W*s / ( G*beta_force_int - W*G*beta_force_int)));
    
    end
end
