%% Function that given the plant G and the feedback loop W=CG/1+CG, returns C
% W=CG/1+CG
% W+WCG=CG
% CG-WCG=W
% C=W/(G-WG)
function C = getC_from_G_and_W(G, W)
    C = zpk(minreal(W/minreal(G-W*G,1e-3),1e-2));
end

