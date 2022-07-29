%% Function that given the plant G and the feedback loop W=CG/1+CG, returns C 
function C = getC_from_G_and_W(G,W)
    C = zpk(minreal(W/(G-W*G)));
end

