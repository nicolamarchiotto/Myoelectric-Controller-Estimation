function C = invC(G,H)
    C = zpk(minreal(H/(G-H*G)));
end

