function f = ponto_estacionario_10(P,Rteta_ant,A,B,teta,K,X)
     Rteta = Rteta_ant + P*teta;
     X = K*Rteta;
     %Cálculo do gradiente
     G = (X.'*A*X)*A*X + A*X+B;
     Gr = K.'*G;
     f = P.'*Gr;
end