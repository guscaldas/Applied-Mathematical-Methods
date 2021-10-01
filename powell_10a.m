%% Método de Powell
% Por Gustavo Caldas (gustavo.caldas@eq.ufrj.br) e Oscar Chamberlain (ochamberlain@eq.ufrj.br)
% Ref: (C.S.  Smith,  1962 e M.J.D. Powell, 1964)

function [X,path,cont,Spath] = powell_10a(X0,A,B,K)
%Sendo X vetor coluna
    X = X0; %chute inicial
    epsilon = 100;  %declarando a tolerância de forma o rodar o script
    cont= 0; % contador de iterações
    R = X(7:8);
    path = X0; %histórico
    Spath = zeros(size(X)); %o histórico do vetor busca
        while(epsilon >= 0.00001 && cont<10000) 
          %Calculando o novo X
          Rant = R;
          X = K*R;
          %Cálculo da Hessiana em k
          Hr = K.'*(2*(A*X)*((A*X).')+(X.'*A*X)*A + A)*K;
          % Cálculo da Base conjugada da Hessiana
          P = conjugado_schmidt(Hr);
          %Estimativa inicial de X(teta):
          %Rteta(:,1) = P(:,1);
          Rteta(:,1) = R(:,1);
          %Condição de pto estacionário
          pe =  @(teta) ponto_estacionario_10(P(:,2),Rteta(:,1),A,B,teta,K,X);
          teta = fzero(pe,0);
          %teta
          Rteta(:,2) = Rant(:,1) + P(:,2)*teta; 
          R = Rteta(:,2);
        % Cálculo do vetor de busca
            S = R-Rant; 
        % Cálculo do vetor de busca S para todo X
        Stotal = K*S;
        Snorm = Stotal/norm(Stotal);
        %Tolerância
            epsilon = norm(R-Rant);
        %Tensor do histórico de R
            path = cat(2,X,path); %path não somente para R, mas para todo X
        %Histórico de S
            Spath = cat(2,Snorm,Spath);
        %Contador
            cont=cont+1;
        end
end