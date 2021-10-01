% Script file: newton_raphson_9.m
%% Método de Newton-Raphson no Espaço Reduzido
% Por Gustavo Caldas (gustavo.caldas@eq.ufrj.br) e Oscar Chamberlain (ochamberlain@eq.ufrj.br)
% 
function [X,path,cont,Spath] = newton_raphson_9(X0,A,B,K)
%Sendo X vetor coluna
    X = X0; %chute inicial
    path = X0; %histórico
    epsilon = 100; %declarando a tolerância de forma o rodar o script
    Spath = zeros(size(X)); %o histórico do vetor busca
    cont= 0; % contador de iterações
    while (epsilon >= 0.001 && cont<10000) %condição de tolerância e condição de n máximo de iterações
        Xant = X; % X^(k) = X^(k-1)
        % Cálculo do gradiente
        G = (X.'*A*X)*A*X + A*X+B;
        % Gradiente reduzido    
        Gr = K.'*G;
        % Cálculo da Hessiana
        H = 2*(A*X)*((A*X).')+(X.'*A*X)*A + A;
        Hr = K.'*H*K;
        % Cálculo do vetor de busca
        S = -Hr\Gr; % É mais rápido que inv(H)*G
        % Cálculo do vetor de busca S para todo X
        Stotal = K*S;
        Snorm = Stotal/norm(Stotal);
        %Limitando o passo t
        t = 1;
        while(norm(S)*t>100)
          t = 100/(norm(S)+0.1);
        end
        %Calculando o novo X
        R = X(7:8,1) + t*S(1:2,1); % S só tem duas posições
        X = K*R;
        %Tolerância
        epsilon = norm(X-Xant);
        %Tensor do histórico de X
        path = cat(2,X,path); %path não somente para R, mas para todo X
        %Histórico de S
        Spath = cat(2,Snorm,Spath);
        %Contador
        cont=cont+1;
    end
end