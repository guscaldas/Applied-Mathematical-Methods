%% Método de Newton-Raphson
% Por Gustavo Caldas (gustavo.caldas@eq.ufrj.br) e Oscar Chamberlain (ochamberlain@eq.ufrj.br)
% 
function [X,path,cont,Spath] = newton_raphson(X0,A,B)
%Sendo X vetor coluna
    X = X0; %chute inicial
    path = X0; %histórico
    epsilon = 100; %declarando a tolerância de forma a rodar o script
    Spath = zeros(size(X)); %o histórico do vetor busca
    cont= 0; % contador de iterações
    while (epsilon >= 0.01 && cont<10000) %condição de tolerância e condição de n máximo de iterações
        Xant = X; % X^(k) = X^(k-1)
        % Cálculo de $F^(n)$
        %F = (1/4)*(X.'*A*X)^2 +(0.5)*(X.'*A*X)+B.'*X + C;
        % Cálculo do gradiente
        G = (X.'*A*X)*A*X + A*X+B;
        % Cálculo da Hessiana
        H = 2*(A*X)*((A*X).')+(X.'*A*X)*A + A;
        % Cálculo do vetor de busca
        S = -H\G; % É mais rápido que inv(H)*G
        %Vetor S unitário para traçar graficamente e calcular o freio t
        Snorm = S/norm(S);
        %Limitando o passo t
        t = 1;
        while(norm(S)*t>100)
          t = 100/(norm(S)+0.1);
        end
        %Calculando o novo X
        X = X + t*S;
        %Tolerância
        epsilon = norm(X-Xant);
        %Tensor do histórico de X
        path = cat(2,X,path);
        %Histórico de S
        Spath = cat(2,Snorm,Spath);
        %Contador
        cont=cont+1;
    end
end