% Objetivo: Este programa executa a Construção de Base Conjugada P1 , P2 , ..., Pn por Matriz Simétrica e Definida, 
% pelo método de Schmidt.
% Referência: S. Chapman, "Programação para Engenheiros"

% Definição do tamanho da matriz
function P = conjugado_schmidt(M)

n = length(M); %dimensão de M

P = zeros(n); % declaração da dimensão de P

W = M;

P(:,1) = W(:,1); %Fazer P1.

alfa = zeros(n);

    for k=2:n % ref. p.141 
        for i=1:(k-1)
            alfa(k,i) =   - ( (P(:,i)).'*M*W(:,k) ) / ( (P(:,i)).'*M*P(:,i)  )    ; % ref. p. 46   
        end
        soma = (P(:,:))*(alfa(k,:).');
        P(:,k) = W(:,k)+soma;
    end
        
    for k=1:n
        P(:, k) = (P(:,k))/(norm(P(:,k)));
    end
end



