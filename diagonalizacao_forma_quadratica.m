%% Script file: diagonalizacao_forma_quadratica.m
% Objetivo: Este programa executa a Diagonalização de Forma Quadrática com Transformação de
%Congruência de uma matriz simétrica pelo método de Schmidt.
% Referência: S. Chapman, "Programação para Engenheiros"

function U = diagonalizacao_forma_quadratica(M)
n = length(M);
I = eye(n); % ref. p. 29
U = I;
alfa = zeros(n); %matriz do alfa é de zeros
soma = 0;
    for k=2:n % ref. p.141 
        for i=1:(k-1)
            %Calculando os alfas
            alfa(k,i) = -( (U(:,i)).'*M*I(:,k) )/((U(:,i)).'*M*U(:,i)); % ref. p. 46
        end
    %Calculando a soma
    soma = U(:,:)*alfa(k,:).';
    U(:,k) = U(:,k)+soma;
    end    
       
end



