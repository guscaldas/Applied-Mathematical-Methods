% Metodologia para transformação da matriz
% Considerando a Matriz M (nxm)
function [Mt] = Trator(M)
Mt = M;

%Tamanho da matriz
t = size(M);
L=t(1);
C=t(2);

%Transformação das linhas
for i=1:L
    Mt(i,:) = Mt(i,:)/Mt(i,i);
    for k=1:L
        if k~=i
            Mt(k,:)= Mt(k,:)-Mt(i,:)*Mt(k,i);
        end
    end   
end
end