% Cálculo dos vetores com grau de liberdade livre
Mt = Trator(E);

% Cálculo da matriz constante K
K = Mt(:,7:8)
K(7,1)=1
K(8,2)=1

% O vetor R corresponde a X7 e X8