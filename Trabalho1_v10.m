%% Trabalho 1 - 2021/1
% EPQB/UFRJ
% EQE 703 - Métodos Matemáticos Aplicados
% Por Gustavo Caldas (gustavo.caldas@eq.ufrj.br) e Oscar Chamberlain (ochamberlain@eq.ufrj.br)

A = [ 1 1 1 1 1 1 1 1;
    1 2 3 4 5 6 7 8;
    1 3 6 10 15 21 28 36;
    1 4 10 20 35 56 84 120;
    1 5 15 35 70 126 210 330;
    1 6 21 56 126 252 462 792;
    1 7 28 84 210 462 924 1716;
    1 8 36 120 330 792 1716 3432 ];
%% Item ii
P = conjugado_schmidt(A);
%% Item iii
D = (P.'*A)*P;

%% Item iv
U = diagonalizacao_forma_quadratica(A);
Q = U.'*A*U;
autovalores = eig(A);
%% Item v - caráter da Matriz Hessiana

%H = 2*(A*X)*((A*X).')+(X.'*A*X)*A + A;

%H*P = P*lambda

% H*P = [2*(A*X)*((A*X)')+(X'*A*X)*A + A]*P = lambda*P
% Considerando que X'AX é um constante K e que (AX)' = X'A'=X'A
% H*P = [2*(A*X)X'+ K + I]AP = lambda*P
% O caracter de A definirá o caracter da hessiana.

% Para conhecer o caráter da matriz Hessiana, precisamos obter o autovalor.

%% Item vi - Newton-Raphson
X0  = [0;0;0;0;0;0;0;0];
B=[11111;
22222;
33333;
44444;
44444;
33333;
22222;
11111];


% X é valor X^(n), onde n é a última iteração
% Iterations é o número de iterações
% path é o tensor histórico com os valores de X^(k): 8 x Iterations
% Spath é o tensor histórico de S^k: 8 x Iterations
[X,path,iterations,Spath] = newton_raphson(X0,A,B);

%Calculando os valores de F
C = 1e5;
F = (1/4)*(X.'*A*X)^2 +(0.5)*(X.'*A*X)+B.'*X + C; % último valor

% Cálculo do valor histórico de f
x = path;
f = zeros(size(x,2),1);
for i=1:size(x,2)
    f(i) = (1/4)*(x(:,i).'*A*x(:,i))^2 +(0.5)*(x(:,i).'*A*x(:,i))+B.'*x(:,i) + C;
end

% Plotando!
figure(1);
quiver3(x(1,:),x(2,:),x(3,:),Spath(1,:),Spath(2,:),Spath(3,:),0.2);
grid on
hold on
plot3(x(1,:),x(2,:),x(3,:),'o','MarkerFaceColor','r');
xlabel({'$X_1$'},'Interpreter','latex');
ylabel({'$X_2$'},'Interpreter','latex');
zlabel({'$X_3$'},'Interpreter','latex');
title({'Trajetória X_1 vs. X_2 vs. X_3 por Método de Newton-Raphson'},'Interpreter','tex');
print('x1x2x3-NR.jpg','-djpeg','-r300')

figure(2);
quiver3(x(4,:),x(5,:),x(6,:),Spath(4,:),Spath(5,:),Spath(6,:),0.2);
grid on
hold on
plot3(x(4,:),x(5,:),x(6,:),'o','MarkerFaceColor','r')
xlabel({'$X_4$'},'Interpreter','latex');
ylabel({'$X_5$'},'Interpreter','latex');
zlabel({'$X_6$'},'Interpreter','latex');
title({'Trajetória X_4 vs. X_5 vs. X_6 por Método de Newton-Raphson'},'Interpreter','tex');
print('x4x5x6-NR.jpg','-djpeg','-r300')

figure(3);
quiver(x(7,:),x(8,:),Spath(7,:),Spath(8,:),0.1)
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
grid on
hold on
plot(x(7,:),x(8,:),'o','MarkerFaceColor','r')
title({'Trajetória X_7 vs. X_8 por Método de Newton-Raphson'},'Interpreter','tex');
print('x7vsx8-NR.jpg','-djpeg','-r300')

figure(4);

plot3(x(7,:),x(8,:),f,'o','MarkerFaceColor','r')
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
zlabel({'F'},'Interpreter','latex');
title({'Trajetória X_7 vs. X_8 vs. Valor F por Método de Newton-Raphson'},'Interpreter','tex');
print('x7vsx8vsF-NR.jpg','-djpeg','-r300')

%Calculando o traço de H;
H = 2*(A*X)*((A*X).')+(X.'*A*X)*A + A;
traco = trace(H);

%% Item vii - Método de Powell
%X0  = [100;10;101;10;10;10;200;5000];
[X_powell,path_powell,iterations_powell,Spath_powell] = powell(X0,A,B);

%Calculando os valores de F
C = 1e5;
F_powell = (1/4)*(X_powell.'*A*X_powell)^2 +(0.5)*(X_powell.'*A*X_powell)+B.'*X_powell + C; % último valor

% Cálculo do valor histórico de f
x_powell = path_powell;
f_powell = zeros(size(x_powell,2),1);
for i=1:size(x_powell,2)
    f_powell(i) = (1/4)*(x_powell(:,i).'*A*x_powell(:,i))^2 +(0.5)*(x_powell(:,i).'*A*x_powell(:,i))+B.'*x_powell(:,i) + C;
end

% Plotando!
figure(5);
quiver3(x_powell(1,:),x_powell(2,:),x_powell(3,:),Spath_powell(1,:),Spath_powell(2,:),Spath_powell(3,:));
grid on
hold on
plot3(x_powell(1,:),x_powell(2,:),x_powell(3,:),'o','MarkerFaceColor','r');
%contour(x(1,:),x(2,:),f);
%quiver(Spath(1,:),Spath(2,:),Spath(3,:))
xlabel({'$X_1$'},'Interpreter','latex');
ylabel({'$X_2$'},'Interpreter','latex');
zlabel({'$X_3$'},'Interpreter','latex');
title({'Trajetória X_1 vs. X_2 vs. X_3 por Método de Powell'},'Interpreter','tex');
print('x1x2x3-powell.jpg','-djpeg','-r300')

figure(6);
quiver3(x_powell(4,:),x_powell(5,:),x_powell(6,:),Spath_powell(4,:),Spath_powell(5,:),Spath_powell(6,:),1.5);
grid on
hold on
plot3(x_powell(4,:),x_powell(5,:),x_powell(6,:),'o','MarkerFaceColor','r')
xlabel({'$X_4$'},'Interpreter','latex');
ylabel({'$X_5$'},'Interpreter','latex');
zlabel({'$X_6$'},'Interpreter','latex');
title({'Trajetória X_4 vs. X_5 vs. X_6 por Método de Powell'},'Interpreter','tex');
print('x4x5x6-powell.jpg','-djpeg','-r300')

figure(7);

quiver(x_powell(7,:),x_powell(8,:),Spath_powell(7,:),Spath_powell(8,:),1.5)
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
grid on
hold on
plot(x_powell(7,:),x_powell(8,:),'o','MarkerFaceColor','r')
title({'Trajetória X_7 vs. X_8 por Método de Powell'},'Interpreter','tex');
print('x7vsx8-powell.jpg','-djpeg','-r300')

figure(8);

plot3(x_powell(7,:),x_powell(8,:),f_powell,'o','MarkerFaceColor','r')
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
zlabel({'F'},'Interpreter','latex');
title({'Trajetória X_1 vs. X_2 vs. Valor F por Método de Powell'},'Interpreter','tex');
print('x7vsx8vsF-powell.jpg','-djpeg','-r300')
%% Item viii - Espaço Reduzido

E = [4 0 4 0 -3 -4 3 1
-2 -4 2 3 -2 2 0 -1
1 3 -3 -4 -3 0 -2 3
0 0 0 -1 1 4 1 0
3 1 4 3 -2 0 3 2
2 2 4 -4 -3 0 -4 0 ];

% Cálculo dos vetores com grau de liberdade livre
Mt = Trator(E);

% Cálculo da matriz constante K
K = -Mt(:,7:8);
K(7,1)=1;
K(8,2)=1;

% Plotando a projeção de F em função de F(X(R))
R(:,1) = linspace(-100,100);
R(:,2) = linspace(-100,100);

%Cálculo de X no Espaço Reduzido
X_reduz = K*R.';

% F no espaço reduzido
F_reduz = (1/4)*(X_reduz.'*A*X_reduz)^2 +(0.5)*(X_reduz.'*A*X_reduz)+B.'*X_reduz + C;

%Projeção no espaço reduzido

%Superfície
figure(9);

surf(R(:,1),R(:,2),F_reduz);
colorbar
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
zlabel({'F'},'Interpreter','latex');
title({'Superfície X_7 vs. X_8 vs. Valor F no Espaço Reduzido'},'Interpreter','tex');
print('x7vsx8-surface_reduced.jpg','-djpeg','-r300')

% Contorno
figure(10);

contour(R(:,1),R(:,2),F_reduz);
colorbar
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
%zlabel({'F'},'Interpreter','latex');
title({'Contorno X_7 vs. X_8 vs. Valor F no Espaço Reduzido'},'Interpreter','tex');
print('x7vsx8vsF-contour_reduced.jpg','-djpeg','-r300')

%% Item ix - Espaço Reduzido com Newton-Rapson

X0  = [0;0;0;0;0;0;1;1];
% X é valor X^(n), onde n é a última iteração
% Iterations é o número de iterações
% path é o tensor histórico com os valores de X^(k): 8 x Iterations
% Spath é o tensor histórico de S^k: 8 x Iterations
[X_nr_er,path_nr_er,iterations_nr_er,Spath_nr_er] = newton_raphson_9(X0,A,B,K);

%Calculando os valores de F
C = 1e5;
F_nr_er = (1/4)*(X_nr_er.'*A*X_nr_er)^2 +(0.5)*(X_nr_er.'*A*X_nr_er)+B.'*X_nr_er + C; % último valor

% Cálculo do valor histórico de f
x_nr_er = path_nr_er;
f_nr_er  = zeros(size(x_nr_er,2),1);
for i=1:size(x_nr_er,2)
    f_nr_er(i) = (1/4)*(x_nr_er(:,i).'*A*x_nr_er(:,i))^2 +(0.5)*(x_nr_er(:,i).'*A*x_nr_er(:,i))+B.'*x_nr_er(:,i) + C;
end

i_nr_er=size(x_nr_er,2);

% Plotando!
figure(11);
quiver3(x_nr_er(1,:),x_nr_er(2,:),x_nr_er(3,:),Spath_nr_er(1,:),Spath_nr_er(2,:),Spath_nr_er(3,:),0.2);
grid on
hold on
plot3(x_nr_er(1,:),x_nr_er(2,:),x_nr_er(3,:),'o','MarkerFaceColor','r');
xlabel({'$X_1$'},'Interpreter','latex');
ylabel({'$X_2$'},'Interpreter','latex');
zlabel({'$X_3$'},'Interpreter','latex');
title({'Trajetória X_1 vs. X_2 vs. X_3 no Esp. Reduzido por Método de Newton-Raphson'},'Interpreter','tex');
print('x1x2x3-NR_er.jpg','-djpeg','-r300')

figure(12);

quiver3(x_nr_er(4,:),x_nr_er(5,:),x_nr_er(6,:),Spath_nr_er(4,:),Spath_nr_er(5,:),Spath_nr_er(6,:),0.2);
grid on
hold on
plot3(x_nr_er(4,:),x_nr_er(5,:),x_nr_er(6,:),'o','MarkerFaceColor','r')
xlabel({'$X_4$'},'Interpreter','latex');
ylabel({'$X_5$'},'Interpreter','latex');
zlabel({'$X_6$'},'Interpreter','latex');
title({'Trajetória X_4 vs. X_5 vs. X_6 no Esp. Reduzido por Método de Newton-Raphson'},'Interpreter','tex');
print('x4x5x6-NR_er.jpg','-djpeg','-r300')

figure(13);

quiver(x_nr_er(7,:),x_nr_er(8,:),Spath_nr_er(7,:),Spath_nr_er(8,:),0.2)
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
grid on
hold on
plot(x_nr_er(7,:),x_nr_er(8,:),'o','MarkerFaceColor','r')
title({'Trajetória X_7 vs. X_8 no Espaço Reduzido por Método de Newton-Raphson'},'Interpreter','tex');
print('x7vsx8-NR_er.jpg','-djpeg','-r300')

% Plotando a projeção de F em função de F(X(R))
R_nr(:,1) = linspace(0,4);
R_nr(:,2) = linspace(-2,2);
%Cálculo de X no Espaço Reduzido
X_reduz_nr = K*R_nr.'; %valores de X para projeção
% F no espaço reduzido
F_reduz_nr = (1/4)*(X_reduz_nr.'*A*X_reduz_nr)^2 +(0.5)*(X_reduz_nr.'*A*X_reduz_nr)+B.'*X_reduz_nr + C;

%Projeção no espaço reduzido

figure(14);
surf(R_nr(:,1),R_nr(:,2),F_reduz_nr);
colorbar
grid on
hold on
plot3(x_nr_er(7,:),x_nr_er(8,:),f_nr_er,'o','MarkerFaceColor','r')
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
zlabel({'F'},'Interpreter','latex');
title({'Trajetória X_7 vs. X_8 vs.F no Esp. Red. por Método de NR'},'Interpreter','tex');
print('x7vsx8vsF-NR_er.jpg','-djpeg','-r300')
%% Item x - Espaço Reduzido com Método de Powell

E = [ 4 0 4 0 -3 -4 3 1
-2 -4 2 3 -2 2 0 -1
1 3 -3 -4 -3 0 -2 3
0 0 0 -1 1 4 1 0
3 1 4 3 -2 0 3 2
2 2 4 -4 -3 0 -4 0 ];

% Cálculo dos vetores com grau de liberdade livre
Mt = Trator(E);

% Cálculo da matriz constante K
K = -Mt(:,7:8);
K(7,1)=1;
K(8,2)=1;

%X0  = -[100;10;101;10;10;10;200;5000];
X0 = 0.016* [1 1 1 1 1 1 1 1]';% 0.3536
[X_powell_er,path_powell_er,iterations_powell_er,Spath_powell_er] = powell_10a(X0,A,B,K);

%Calculando os valores de F
C = 1e5;
F_powell_er = (1/4)*(X_powell_er.'*A*X_powell_er)^2 +(0.5)*(X_powell_er.'*A*X_powell_er)+B.'*X_powell_er + C; % último valor

% Cálculo do valor histórico de f
x_powell_er = path_powell_er;
f_powell_er = zeros(size(x_powell_er,2),1);
for i=1:size(x_powell_er,2)
    f_powell_er(i) = (1/4)*(x_powell_er(:,i).'*A*x_powell_er(:,i))^2 +(0.5)*(x_powell_er(:,i).'*A*x_powell_er(:,i))+B.'*x_powell_er(:,i) + C;
end

i_powell_er=size(x_powell_er,2);

% Plotando!
figure(15);
quiver3(x_powell_er(1,:),x_powell_er(2,:),x_powell_er(3,:),Spath_powell_er(1,:),Spath_powell_er(2,:),Spath_powell_er(3,:),0.2);
grid on
hold on
plot3(x_powell_er(1,:),x_powell_er(2,:),x_powell_er(3,:),'o','MarkerFaceColor','r');
xlabel({'$X_1$'},'Interpreter','latex');
ylabel({'$X_2$'},'Interpreter','latex');
zlabel({'$X_3$'},'Interpreter','latex');
title({'Trajetória X_1 vs. X_2 vs. X_3 no Espaço Reduzido por Método de Powell'},'Interpreter','tex');
print('x1x2x3-powell_er.jpg','-djpeg','-r300')

figure(16);
quiver3(x_powell_er(4,:),x_powell_er(5,:),x_powell_er(6,:),Spath_powell_er(4,:),Spath_powell_er(5,:),Spath_powell_er(6,:),0.2);
grid on
hold on
plot3(x_powell_er(4,:),x_powell_er(5,:),x_powell_er(6,:),'o','MarkerFaceColor','r')
xlabel({'$X_4$'},'Interpreter','latex');
ylabel({'$X_5$'},'Interpreter','latex');
zlabel({'$X_6$'},'Interpreter','latex');
title({'Trajetória X_4 vs. X_5 vs. X_6 no Espaço Reduzido por Método de Powell'},'Interpreter','tex');
print('x4x5x6-powell_er.jpg','-djpeg','-r300')

figure(17);

quiver(x_powell_er(7,:),x_powell_er(8,:),Spath_powell_er(7,:),Spath_powell_er(8,:),0.2)
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
grid on
hold on
plot(x_powell_er(7,:),x_powell_er(8,:),'o','MarkerFaceColor','r')
title({'Trajetória X_7 vs. X_8 no Espaço Reduzido por Método de Powell'},'Interpreter','tex');
print('x7vsx8-powell_er.jpg','-djpeg','-r300')

% Plotando a projeção de F em função de F(X(R))
R_powell(:,1) = linspace(-1,5);
R_powell(:,2) = linspace(-1,1);
%Cálculo de X no Espaço Reduzido
X_reduz_powell = K*R_powell.'; %valores de X para projeção
% F no espaço reduzido
F_reduz_powell = (1/4)*(X_reduz_powell.'*A*X_reduz_powell)^2 +(0.5)*(X_reduz_powell.'*A*X_reduz_powell)+B.'*X_reduz_powell + C;

figure(18);
surf(R_powell(:,1),R_powell(:,2),F_reduz_powell);
colorbar
grid on
hold on
plot3(x_powell_er(7,:),x_powell_er(8,:), f_powell_er,'o','MarkerFaceColor','r')
xlabel({'$X_7$'},'Interpreter','latex');
ylabel({'$X_8$'},'Interpreter','latex');
zlabel({'F'},'Interpreter','latex');
title({'Trajetória X_7 vs. X_8 vs. Valor F no Esp. Red. por Método de Powell'},'Interpreter','tex');
print('x7vsx8vsF-powell_er.jpg','-djpeg','-r300')