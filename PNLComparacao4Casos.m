% =========================================================================
% OTIMIZAÇÃO NÃO LINEAR (PNL) - FRENAGEM REGENERATIVA FERROVIÁRIA
% Simulação Multicaso: Cargas e Velocidades Iniciais Diferentes
% =========================================================================
clear; clc; close all;

%% 1. DEFINIÇÃO DOS CENÁRIOS (CASOS)
% Vamos estruturar os 4 casos em uma matriz ou célula:
% Casos: 
% 1) Carregado, V_ini = 22.22 m/s (80 km/h)
% 2) Carregado, V_ini = 11.11 m/s (40 km/h)
% 3) Vazio,     V_ini = 22.22 m/s (80 km/h)
% 4) Vazio,     V_ini = 11.11 m/s (40 km/h)

cenarios(1) = struct('nome', 'Carregado (22.22 m/s)', 'M_prime', 37482000, 'v_ini', 22.22, 'A', 5000, 'B', 150, 'C', 1.2);
cenarios(2) = struct('nome', 'Carregado (11.11 m/s)', 'M_prime', 37482000, 'v_ini', 11.11, 'A', 5000, 'B', 150, 'C', 1.2);
cenarios(3) = struct('nome', 'Vazio (22.22 m/s)',     'M_prime', 6492000, 'v_ini', 22.22, 'A', 5000, 'B', 150,  'C', 1.2); % Valores de exemplo para trem vazio
cenarios(4) = struct('nome', 'Vazio (11.11 m/s)',     'M_prime', 6492000, 'v_ini', 11.11, 'A', 5000, 'B', 150,  'C', 1.2); % Valores de exemplo para trem vazio

num_casos = length(cenarios);

% Parâmetros comuns
a_br_max = 0.3;        % Taxa máxima de desaceleração [m/s^2]
D_total = 2000;        % Distância total de frenagem [m]
v_final = 0;           % Velocidade final [m/s]

N = 200;                                % Número de nós espaciais
ds = D_total / (N - 1);                 % Passo de distância fixo [m]
s_nodes = linspace(0, D_total, N)';     % Posição acumulada na via [m]

% Estrutura para salvar resultados de cada caso
resultados = struct();

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 50000, ...
    'MaxIterations', 1000, ...
    'ConstraintTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'Display', 'off'); % Desligado para não poluir o prompt com 4 execuções

%% 2. LOOP DE EXECUÇÃO DOS 4 CASOS
for i = 1:num_casos
    fprintf('Executando Caso %d de %d: %s...\n', i, num_casos, cenarios(i).nome);
    
    M_prime = cenarios(i).M_prime;
    v_ini = cenarios(i).v_ini;
    A_dav = cenarios(i).A; 
    B_dav = cenarios(i).B; 
    C_dav = cenarios(i).C;
    
    y_ini = v_ini^2;
    y_final = v_final^2;
    
    % Chute Inicial: perfil linear de y(s)
    y0 = linspace(y_ini, y_final, N)';
    
    % Limites de busca
    lb = zeros(N, 1);
    lb(1) = y_ini;  lb(N) = y_final;
    ub = y_ini * ones(N, 1);
    ub(1) = y_ini;  ub(N) = y_final;
    
    % Função objetivo e restrições locais
    fun_obj = @(y) obj_frenagem_pnl(y, N, ds, M_prime, A_dav, B_dav, C_dav);
    nonlcon = @(y) nonlcon_frenagem_pnl(y, N, ds, M_prime, a_br_max, A_dav, B_dav, C_dav);
    
    tic;
    [y_opt, fval, exitflag, ~] = fmincon(fun_obj, y0, [], [], [], [], lb, ub, nonlcon, options);
    t_exec = toc;
    
    % Pós-processamento do caso
    v_opt = sqrt(max(0, y_opt)); 
    v_avg = sqrt((y_opt(1:N-1) + y_opt(2:N)) / 2);
    F_d_avg = A_dav + B_dav * v_avg + C_dav * (v_avg.^2);
    F_total = M_prime * (y_opt(1:N-1) - y_opt(2:N)) / (2 * ds);
    F_disponivel = max(0, F_total - F_d_avg);
    F_e_max = max(200000, 1400000 - 15000 * v_avg);
    F_e_opt = min(F_disponivel, F_e_max);
    F_mech_opt = F_disponivel - F_e_opt;
    
    E_regenerada_kwh = (-fval) / 3.6e6;
    E_cinetica_ini_kwh = (0.5 * M_prime * v_ini^2) / 3.6e6;
    taxa_recuperacao = (E_regenerada_kwh / E_cinetica_ini_kwh) * 100;
    
    % Salvar na estrutura de resultados
    resultados(i).nome = cenarios(i).nome;
    resultados(i).v_opt = v_opt;
    resultados(i).F_e_opt = F_e_opt;
    resultados(i).F_mech_opt = F_mech_opt;
    resultados(i).E_reg = E_regenerada_kwh;
    resultados(i).taxa = taxa_recuperacao;
    resultados(i).tempo = t_exec;
end

fprintf('Otimizações concluídas com sucesso!\n\n');

%% 3. PLOTAGEM COMPARATIVA DOS 4 CASOS
figure('Name', 'Comparação dos 4 Casos - Otimização PNL', 'Position', [50, 50, 1200, 700]);

coordenadas_cores = {'b', 'r', 'g', 'm'};
estilos_linhas = {'-', '-', '--', '--'};

% Gráfico 1: Perfis de Velocidade
subplot(2,1,1);
hold on;
for i = 1:num_casos
    plot(s_nodes, resultados(i).v_opt, 'Color', coordenadas_cores{i}, ...
        'LineStyle', estilos_linhas{i}, 'LineWidth', 2, 'DisplayName', resultados(i).nome);
end
grid on;
xlabel('Distância na Via [m]');
ylabel('Velocidade [m/s]');
title('Comparação do Perfil de Velocidade Otimizado v(s) para os 4 Casos');
legend('Location', 'northeast');

% Gráfico 2: Força Elétrica (Regenerativa) Comparativa
subplot(2,1,2);
hold on;
for i = 1:num_casos
    stairs(s_nodes(1:end-1), resultados(i).F_e_opt, 'Color', coordenadas_cores{i}, ...
        'LineWidth', 1.5, 'DisplayName', [resultados(i).nome, ' (Fe)']);
end
grid on;
xlabel('Distância na Via [m]');
ylabel('Força de Frenagem Elétrica [N]');
title('Comparação da Força de Frenagem Elétrica (F_e) entre os Casos');
legend('Location', 'northeast');

%% =========================================================================
%% FUNÇÕES AUXILIARES
%% =========================================================================
function E_negativa = obj_frenagem_pnl(y, N, ds, M_prime, A_dav, B_dav, C_dav)
    v_avg = sqrt(max(0, (y(1:N-1) + y(2:N)) / 2));
    F_d = A_dav + B_dav * v_avg + C_dav * (v_avg.^2);
    F_total = M_prime * (y(1:N-1) - y(2:N)) / (2 * ds);
    
    F_disponivel = max(0, F_total - F_d);
    F_e_max = max(200000, 1400000 - 15000 * v_avg);
    
    F_e = min(F_disponivel, F_e_max);
    E_negativa = -sum(F_e * ds);
end

function [c, ceq] = nonlcon_frenagem_pnl(y, N, ds, M_prime, a_br_max, A_dav, B_dav, C_dav)
    c_mono = y(2:N) - y(1:N-1);
    c_accel = (y(1:N-1) - y(2:N)) - (2 * a_br_max * ds);
    
    v_avg = sqrt(max(0, (y(1:N-1) + y(2:N)) / 2));
    F_d = A_dav + B_dav * v_avg + C_dav * (v_avg.^2);
    F_total = M_prime * (y(1:N-1) - y(2:N)) / (2 * ds);
    c_no_motor = -(F_total - F_d);
    
    c = [c_mono; c_accel; c_no_motor];
    ceq = [];
end