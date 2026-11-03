% =========================================================================
% OTIMIZAÇÃO NÃO LINEAR (PNL) - FRENAGEM REGENERATIVA FERROVIÁRIA
% Formulacao por Variavel de Estado y(s) = v(s)^2 com Frenagem Mista
% =========================================================================
clear; clc; close all;

%% 1. PARÂMETROS DO SISTEMA (TREM CARREGADO - CASO A)
M_prime = 6492000;    % Massa efetiva do trem [kg]
a_br_max = 0.3;        % Taxa máxima de desaceleração [m/s^2]
D_total = 2000;        % Distância total de frenagem [m]
v_ini = 30/3.6;            % Velocidade inicial [m/s]
v_final = 0;           % Velocidade final [m/s]

% Coeficientes de Davis (Resistência ao avanço)
A_dav = 5000; B_dav = 150; C_dav = 1.2;

%% 2. DISCRETIZAÇÃO ESPACIAL E VARIÁVEIS DE DECISÃO
N = 200;                                % Número de nós espaciais
ds = D_total / (N - 1);                 % Passo de distância fixo [m]
s_nodes = linspace(0, D_total, N)';     % Posição acumulada na via [m]

% Variável de decisão: y(s) = v(s)^2 nos N nós espaciais [m^2/s^2]
y_ini = v_ini^2;
y_final = v_final^2;

% Chute Inicial: perfil linear de y(s) (100% factível fisicamente)
y0 = linspace(y_ini, y_final, N)';

% Limites de busca (lb <= y <= ub)
lb = zeros(N, 1);
lb(1) = y_ini;          % Condição de contorno inicial fixa
lb(N) = y_final;        % Condição de contorno final fixa

ub = y_ini * ones(N, 1);
ub(1) = y_ini;
ub(N) = y_final;

%% 3. EXECUÇÃO DO SOLVER fmincon
fun_obj = @(y) obj_frenagem_pnl(y, N, ds, M_prime, A_dav, B_dav, C_dav);

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 50000, ...
    'MaxIterations', 1000, ...
    'ConstraintTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, ...
    'Display', 'iter');

fprintf('Iniciando Otimização Não Linear (PNL - Formulação Mista)...\n');
tic;
[y_opt, fval, exitflag, output] = fmincon(...
    fun_obj, y0, [], [], [], [], lb, ub, ...
    @(y) nonlcon_frenagem_pnl(y, N, ds, M_prime, a_br_max, A_dav, B_dav, C_dav), ...
    options);
t_exec = toc;

%% 4. PÓS-PROCESSAMENTO E EXIBIÇÃO DE RESULTADOS
v_opt = sqrt(max(0, y_opt)); 

v_avg = sqrt((y_opt(1:N-1) + y_opt(2:N)) / 2);
F_d_avg = A_dav + B_dav * v_avg + C_dav * (v_avg.^2);
F_total = M_prime * (y_opt(1:N-1) - y_opt(2:N)) / (2 * ds);

F_disponivel = max(0, F_total - F_d_avg);
F_e_max = max(200000, 1400000 - 15000 * v_avg);
F_e_opt = min(F_disponivel, F_e_max);
F_mech_opt = F_disponivel - F_e_opt;

E_regenerada_joules = -fval;
E_regenerada_kwh = E_regenerada_joules / 3.6e6;
E_cinetica_ini_kwh = (0.5 * M_prime * v_ini^2) / 3.6e6;
taxa_recuperacao = (E_regenerada_kwh / E_cinetica_ini_kwh) * 100;

disp('==================================================');
disp('      RESULTADOS DA OTIMIZAÇÃO NÃO LINEAR (PNL)   ');
disp('==================================================');
disp(['Exitflag: ', num2str(exitflag)]);
disp(['Tempo de Execução: ', num2str(t_exec, '%.4f'), ' s']);
disp(['Energia Cinética Inicial: ', num2str(E_cinetica_ini_kwh, '%.2f'), ' kWh']);
disp(['Energia Elétrica Regenerada: ', num2str(E_regenerada_kwh, '%.2f'), ' kWh']);
disp(['Taxa de Recuperação: ', num2str(taxa_recuperacao, '%.2f'), ' %']);

%% 5. PLOTAGEM DOS RESULTADOS
figure('Name', 'Otimização PNL - Frenagem Regenerativa e Mecânica', 'Position', [100, 100, 800, 600]);

subplot(2,1,1);
plot(s_nodes, v_opt, 'b-', 'LineWidth', 2);
grid on;
xlabel('Distância na Via [m]');
ylabel('Velocidade [m/s]');
title('Perfil Não Linear de Velocidade Otimizado v(s)');

subplot(2,1,2);
stairs(s_nodes(1:end-1), F_e_opt, 'g-', 'LineWidth', 1.5); hold on;
stairs(s_nodes(1:end-1), F_mech_opt, 'r--', 'LineWidth', 1.2);
grid on;
xlabel('Distância na Via [m]');
ylabel('Força de Frenagem [N]');
title('Distribuição de Esforços: Elétrico (Regenerativo) vs. Mecânico');
legend('Força Elétrica (F_e)', 'Força Mecânica (F_{mech})', 'Location', 'northeast');

%% =========================================================================
%% FUNÇÕES AUXILIARES
%% =========================================================================
function E_negativa = obj_frenagem_pnl(y, N, ds, M_prime, A_dav, B_dav, C_dav)
    v_avg = sqrt(max(0, (y(1:N-1) + y(2:N)) / 2));
    F_d = A_dav + B_dav * v_avg + C_dav * (v_avg.^2);
    F_total = M_prime * (y(1:N-1) - y(2:N)) / (2 * ds);
    
    F_disponivel = max(0, F_total - F_d);
    F_e_max = max(200000, 1400000 - 15000 * v_avg);
    
    % Saturação física no teto do inversor
    F_e = min(F_disponivel, F_e_max);
    
    E_negativa = -sum(F_e * ds);
end

function [c, ceq] = nonlcon_frenagem_pnl(y, N, ds, M_prime, a_br_max, A_dav, B_dav, C_dav)
    % 1. Monotonicidade (velocidade decrescente)
    c_mono = y(2:N) - y(1:N-1);
    
    % 2. Limite mecânico de desaceleração
    c_accel = (y(1:N-1) - y(2:N)) - (2 * a_br_max * ds);
    
    % 3. Composição em frenagem (sem força de tração motora)
    v_avg = sqrt(max(0, (y(1:N-1) + y(2:N)) / 2));
    F_d = A_dav + B_dav * v_avg + C_dav * (v_avg.^2);
    F_total = M_prime * (y(1:N-1) - y(2:N)) / (2 * ds);
    c_no_motor = -(F_total - F_d);
    
    c = [c_mono; c_accel; c_no_motor];
    ceq = [];
end