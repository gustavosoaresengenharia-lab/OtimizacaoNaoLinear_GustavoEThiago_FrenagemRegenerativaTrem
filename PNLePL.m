% =========================================================================
% COMPARAÇÃO: OTIMIZAÇÃO LINEAR (PL) vs. NÃO LINEAR (PNL)
% Frenagem Regenerativa Ferroviária
% =========================================================================
clear; clc; close all;

%% =========================================================================
%% 1. PARÂMETROS GERAIS (MODIFIQUE AQUI PARA VARIAR OS TESTES)
%% =========================================================================
M_prime  = 37482000;   % Massa efetiva [kg]: 37482000 (carregado) / 6492000 (vazio)
v_ini    = 11.11;      % Velocidade inicial [m/s]: 22.22 (80km/h) / 11.11 (40km/h)
v_final  = 0;          % Velocidade final [m/s]
D_total  = 2000;       % Distância total disponível para frenagem [m]
a_br_max = 0.3;        % Taxa máxima de desaceleração [m/s^2]

% Coeficientes de Davis (Resistência ao avanço)
A_dav = 5000; 
B_dav = 150; 
C_dav = 1.2;

% Discretizações (N)
N_pl  = 400;           % Número de nós para PL (discretização de velocidade)
N_pnl = 200;           % Número de nós para PNL (discretização de espaço)

fprintf('Iniciando Comparação PL vs PNL...\n');
fprintf('Cenário: Massa = %.0f kg | V_ini = %.2f m/s | Distância = %.0f m\n\n', M_prime, v_ini, D_total);

%% =========================================================================
%% 2. MODELO DE PROGRAMAÇÃO LINEAR (PL) - LU ET AL. (Discretização V)
%% =========================================================================
fprintf('Executando Programação Linear (PL)...\n');
tic;

v_nodes = linspace(v_ini, v_final, N_pl); 
num_intervals = N_pl - 1;
num_vars = 2 * num_intervals;

f_pl = [zeros(num_intervals, 1); -ones(num_intervals, 1)];

d_min = zeros(num_intervals, 1);
F_d_avg_pl = zeros(num_intervals, 1);
F_e_max_d = zeros(num_intervals, 1);
v_avg_pl_array = zeros(num_intervals, 1);

for i = 1:num_intervals
    v_avg = (v_nodes(i) + v_nodes(i+1)) / 2;
    v_avg_pl_array(i) = v_avg;
    d_min(i) = (v_nodes(i)^2 - v_nodes(i+1)^2) / (2 * a_br_max); 
    F_d_avg_pl(i) = A_dav + B_dav * v_avg + C_dav * (v_avg^2);     
    F_e_max = max(200000, 1400000 - 15000 * v_avg);
    F_e_max_d(i) = F_e_max;
end

A_ineq = []; b_ineq = [];
for i = 1:num_intervals
    % Restricao 1: d_i >= d_min
    row_d = zeros(1, num_vars); row_d(i) = -1;
    A_ineq = [A_ineq; row_d];
    b_ineq = [b_ineq; -d_min(i)];
    
    % Restricao 2: Limite dinamico
    row_e1 = zeros(1, num_vars); row_e1(i) = F_d_avg_pl(i); row_e1(num_intervals + i) = 1;
    A_ineq = [A_ineq; row_e1];
    b_ineq = [b_ineq; 0.5 * M_prime * (v_nodes(i)^2 - v_nodes(i+1)^2)];
    
    % Restricao 3: Limite do inversor
    row_e2 = zeros(1, num_vars); row_e2(i) = -F_e_max_d(i); row_e2(num_intervals + i) = 1;
    A_ineq = [A_ineq; row_e2];
    b_ineq = [b_ineq; 0];
end

Aeq = zeros(1, num_vars); Aeq(1:num_intervals) = 1;
beq = D_total;
lb_pl = zeros(num_vars, 1); 
ub_pl = inf(num_vars, 1);

options_pl = optimoptions('linprog', 'Display', 'off');
[x_opt_pl, fval_pl, exitflag_pl] = linprog(f_pl, A_ineq, b_ineq, Aeq, beq, lb_pl, ub_pl, options_pl);

t_exec_pl = toc;

% Pós-processamento PL
d_otimo_pl = x_opt_pl(1:num_intervals);
E_eb_otimo_pl = x_opt_pl(num_intervals+1:end);
distancia_acumulada_pl = [0; cumsum(d_otimo_pl)];
F_eb_otimo_pl = E_eb_otimo_pl ./ max(d_otimo_pl, eps); 

% Cálculo de Tempo PL (dt = d / v_avg)
dt_pl = d_otimo_pl ./ v_avg_pl_array;
tempo_frenagem_pl = sum(dt_pl);

rbe_total_kwh_pl = (-fval_pl) / 3.6e6;


%% =========================================================================
%% 3. MODELO DE PROGRAMAÇÃO NÃO LINEAR (PNL) - (Discretização S)
%% =========================================================================
fprintf('Executando Programação Não Linear (PNL)...\n');
tic;

ds = D_total / (N_pnl - 1);                 
s_nodes_pnl = linspace(0, D_total, N_pnl)';     
y_ini = v_ini^2; y_final = v_final^2;
y0 = linspace(y_ini, y_final, N_pnl)';

lb_pnl = zeros(N_pnl, 1); lb_pnl(1) = y_ini; lb_pnl(N_pnl) = y_final;
ub_pnl = y_ini * ones(N_pnl, 1); ub_pnl(1) = y_ini; ub_pnl(N_pnl) = y_final;

fun_obj = @(y) obj_frenagem_pnl(y, N_pnl, ds, M_prime, A_dav, B_dav, C_dav);
nonlcon = @(y) nonlcon_frenagem_pnl(y, N_pnl, ds, M_prime, a_br_max, A_dav, B_dav, C_dav);

options_pnl = optimoptions('fmincon', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 50000, 'MaxIterations', 1000, ...
    'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'Display', 'off');

[y_opt_pnl, fval_pnl, exitflag_pnl] = fmincon(fun_obj, y0, [], [], [], [], lb_pnl, ub_pnl, nonlcon, options_pnl);

t_exec_pnl = toc;

% Pós-processamento PNL
v_opt_pnl = sqrt(max(0, y_opt_pnl)); 
v_avg_pnl = sqrt((y_opt_pnl(1:N_pnl-1) + y_opt_pnl(2:N_pnl)) / 2);
dt_pnl = ds ./ v_avg_pnl;                 
tempo_frenagem_pnl = sum(dt_pnl);         

F_d_avg_pnl = A_dav + B_dav * v_avg_pnl + C_dav * (v_avg_pnl.^2);
F_total_pnl = M_prime * (y_opt_pnl(1:N_pnl-1) - y_opt_pnl(2:N_pnl)) / (2 * ds);
F_disponivel_pnl = max(0, F_total_pnl - F_d_avg_pnl);
F_e_max_pnl = max(200000, 1400000 - 15000 * v_avg_pnl);
F_e_opt_pnl = min(F_disponivel_pnl, F_e_max_pnl);

rbe_total_kwh_pnl = (-fval_pnl) / 3.6e6;


%% =========================================================================
%% 4. COMPARAÇÃO DE RESULTADOS (CONSOLE)
%% =========================================================================
E_cinetica_inicial_kwh = (0.5 * M_prime * v_ini^2) / 3.6e6;
taxa_recuperacao_pl = (rbe_total_kwh_pl / E_cinetica_inicial_kwh) * 100;
taxa_recuperacao_pnl = (rbe_total_kwh_pnl / E_cinetica_inicial_kwh) * 100;

disp('=================================================================');
disp('                TABELA COMPARATIVA: PL vs PNL                    ');
disp('=================================================================');
fprintf('Energia Cinética Inicial:        %.2f kWh\n', E_cinetica_inicial_kwh);
disp('-----------------------------------------------------------------');
fprintf('MÉTRICA                          | PL (Linear)  | PNL (Não Linear)\n');
disp('-----------------------------------------------------------------');
fprintf('Energia Regenerada [kWh]         | %-12.2f | %-12.2f\n', rbe_total_kwh_pl, rbe_total_kwh_pnl);
fprintf('Taxa de Recuperação [%%]          | %-12.2f | %-12.2f\n', taxa_recuperacao_pl, taxa_recuperacao_pnl);
fprintf('Tempo de Frenagem do Trem [s]    | %-12.2f | %-12.2f\n', tempo_frenagem_pl, tempo_frenagem_pnl);
fprintf('Tempo de Execução do Código [s]  | %-12.4f | %-12.4f\n', t_exec_pl, t_exec_pnl);
fprintf('Status de Saída (Exitflag)       | %-12d | %-12d\n', exitflag_pl, exitflag_pnl);
disp('=================================================================');

%% =========================================================================
%% 5. PLOTAGEM COMPARATIVA
%% =========================================================================
figure('Name', 'Comparação de Otimização: PL vs PNL', 'Position', [100, 100, 900, 600]);

% Gráfico 1: Perfil de Velocidade
subplot(2, 1, 1);
plot(distancia_acumulada_pl, v_nodes, 'r--', 'LineWidth', 2); hold on;
plot(s_nodes_pnl, v_opt_pnl, 'b-', 'LineWidth', 2);
grid on;
xlabel('Distância na Via [m]');
ylabel('Velocidade [m/s]');
title('Comparação do Perfil de Velocidade: PL vs PNL');
legend('PL (Discretização de Velocidade)', 'PNL (Discretização Espacial)', 'Location', 'northeast');

% Gráfico 2: Força Elétrica (Regenerativa)
subplot(2, 1, 2);
stairs([distancia_acumulada_pl(1:end-1); distancia_acumulada_pl(end)], [F_eb_otimo_pl; F_eb_otimo_pl(end)], 'r--', 'LineWidth', 1.5); hold on;
stairs(s_nodes_pnl(1:end-1), F_e_opt_pnl, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Distância na Via [m]');
ylabel('Força de Frenagem Elétrica [N]');
title('Comparação da Força Elétrica Regenerativa');
legend('Força Elétrica - PL', 'Força Elétrica - PNL', 'Location', 'northeast');

%% =========================================================================
%% FUNÇÕES AUXILIARES (APENAS PARA O PNL)
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