clc;
clear all;
close all;

SNR_list = -10:1:40;
N_list = [1 2 4 8 16];
M = 4;
SER_ana_lst = zeros(length(N_list), length(SNR_list));
SER_mc_lst = zeros(length(N_list), length(SNR_list));

progress = 0;
for i=1:1:length(SNR_list)
    for j=1:1:length(N_list)
        [SER_ana_lst(j, i), SER_mc_lst(j, i)] = runSystem( N_list(j), M, SNR_list(i));
        progress = (progress+1);
        progress_perc = progress/(length(SNR_list)*length(N_list))*100
    end
end

figure;

semilogy(SNR_list, SER_mc_lst(1,:), 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 10, 'Color', 'black', 'Linewidth', 1);hold on;
semilogy(SNR_list, SER_ana_lst(1,:), 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'black', 'Linewidth', 1);

semilogy(SNR_list, SER_mc_lst(2,:), 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 10, 'Color', 'red', 'Linewidth', 1);
semilogy(SNR_list, SER_ana_lst(2,:), 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'red', 'Linewidth', 1);

semilogy(SNR_list, SER_mc_lst(3,:), 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 10, 'Color', 'blue', 'Linewidth', 1);
semilogy(SNR_list, SER_ana_lst(3,:), 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'blue', 'Linewidth', 1);

semilogy(SNR_list, SER_mc_lst(4,:), 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 10, 'Color', 'magenta', 'Linewidth', 1);hold on;
semilogy(SNR_list, SER_ana_lst(4,:), 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'magenta', 'Linewidth', 1);

semilogy(SNR_list, SER_mc_lst(5,:), 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 10, 'Color', 'green', 'Linewidth', 1);hold on;
semilogy(SNR_list, SER_ana_lst(5,:), 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 6, 'Color', 'green', 'Linewidth', 1);

hold off;
grid on
legend("$N = 1$, Trad. ASK (sim.)", "$N = 1$, Trad. ASK (comp.)", ...
    "$N = 2$, Trad. ASK (sim.)", "$N = 2$, Trad. ASK (comp.)", ...
    "$N = 4$, Trad. ASK (sim.)", "$N = 4$, Trad. ASK (comp.)", ...
    "$N = 8$, Trad. ASK (sim.)", "$N = 8$, Trad. ASK (comp.)",  ...
    "$N = 16$, Trad. ASK (sim.)", "$N = 16$, Trad. ASK (comp.)", ...
    "Location","southwest", "Interpreter", "Latex", "FontSize",14);
xlim([0 30]);
ylim([5e-2 1]);
xlabel("average SNR per symbol, $\Gamma_{av}$ (dB)", "FontSize",14, "Interpreter","latex");
ylabel("SEP","FontSize",14, "Interpreter", "latex");
title("$M = "+num2str(M)+"$","FontSize",14, "Interpreter","latex");

function [SER_ana, SER_mc] = runSystem(N, M, SNR)
    mu_1 = 0.8 + 0.5j;
    mu_2 = 0.3 + 0.4j;
    sigma_h = 1;
    num_symbols_tx = 1e5;
    K_1 = (abs(mu_1)^2)/(sigma_h^2);
    K_2 = (abs(mu_2)^2)/(sigma_h^2);
    
    alpha = (N * pi/4) * (laguerreL(1/2, -K_1)) * (laguerreL(1/2, -K_2));
    beta = N * ((1 + K_1) * (1 + K_2) - ((pi^2)/16) * ((laguerreL(1/2, -K_1))^2) * ((laguerreL(1/2, -K_2))^2));

    % Channel Parameters
    h_1_r = sqrt((sigma_h^2)/2).*randn(num_symbols_tx, N) + real(mu_1);
    h_1_i = sqrt((sigma_h^2)/2).*randn(num_symbols_tx, N) + imag(mu_1);
    h_1 = h_1_r + 1j.*h_1_i;

    h_2_r = sqrt((sigma_h^2)/2).*randn(num_symbols_tx, N) + real(mu_2);
    h_2_i = sqrt((sigma_h^2)/2).*randn(num_symbols_tx, N) + imag(mu_2);
    h_2 = h_2_r + 1j.*h_2_i;

    h_T_actual = sum(abs(h_1).*abs(h_2),2);

    mu_f = alpha * (sigma_h^2);
    sigma_f = sqrt(beta * (sigma_h^4));

    h_T_approx = sigma_f*randn(num_symbols_tx,1) + mu_f;

    const_amp = (2/3)*(M-1)*(2*M-1);
    sigma = sqrt(2*(beta + alpha^2)*(sigma_h^4)*const_amp/(10^(SNR/10)));
    n = sigma*randn(num_symbols_tx,1);

    % % Linearising SNR and calculating const terms
    % SNR_lin = 10^(SNR/10);
    % const_amp = SNR_lin/(2*(alpha^2 + beta + 1));
    % A = sqrt((6*SNR_lin)/(2*(alpha^2 + beta + 1)*(M - 1)*(2*M - 1)));
    % 
    % % computing traditional ASK symbols
    % symbol_dict = (0:1:(M-1)).*A;

    % TRADITIONAL ASK
    E_h_T_sq = (sigma_h^4)*((alpha^2) + beta);
    P_T = [0:2:2*(M-1)].^2;
    D_T = [-1:2:(2*M-1)].^2;
    y_T_actual = h_T_actual.* sqrt(P_T) + n;
    y_T_approx = abs(h_T_approx).* sqrt(P_T) + n;
    en_y_T_actual = (y_T_actual.^2)./E_h_T_sq;
    en_y_T_approx = (y_T_approx.^2)./E_h_T_sq;
    
    dl_T = D_T(1:M);
    dr_T = D_T(2:M+1);
    dl_T(1) = -1;
    dr_T(M) = -1;
    
    ser_t_actual = zeros(1,M);
    for l=1:1:M
        for k=1:1:num_symbols_tx
            if ((l == 1)) 
                if((en_y_T_actual(k,l) < dr_T(l)))
                    ser_t_actual(l) = ser_t_actual(l) + 1;
                end
            elseif ((l == M))
                if((en_y_T_actual(k,l) > dl_T(l)))
                    ser_t_actual(l) = ser_t_actual(l) + 1;
                end
            else
                if((en_y_T_actual(k,l) > dl_T(l)) && (en_y_T_actual(k,l) < dr_T(l)))
                    ser_t_actual(l) = ser_t_actual(l) + 1;
                end
            end
        end
    end
    ser_t_actual = 1 - ser_t_actual./num_symbols_tx;
    SER_mc = sum(ser_t_actual)/M;

    ser_t_approx = zeros(1,M);
    for l=1:1:M
        for k=1:1:num_symbols_tx
            if ((l == 1)) 
                if((en_y_T_approx(k,l) < dr_T(l)))
                    ser_t_approx(l) = ser_t_approx(l) + 1;
                end
            elseif ((l == M))
                if((en_y_T_approx(k,l) > dl_T(l)))
                    ser_t_approx(l) = ser_t_approx(l) + 1;
                end
            else
                if((en_y_T_approx(k,l) > dl_T(l)) && (en_y_T_approx(k,l) < dr_T(l)))
                    ser_t_approx(l) = ser_t_approx(l) + 1;
                end
            end
        end
    end
    ser_t_approx = 1 - ser_t_approx./num_symbols_tx;
    SER_ana = sum(ser_t_approx)/M;

end