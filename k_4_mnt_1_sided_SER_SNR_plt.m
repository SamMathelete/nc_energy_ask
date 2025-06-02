clc;
clear all;
close all;
SNR = -10:1:40;
N_list = [32 64 128 256 512];
M = 4;
SER_list = zeros(5, length(SNR));
SER_T_list = zeros(5, length(SNR));
tic
for j = 1:1:5
    N = N_list(j);
    parfor (i = 1:1:length(SNR), 16)
        [SER_list(j, i), SER_T_list(j, i)] = runSystem(SNR(i), N, M);
    end
    SER_list
end
toc
figure("Name","SER vs SNR plot")
semilogy(SNR, smooth(SER_list(1,:)), '-o', 'Color', '#D95319'); hold on;
semilogy(SNR, smooth(SER_T_list(1,:)), ':o', 'Color', '#D95319');
semilogy(SNR, smooth(SER_list(2,:)), '-^', 'Color', '#A2142F');
semilogy(SNR, smooth(SER_T_list(2,:)), ':^', 'Color', '#A2142F');
semilogy(SNR, smooth(SER_list(3,:)), '-v', 'Color', '#0072BD');
semilogy(SNR, smooth(SER_T_list(3,:)), ':v', 'Color', '#0072BD');
semilogy(SNR, smooth(SER_list(4,:)), '-*', 'Color', '#7E2F8E');
semilogy(SNR, smooth(SER_T_list(4,:)), ':*', 'Color', '#7E2F8E');
semilogy(SNR, smooth(SER_list(5,:)), '-s', 'Color','#77AC30');
semilogy(SNR, smooth(SER_T_list(5,:)), ':s', 'Color','#77AC30');
hold off;
grid on
legend("$N = 32$ (Opt.)", "$N = 32$ (Trad.)", "$N = 64$ (Opt.)", "$N = 64$ (Trad.)", "$N = 128$ (Opt.)", "$N = 128$ (Trad.)", "$N = 256$ (Opt.)", "$N = 256$ (Trad.)", "$N = 512$ (Opt.)", "$N = 512$ (Trad.)", "Location","southwest", "Interpreter", "Latex");
xlim([5 40]);
ylim([1e-3 1e0]);
xlabel("SNR (dB)");
ylabel("SER");

%title("Plot of SER vs SNR when first four moments of channel are perfectly known for 4 ASK")
function [SER, SER_T] = runSystem(SNR, N, M)

% SYSTEM PARAMETERS
mu_1 = 0.1;
mu_2 = 0.2;
sigma_h = 0.6;
num_symbols_tx = 10^5;
K_1 = (mu_1^2)/(sigma_h^2);
K_2 = (mu_2^2)/(sigma_h^2);

alpha = (N * pi/2) * (laguerreL(1/2, -K_1)) * (laguerreL(1/2, -K_2));
beta = N * ((1 + K_1) * (1 + K_2) - ((pi^2)/4) * ((laguerreL(1/2, -K_1))^2) * ((laguerreL(1/2, -K_2))^2));

const_amp = (2/3)*(M-1)*(2*M-1);

sigma = sqrt((2*beta + alpha^2)*(sigma_h^4)*const_amp/(10^(SNR/10)));
mu_f = alpha * (sigma_h^2);
sigma_f = sqrt(beta * (sigma_h^4));

% RIS ASSISTED CHANNEL PARAMETERS
h_T = sigma_f*randn(num_symbols_tx,1) + mu_f;
n = sigma*randn(num_symbols_tx,1);
E_h_T_sq = (sigma_h^4)*((alpha^2) + beta);

% OPTIMISATION OF POWER CODEBOOK
[t_f, P_f, D_f, S_t] = Bisection(M, N, h_T/sqrt(E_h_T_sq), sigma/sqrt(E_h_T_sq));
P_f = P_f .* const_amp;
D_f = D_f .* const_amp;

% FINAL OUTCOME
y_f = abs(h_T).*sqrt(P_f) + n;
en_y_f = (y_f.^2)./E_h_T_sq;
% FINAL CONSTELLATION
r = (E_h_T_sq.*P_f + sigma^2)./E_h_T_sq;
dl_f = r - D_f(:,1)';
dr_f = r + D_f(:,2)';
dr_f(1) = dl_f(2);
dl_f(1) = -1;
dr_f(M) = -1;

% TRADITIONAL ASK
P_T = [0:2:2*(M-1)].^2;
D_T = [-1:2:(2*M-1)].^2;
y_T = abs(h_T).* sqrt(P_T) + n;
en_y_T = (y_T.^2)./E_h_T_sq;

dl_T = D_T(1:M);
dr_T = D_T(2:M+1);
dl_T(1) = -1;
dr_T(M) = -1;

ser_t = zeros(1,M);
for l=1:1:M
    for k=1:1:num_symbols_tx
        if ((l == 1)) 
            if((en_y_T(k,l) < dr_T(l)))
                ser_t(l) = ser_t(l) + 1;
            end
        elseif ((l == M))
            if((en_y_T(k,l) > dl_T(l)))
                ser_t(l) = ser_t(l) + 1;
            end
        else
            if((en_y_T(k,l) > dl_T(l)) && (en_y_T(k,l) < dr_T(l)))
                ser_t(l) = ser_t(l) + 1;
            end
        end
    end
end
ser_t = 1 - ser_t./num_symbols_tx;
SER_T = sum(ser_t)/M;
% PLOTTING THE CONSTELLATION
% figure("Name","Constellation");
% hold on
% plot(r, zeros(1,M), "*", "color", "b")
% plot(dl_f, zeros(1,M), "o", "color","red")
% plot(dr_f, zeros(1,M), "x", "color","yellow")
% title("Optimal constellation points")
% xlabel("(|y|^2)/(\sigma_h^4(\beta + \alpha^2))")
% xlim([-0.5 ceil(P_f(M))+1])
% hold off
% SER CALCULATION
ser = zeros(1,M);
for l=1:1:M
    for k=1:1:num_symbols_tx
        if ((l == 1)) 
            if((en_y_f(k,l) < dr_f(l)))
                ser(l) = ser(l) + 1;
            end
        elseif ((l == M))
            if((en_y_f(k,l) > dl_f(l)))
                ser(l) = ser(l) + 1;
            end
        else
            if((en_y_f(k,l) > dl_f(l)) && (en_y_f(k,l) < dr_f(l)))
                ser(l) = ser(l) + 1;
            end
        end
    end
end
ser = 1 - ser./num_symbols_tx;
SER = sum(ser)/M;
end
function [t_f, P_f, D_f, S_t] = Bisection(M, N, h_T, sigma)
    t = -1;
    t_l = 0;
    t_u = N/10;
    eps_t = 1e-4;
    eps_s = 1e-4;
    S_t_calc = 1e9;
    while (abs(t_u - t_l) > eps_t && abs(S_t_calc - 1) > eps_s)
        t = (t_l + t_u) / 2;
        [P_calc, D_calc, S_t_calc] = ConstellationDesign(M, h_T, sigma, t);
        if (S_t_calc < 1)
            t_l = t;
        else
            t_u = t;
        end
    end
    t_f = t;
    P_f = P_calc;
    D_f = D_calc;
    S_t = S_t_calc;
end
function [P_calc, D_calc, S_t_calc] = ConstellationDesign(M, h_T, sigma, t)
    alpha_1 = mean(h_T.^4) - 1;
    alpha_2 = 2*(sigma^2);
    alpha_3 = sigma^4;
    s = @(p_k) (alpha_1*(p_k^2) + alpha_2*(p_k) + alpha_3);
    P = zeros(1,M);
    D = zeros(M,2);
    d_l = zeros(M,1);
    d_r = ones(M,1).*sqrt(2*t*s(0));
    d_l(1) = 1e5;
    d_r(M) = 1e5;
    
    for k = 2:1:M
        Pkl = P(k-1);
        Pku = M;
        P(k) = (Pkl + Pku)/2;
        err = ((P(k) - P(k - 1))^2)/(2*(sqrt(s(P(k))) + sqrt(s(P(k-1))))^2) - t;
        while (abs(Pkl - Pku) > 1e-8 && abs(err) > 1e-8)
            if (err > 0)
                    Pku = P(k);
                else
                    Pkl = P(k);
            end
            P(k) = (Pkl + Pku)/2;
            err = ((P(k) - P(k - 1))^2)/(2*(abs(sqrt(s(P(k)))) + abs(sqrt(s(P(k-1)))))^2) - t;
        end
        if (k < M)
            d_r(k) = sqrt(2*t*s(P(k)));
        end
        d_l(k) = sqrt(2*t*s(P(k)));
        D(k,1) = d_l(k);
        D(k,2) = d_r(k);
    end
    P_calc = P;
    D_calc = D;
    D_calc(1,1) = -1;
    D_calc(M,2) = -1;
    S_t_calc = sum(P)/M;
end