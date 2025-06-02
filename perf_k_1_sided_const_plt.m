clc;
clear all;
close all;

SNR = [10 40];
N_list = [128 512];
M = [4 8];
cnt = 1;
figure("Name", "Constellation");
for k = 1:1:2
    for j = 1:1:2
        for i = 1:1:2
            subplot(4, 2, cnt);
            [A, A_T] = runSystem(SNR(i), N_list(j), M(k));
            hold on
            plot(A_T, zeros(1, M(k)), "o", "Color","red");
            plot(A, zeros(1, M(k)), "*", "Color", "blue");
            hold off
            grid on
            title("$M = $ "+num2str(M(k))+", $N = $ "+num2str(N_list(j))+", SNR = "+num2str(SNR(i))+"dB", "Interpreter","latex","FontSize",13);
            cnt = cnt + 1
        end
    end
end
legend("Trad.", "Opt.");
sgtitle("Constellation Plot for One-Sided ASK");

function [A, A_T] = runSystem(SNR, N, M)

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

sigma = sqrt(2*(beta + alpha^2)*(sigma_h^4)*const_amp/(10^(SNR/10)));

mu_f = alpha * (sigma_h^2);
sigma_f = sqrt(beta * (sigma_h^4));

% RIS ASSISTED CHANNEL PARAMETERS
h_T = sigma_f*randn(num_symbols_tx,1) + mu_f;
n = sigma*randn(num_symbols_tx,1);
E_h_T_sq = (sigma_h^4)*((alpha^2) + beta);

% TRANSMITTER SIDE
P_i = ([0:2:2*(M-1)].^2)./const_amp;

% RECEIVER SIDE
y = abs(h_T).*sqrt(P_i) + n;
en_y = (abs(y).^2)./E_h_T_sq;
E_en_y = mean(en_y);

% OPTIMISATION OF POWER CODEBOOK
u = en_y - E_en_y;
Mgf = cell(1,M);
for i=1:1:M
    Mgf{i} = @(theta) mean(exp(theta .* u(:,i)'));
end
Il = cell(1, M);
for i=1:1:M
    Il{i} = @(d) fminbnd(@(theta) -(theta*d - log(Mgf{i}(-theta))), 0, 100);
end
Ir = cell(1, M);
for i=1:1:M
    Ir{i} = @(d) fminbnd(@(theta) -(theta*d - log(Mgf{i}(theta))), 0, 100);
end
[t_f, P_f, D_f, S_t] = Bisection(M, Il, Ir, Mgf);
P_f = P_f .* const_amp;
D_f = D_f .* const_amp;

% FINAL OUTCOME
y_f = abs(h_T).*sqrt(P_f) + n;
en_y_f = (y_f.^2)./E_h_T_sq;

A = sqrt(P_f);

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

A_T = sqrt(P_T);

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

function [t_f, P_f, D_f, S_t] = Bisection(M, Il, Ir, Mgf)
    t = -1;
    t_l = 0;
    t_u = 1e1;
    eps_t = 1e-2;
    eps_s = 1e-2;
    S_t_calc = 1e9;
    while (abs(t_u - t_l) > eps_t && abs(S_t_calc - 1) > eps_s)
        t = (t_l + t_u) / 2;
        [P_calc, D_calc, S_t_calc] = ConstellationDesign(M, t, Il, Ir, Mgf);
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

function [P_calc, D_calc, S_t_calc] = ConstellationDesign(M, t, Il, Ir, Mgf)
    P = zeros(1,M);
    D = zeros(M,2);
    d_l = zeros(M,1);
    d_r = zeros(M,1);
    d_l(1) = 1e5;
    d_r(M) = 1e5;
    for k=1:1:M
        if (k == 1)
            P(1) = 0;
        else
            Pkl = P(k-1) + d_r(k-1);
            Pku = M;
            P(k) = (Pkl + Pku)/2;
            klm = ((P(k) - P(k-1) - d_r(k-1))*Il{k}(P(k) - P(k-1) - d_r(k-1)) - log(Mgf{k}(-Il{k}(P(k) - P(k-1) - d_r(k-1)))));
            err = klm - t;
            while (abs(Pkl - Pku) > 1e-5 && abs(err) > 1e-2)
                if (err > 0)
                    Pku = P(k);
                else
                    Pkl = P(k);
                end
                P(k) = (Pkl + Pku)/2;
                klm = ((P(k) - P(k-1) - d_r(k-1))*Il{k}(P(k) - P(k-1) - d_r(k-1)) - log(Mgf{k}(-Il{k}(P(k) - P(k-1) - d_r(k-1)))));
                err = klm - t;
            end
            d_l(k) = P(k) - P(k-1) - d_r(k-1);
        end
        if (k < M)
            dru = 1;
            drl = 0;
            d_r(k) = (dru + drl)/2;
            kl = (d_r(k)*Ir{k}(d_r(k)) - log(Mgf{k}(Ir{k}(d_r(k)))));
            err = t - kl;
            while(abs(dru - drl) > 1e-8 && abs(err) > 1e-2)
                % t
                % d_r(k)
                % (d_r(k)*Ir{k}(d_r(k)) - log(Mgf{k}(Ir{k}(d_r(k)))))
                % t - (d_r(k)*Ir{k}(d_r(k)) - log(Mgf{k}(Ir{k}(d_r(k)))))
                if (err > 0)
                    drl = d_r(k);
                else
                    dru = d_r(k);
                end
                d_r(k) = (dru + drl)/2;
                kl = (d_r(k)*Ir{k}(d_r(k)) - log(Mgf{k}(Ir{k}(d_r(k)))));
                err = t - kl;
            end
        end
        % k
        % errf = err
        % t - (d_r(k)*Ir{k}(d_r(k)) - log(Mgf{k}(Ir{k}(d_r(k)))))
        D(k,1) = d_l(k);
        D(k,2) = d_r(k);
    end
    P_calc = P;
    D_calc = D;
    D_calc(1,1) = -1;
    D_calc(M,2) = -1;
    S_t_calc = sum(P)/M;
end