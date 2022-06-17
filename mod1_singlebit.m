% Title: Simplified 1st-order single-bit incremental modulator
% Author: Bo Wang; Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: Matlab simulation of the output quantization error using different decoders for a single-bit 1st-order modulator, including CoI1, CoI2, sinc2, optimal, and the proposed near-optimal decoder
% Pure differential NTF is used for the modulator model
% Cite as:  B. Wang, M. -K. Law, S. B. Belhaouari and A. Bermak, "Near-Optimal Decoding of Incremental Delta-Sigma ADC Output," IEEE Trans. Circuits Syst. I, vol. 67, no. 11, pp. 3670-3680, Nov. 2020.

clear;
format long;

%% simulation parameter setup
OSR = 1024;             %oversampling ratio of the converter
noise_flag = 0;         %1-turn on thermal noise simulation; 0-turnoff noise
noise_seed = normrnd(0,64e-6,[40000,1]);%thermal noise seed
step_size = 17.9031e-6; %simulation steps, use small, irregular number
start_point = -1;       %minimum input
end_point = 1;          %maximum input
input_scaling = 1;      %input scaling for maximum integrator output control
vinsp = [start_point:step_size:end_point];
N = length(vinsp);
full_scale = vinsp(N)-vinsp(1);
q = zeros(N,OSR);       %stores the quantizer output

%% simulate the modulator with different input
for i = 1:N
    U1 = zeros(1,OSR);
    for j = 1:OSR
        noise_in = 0;%noise_flag*noise_seed(unidrnd(40000),1);
        if (j==1)
            U1(1,j) = input_scaling*(vinsp(i)+noise_in);
            %no feedback in the fist cycle
        else
            U1(1,j) = U1(1,j-1)+input_scaling*(vinsp(i)+noise_in-q(i,j-1));
        end
        if U1(1,j) >= 0
            q(i,j) = 1;
        else
            q(i,j) = -1;
        end
    end
    if i == floor(N/2) %to visualize the simulation progress
        fprintf('---50%%--- finish\n');
    end
end
fprintf('--100%%--- finish\n');

%% output decoding by different filters
vo_coi1 = zeros(N,1);
vo_coi2 = zeros(N,1);
vo_sinc2 = zeros(N,1);
vo_opt = zeros(N,1);
vo_propose = zeros(N,1);
parfor i = 1:N
    temp = q(i,:);
    %coi1 filter
    vo_coi1(i,1) = sum(temp)/OSR;
    %coi2 filter
    vo_coi2(i,1) = sum(temp.*(OSR:-1:1))/sum(1:1:OSR);
    %sinc2 filter
    coefficient = ones(1,floor((OSR+1)/2));
    coefficient_sinc2 = conv(coefficient,coefficient);
    vo_sinc2(i,1) = sum(temp(1:length(coefficient_sinc2)).*coefficient_sinc2)/sum(coefficient_sinc2);
    %optimal
    ub = 1;
    lb = -1;
    sn = 0;
    pub = 1;
    plb = 1;
    for n = 1:OSR
        if temp(n) == -1
            ub = min(ub,sn/n);
            if (ub == sn/n)
                pub = n;
            end
        else
            lb = max(lb,sn/n);
            if (lb == sn/n)
                plb = n;
            end
        end
        sn = sn + temp(n);
    end
    %optimal
    vo_opt(i,1) = (ub+lb)/2;
    %proposed
    if (plb > pub)
        vo_propose(i,1) = lb;
    else
        vo_propose(i,1) = ub;
    end
end

%% plot the output curve
figure
plot(vinsp,vo_coi1,'k')
hold on
plot(vinsp,vo_coi2,'g')
plot(vinsp,vo_sinc2,'b')
plot(vinsp,vo_opt,'c')
plot(vinsp,vo_propose,'r')
legend('CoI1','CoI2','sinc2','opt','propose')

%% calculate the MSE of the errors
qpower_coi1 = mse(vinsp' - vo_coi1);
qpower_coi2 = mse(vinsp' - vo_coi2);
qpower_sinc2 = mse(vinsp' - vo_sinc2);
qpower_opt = mse(vinsp' - vo_opt);
qpower_propose = mse(vinsp' - vo_propose)