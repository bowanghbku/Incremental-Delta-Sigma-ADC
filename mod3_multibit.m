% Title: Simplified 3rd-order multi-bit incremental modulator
% Author: Bo Wang; Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: Matlab simulation setup of the output quantization error using different decoders for the single-ended multi-bit third-order modulator, including CoI3, CoI4, sinc3, sinc4, match, optimal, and the near-optimal decoder
% Pure differential NTF is used for the modulator model
% Cite as:  B. Wang, M. -K. Law, S. B. Belhaouari and A. Bermak, "Near-Optimal Decoding of Incremental Delta-Sigma ADC Output," IEEE Trans. Circuits Syst. I, vol. 67, no. 11, pp. 3670-3680, Nov. 2020.

clear;
format long;

%% simulation parameter setup
OSR = 128;              %oversampling ratio
noise_flag = 0;         %1-turn on noise; 0-turnoff noise
noise_seed = normrnd(0,60e-6,[40000,1]); %thermal noise seed
step_size = 0.7031e-6;  %simulation step
start_point = 1/4;      %minimum stable input, change according to modulator configuration
end_point = 3/4;        %maximum stable input, change according to modulator configuration
vinsp = [start_point:step_size:end_point];
N = length(vinsp);
full_scale = vinsp(N)-vinsp(1);
q = zeros(N,OSR);       %stores the quantizer output
q_level = 17;           %number of internal quantizer levels (>2 and preferred to be 2^n+1 for easier implementation)
LSB = 1/(q_level-1);    %LSB of the midtread quantizer

%% simulate the modulator with different input
U1_max = zeros(1,N); %monitor the integrator output
U1_min = zeros(1,N); %monitor the integrator output
U2_max = zeros(1,N); %monitor the integrator output
U2_min = zeros(1,N); %monitor the integrator output
U3_max = zeros(1,N); %monitor the integrator output
U3_min = zeros(1,N); %monitor the integrator output
for i = 1:N
    U1 = zeros(1,OSR);
    U2 = zeros(1,OSR);
    U3 = zeros(1,OSR);
    Y = zeros(1,OSR);
    for j = 1:OSR
        noise_in = 0;%noise_flag*noise_seed(unidrnd(40000),1);
        if (j==1)
            U1(1,j) = vinsp(i) + noise_in;
            U2(1,j) = 0;
            U3(1,j) = 0;
            q(i,j) = 0;
        else
            if (j==2)
                U1(1,j) = U1(1,j-1) + vinsp(i) + noise_in;
                U2(1,j) = U2(1,j-1) + U1(1,j-1);
                U3(1,j) = U3(1,j-1) + U2(1,j-1);
                q(i,j) = 0;
            else
                U1(1,j) = U1(1,j-1) + vinsp(i) + noise_in - q(i,j-1)*LSB;
                U2(1,j) = U2(1,j-1) + U1(1,j-1) - 3*q(i,j-1)*LSB;
                U3(1,j) = U3(1,j-1) + U2(1,j-1) - 3*q(i,j-1)*LSB;
                q(i,j) = min(max(round(U3(1,j)/LSB),0),q_level-1);%multibit ad
            end
        end
    end
    U1_max(1,i) = max(U1); %maximum integrator output
    U1_min(1,i) = min(U1); %minimum integrator output
    U2_max(1,i) = max(U2);
    U2_min(1,i) = min(U1);
    U3_max(1,i) = max(U3);
    U3_min(1,i) = min(U3);
    if i == floor(N/2)
        fprintf('---50%%--- finish\n'); % to display the simulation progress
    end
end
fprintf('--100%%--- finish\n');

%% check the quantizer stability
if (max(U3_max)>1+LSB/2)
    warning('quantizer out-of-range, please check U3_max for details. For implementable modulator design, please increase the quantizer level or reduce the input range (decrease end_point)');
end
if (min(U3_min)<-LSB/2)
    warning('quantizer out-of-range, please check U3_min for details. For implementable modulator design, please increase the quantizer level or reduce the input range (increase start_point)');
end

%% output decoding by different filters
vo_coi3 = zeros(N,1);
vo_coi4 = zeros(N,1);
vo_sinc3 = zeros(N,1);
vo_sinc4 = zeros(N,1);
vo_match = zeros(N,1);
vo_opt = zeros(N,1);
vo_propose = zeros(N,1);
q(:,1:2) = []; %the first two bits are zero, delete
OSR = OSR-2;
for i = 1:N
    temp = q(i,:);
    %coi3 filter
    coefficient_tmp = conv(1:1:OSR,ones(1,OSR));
    coefficient_coi3 = flip(coefficient_tmp(1:OSR));
    vo_coi3(i,1) = sum(temp.*coefficient_coi3)/sum(coefficient_coi3(1:OSR))/(q_level-1);
    %coi4 filter
    coefficient_tmp = conv(1:1:OSR,1:1:OSR);
    coefficient_coi4 = flip(coefficient_tmp(1:OSR));
    vo_coi4(i,1) = sum(temp.*coefficient_coi4)/sum(coefficient_coi4(1:OSR))/(q_level-1);
    %sinc3 filter
    coefficient = ones(1,floor((OSR+2)/3));
    coefficient_sinc3 = conv(conv(coefficient,coefficient),coefficient);
    vo_sinc3(i,1) =sum(temp(1:length(coefficient_sinc3)).*coefficient_sinc3)/sum(coefficient_sinc3)/(q_level-1);
    %sinc4 filter
    coefficient = ones(1,floor((OSR+3)/4));
    coefficient_sinc4 = conv(conv(coefficient,coefficient),conv(coefficient,coefficient));
    vo_sinc4(i,1) =sum(temp(1:length(coefficient_sinc4)).*coefficient_sinc4)/sum(coefficient_sinc4)/(q_level-1);
    %opt decoding filter
    ub = end_point*(q_level-1);
    lb = start_point*(q_level-1);
    sn = 0;
    rn = 0;
    pub = 1;
    plb = 1;
    for n = 1:OSR
        rn = n*(n+1)*(n+2)/6;
        if (n==1)
            sn = 0;
        else
            if (n==2)
                sn = 3*temp(1,1);
            else
                tmp_co = zeros(1,n-1);
                tmp_co(1,1) = 3;
                for j=2:n-1
                    tmp_co(1,j) = tmp_co(1,j-1) + (3+j-2);
                end
                sn = sum(tmp_co.*temp(1,n-1:-1:1));
            end
        end
        ub = min((sn + min(temp(1,n)+0.5,q_level-1))/rn,ub);
        if(ub == (sn + min(temp(1,n)+0.5,q_level-1))/rn)
            pub = n;
        end
        lb = max((sn + max(temp(1,n)-0.5,0))/rn,lb);
        if(lb == (sn + max(temp(1,n)-0.5,0))/rn)
            plb = n;
        end
        %match filter
        vo_match(i,1) = sn/rn/(q_level-1);
    end
    %optimal
    vo_opt(i,1) = (ub+lb)/2/(q_level-1);
    %proposed
    if (plb > pub)
        vo_propose(i,1) = lb/(q_level-1);
    else
        vo_propose(i,1) = ub/(q_level-1);
    end
end
%gain correction for the match filter
tmp_match = vinsp'-vo_match;
P = polyfit(vinsp', tmp_match,1);
vo_match = vo_match+vinsp'*P(1);

%% plot the output curve
figure
stairs(vinsp,vo_coi3,'k')
hold on
stairs(vinsp,vo_coi4,'g')
stairs(vinsp,vo_sinc3,'b')
stairs(vinsp,vo_sinc4,'m')
stairs(vinsp,vo_match,'y')
stairs(vinsp,vo_opt,'c')
stairs(vinsp,vo_propose,'r')
legend('CoI3','CoI4','sinc3','sinc4','match','opt','propose')

%% calculate the MSE of the conversion errors (the mse results are not meaningful if the modulator is not stable)
if (max(U3_max)>1+LSB/2 || min(U3_min)<-LSB/2)
    warning('the MSE calculation is not meaningful as the quantizer is out-of-range');
end
qpower_coi3 = mse(vinsp'-vo_coi3);
qpower_coi4 = mse(vinsp'-vo_coi4);
qpower_sinc3 = mse(vinsp'-vo_sinc3);
qpower_sinc4 = mse(vinsp'-vo_sinc4);
qpower_match = mse(vinsp'-vo_match);
qpower_opt = mse(vinsp'-vo_opt);
qpower_propose = mse(vinsp'-vo_propose)