% Title: Simplified 2nd-order single-bit incremental modulator
% Author: Bo Wang; Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: Matlab simulation of the output quantization error using different decoders for the single-bit 2nd-order CIFB/CIFF modulator, including CoI2, CoI3, sinc2, sinc3, match, optimal, and the near-optimal decoder
% Pure differential NTF is used for the modulator model
% Cite as:  B. Wang, M. -K. Law, S. B. Belhaouari and A. Bermak, "Near-Optimal Decoding of Incremental Delta-Sigma ADC Output," IEEE Trans. Circuits Syst. I, vol. 67, no. 11, pp. 3670-3680, Nov. 2020.

clear;
format long;

%% simulation parameter setup
OSR = 400;                  %oversampling ratio of the converter
noise_flag = 0;             %1-turn on thermal noise simulation; 0-turnoff noise
noise_seed = normrnd(0,64e-6,[40000,1]); %thermal noise seed
step_size = 7.9031e-6;      %simulate the input with small steps when large OSR is used
start_point = -3/4;         %minimum input, usually -3/4 or -2/3 for 2nd-order SDM
end_point = 3/4;            %maximum input, usually 3/4 or 2/3 for 2nd-order SDM
vinsp = [start_point:step_size:end_point];
N = length(vinsp);
full_scale = vinsp(N)-vinsp(1);
q = zeros(N,OSR);           %stores the quantizer output

%% modulator parameters
signal_scaling = 8;         %signal is scaled by 8 to limit the integrator outputs, the reader can change if other scaling will be used
a1 = 1/signal_scaling;
a2 = 0;
a3 = 0;
b1 = 1/signal_scaling;
b2 = 2/signal_scaling;
c1 = 0;
% use the parameters below for ciff topology
% a1 = 1/signal_scaling;
% a2 = 0;
% a3 = 1/signal_scaling;
% b1 = 1/signal_scaling;
% b2 = 0;
% c1 = 2;

%% simulate the modulator with different input
U1_max = zeros(1,N); %monitor the integrator output
U1_min = zeros(1,N); %monitor the integrator output
U2_max = zeros(1,N); %monitor the integrator output
U2_min = zeros(1,N); %monitor the integrator output
Y_max = zeros(1,N); %monitor the quantizer input
Y_min = zeros(1,N); %monitor the quantizer input
for i = 1:N
    U1 = zeros(1,OSR);
    U2 = zeros(1,OSR);
    Y = zeros(1,OSR);
    for j = 1:OSR
        noise_in = 0;%noise_flag*noise_seed(unidrnd(40000),1);
        if (j==1)
            U1(1,j) = vinsp(i)*a1+ a1*noise_in;
            U2(1,j) = 0;
            Y(1,j) = vinsp(i)*a3+c1*U1(1,j)+U2(1,j);
            q(i,j) = 0; %for classical topology
            %if Y(1,j) >= 0 %for ciff topology
            % q(i,j) = 1; %for ciff topology
            %else %for ciff topology
            % q(i,j) = -1; %for ciff topology
            %end %for ciff topology
        else
            U1(1,j) = U1(1,j-1)+a1*vinsp(i)-b1*q(i,j-1) + a1*noise_in;
            U2(1,j) = U2(1,j-1)+U1(1,j-1)+a2*vinsp(i)-b2*q(i,j-1)+a2*noise_in;
            Y(1,j) = a3*vinsp(i)+c1*U1(1,j)+U2(1,j)+a3*noise_in;
            if Y(1,j) >= 0
                q(i,j) = 1;
            else
                q(i,j) = -1;
            end
        end
    end
    U1_max(1,i) = max(U1); %maximum integrator output
    U1_min(1,i) = min(U1); %minimum integrator output
    U2_max(1,i) = max(U2);
    U2_min(1,i) = min(U2);
    Y_max(1,i) = max(Y);
    Y_min(1,i) = min(Y);
    if i == floor(N/2)
        fprintf('---50%%--- finish\n'); %to display the simulation progress
    end
end
fprintf('--100%%--- finish\n');

%% check the quantizer stability
if (max(Y_max)>1.5) %quantizer full scale is 1 + LSB/2 where LSB=1 for 1bit quantizer
    warning('quantizer might be out-of-range, please check Y_max for details. For implementable modulator design, please reduce the input range (decrease end_point)');
end
if (min(Y_min)<-1.5)
    warning('quantizer might be out-of-range, please check Y_min for details. For implementable modulator design, please reduce the input range (increase start_point)');
end

%% output decoding by different filters
%modulator parameters without signal scaling
a1 = 1;
a2 = 0;
a3 = 0;
b1 = 1;
b2 = 2;
c1 = 0;
% use the parameters below for ciff topology
% a1 = 1;
% a2 = 0;
% a3 = 1;
% b1 = 1;
% b2 = 0;
% c1 = 2;
vo_coi2 = zeros(N,1);
vo_coi3 = zeros(N,1);
vo_sinc2 = zeros(N,1);
vo_sinc3 = zeros(N,1);
vo_match = zeros(N,1);
vo_opt = zeros(N,1);
vo_propose = zeros(N,1);
q_opt = q;
q_linear = q;
q_linear(:,1)=[];       %delete the first column with all zeros, comment for ciff
OSR_liner = OSR - 1;    %data length reduces by 1, comment for ciff
%OSR_liner = OSR; %uncomment for ciff
parfor i = 1:N
    temp = q_linear(i,:);
    % coi2 filter
    vo_coi2(i,1) = sum(temp.*(OSR_liner:-1:1))/sum(1:1:OSR_liner); %for classical topology
    %vo_coi2(i,1) = sum(temp.*(OSR_liner:-1:1))/sum(1:1:OSR_liner+1); %for ciff topology
    % coi3 filter
    coefficient_tmp = conv(1:1:OSR_liner,ones(1,OSR_liner));
    coefficient_coi3 = flip(coefficient_tmp(1:OSR_liner));
    vo_coi3(i,1) =     sum(temp.*coefficient_coi3)/sum(coefficient_tmp(1:OSR_liner)); %for classical topology
    %temp3 = temp(1,1:end-3); %for ciff topology
    %temp3 = [0,0,0,temp3]; %for ciff topology
    %vo_coi3(i,1) = sum(temp3.*coefficient_coi3)/sum(coefficient_tmp(1:OSR_liner-2));
    %for ciff topology
    % sinc2 filter
    coefficient = ones(1,floor((OSR_liner+1)/2));
    coefficient_sinc2 = conv(coefficient,coefficient);
    vo_sinc2(i,1) =     sum(temp(1:length(coefficient_sinc2)).*coefficient_sinc2)/sum(coefficient_sinc2);
    % sinc3 filter
    coefficient = ones(1,floor((OSR_liner+2)/3));
    coefficient_sinc3 = conv(conv(coefficient,coefficient),coefficient);
    vo_sinc3(i,1) =    sum(temp(1:length(coefficient_sinc3)).*coefficient_sinc3)/sum(coefficient_sinc3);
    %opt decoding filter
    temp = q_opt(i,:);
    ub = end_point;
    lb = start_point;
    sn = 0;
    rn = a3;
    pub = 1;
    plb = 1;
    tmp1 = 0;
    for n = 1:OSR
        rn = rn + a1*(n-1)+(a1*c1+a2);
        if temp(1,n) == -1
            ub = min(ub,(sn/rn));
            if (ub==sn/rn)
                pub = n;
            end
        else
            lb = max(lb,(sn/rn));
            if (lb==sn/rn)
                plb = n;
            end
        end
        vo_match(i,1) = sn/rn; %match filter
        sn = sn + (b1*c1+b2)*temp(1,n) + b1*tmp1;
        tmp1 = tmp1 + temp(1,n);
    end
    %optimal
    vo_opt(i,1) = (ub+lb)/2;
    % proposed
    if (plb > pub)
        vo_propose(i,1) = lb;
    else
        vo_propose(i,1) = ub;
    end
end

%% plot the output curve
figure
stairs(vinsp,vo_coi2,'k')
hold on
stairs(vinsp,vo_coi3,'g')
stairs(vinsp,vo_sinc2,'b')
stairs(vinsp,vo_sinc3,'m')
stairs(vinsp,vo_match,'y')
stairs(vinsp,vo_opt,'c')
stairs(vinsp,vo_propose,'r')
legend('CoI2','CoI3','sinc2','sinc3','match','opt','propose')

%% calculate the MSE of the conversion errors (change the input range to make the calculation below meaningful)
if (max(Y_max)>1.5 || min(Y_min)<-1.5)
    warning('the MSE calculation is not meaningful as the quantizer is out-of-range');
end
qpower_coi2 = mse(vinsp'-vo_coi2);
qpower_coi3 = mse(vinsp'-vo_coi3);
qpower_sinc2 = mse(vinsp'-vo_sinc2);
qpower_sinc3 = mse(vinsp'-vo_sinc3);
qpower_match = mse(vinsp'-vo_match);
qpower_opt = mse(vinsp'-vo_opt);
qpower_propose = mse(vinsp'-vo_propose)