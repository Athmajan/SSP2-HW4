clear all;
close all;
clc;

T = 1e3;    % number of training samples
N = 3e3;    % number of data samples
SNRdB = 25;     % SNR in dB value
L = 20; % length for smoothing(L+1)
ChL = 5;  % length of the channel(ChL+1)
mu = 0.01;  % Step size of the LMS algorithm

%% Training phase
% QPSK symbols to transmit
d = round(rand(1, T)) * 2 - 1;  
d = d + 1i * (round(rand(1, T)) * 2 - 1);

% Channel 
Ch = randn(1, ceil(ChL/2)) + 1i * randn(1, ceil(ChL/2)));
Ch = [Ch Ch(end-1:-1:1)];
Ch = Ch / norm(Ch);

% signal filtered by channel
x = filter(Ch, 1, d); 

% Noise
v = randn(1, T);  
v = v / norm(v) * 10^(-SNRdB/20) * norm(x);

% Received signal
u = x + v;                          

%%%%%%%%%%%%%% LMS ALGORITHM FOR EQUALIZATION %%%%%%%%%%%%%%%%%%
EqD = round((L + ChL) / 2); %delay of input desired response

% Initialize filter coefficients
w_hat = zeros(L, 1);

% Initialize error signal
e = zeros(T, 1);

% Apply LMS algorithm
for t = 1:T
    % Extract the desired response
    desired = d(t);

    % Extract input signal vector
    input_signal = flip(u(t:t+L-1));
    
    % Compute filter output
    output = w_hat' * input_signal;
    
    % Compute error
    error = desired - output;
    
    % Update filter coefficients
    w_hat = w_hat + mu * conj(error) * input_signal;
    
    % Store error
    e(t) = error;
end

%% Transmission phase
% Generate QPSK symbols for transmission
d = round(rand(1, N)) * 2 - 1;  
d = d + 1i * (round(rand(1, N)) * 2 - 1);

% Filter symbols with FIR delay line channel
x = filter(Ch, 1, d);

% Generate noise
v = randn(1, N);  
v = v / norm(v) * 10^(-SNRdB/20) * norm(x);

% Received signal
u = x + v;  

% Add proper padding (zeros) to u
u = [zeros(1, L), u]; % Add L zeros to the beginning to align with the filter length

% Initialize matrix to store delayed versions of u
Un = zeros(L, N);

% Apply LMS filter to equalize the received signal
for n = 1:N
    % Construct the input vector for the filter
    input_signal = flip(u(n:n+L-1));
    
    % Apply the filter
    sb(n) = conj(w_hat') * input_signal;
end

% Normalize received symbol estimation
sb = sb / norm(w_hat);

% Decision part (symbol detection)
sb = sign(real(sb)) + 1i * sign(imag(sb));

% Error detection
sb_error = sb - d(1:length(sb));

% Symbol Error Rate (SER) calculation
SER = sum(sb_error ~= 0) / length(sb_error);
disp(['Symbol Error Rate (SER): ', num2str(SER)]);

% Plot of transmitted symbols
subplot(2,2,1);
plot(d, '*');   
grid on;
title('Input symbols');
xlabel('Real part');
ylabel('Imaginary part');
axis([-2 2 -2 2]);

% Plot of received symbols
subplot(2,2,2);
plot(u, 'o');
grid on;
title('Received samples');
xlabel('Real part');
ylabel('Imaginary part');

% Plot of the equalized symbols    
subplot(2,2,3);
plot(sb, 'o');   
grid on;
title('Equalized symbols');
xlabel('Real part');
ylabel('Imaginary part');

% Convergence of LMS algorithm
subplot(2,2,4);
plot(abs(e));   % e is estimation error in LMS algorithm
grid on;
title('Convergence');
xlabel('n');
ylabel('Error signal');
