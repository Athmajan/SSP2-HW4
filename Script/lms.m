clear all;
close all;
clc;
rng(42);
T=1000;    % number of training samples


N=3e3;    % number of data samples
SNRdB=25;     % SNR in dB value
L=20; % length for smoothing(L+1)
ChL=5;  % length of the channel(ChL+1)
mu = 0.01;  % Step size of the LMS algorithm
%SERlist = zeros(size(0.03:0.001:0.05));
ath = 1;
%for mu=0.03:0.001:0.05
    



%% Training phase
% QPSK symbols to transmit
d=round(rand(1,T))*2-1;  
d=d+1i*(round(rand(1,T))*2-1);

% Channel 
Ch=randn(1,ceil(ChL/2))+1i*randn(1,ceil(ChL/2));
Ch=[Ch Ch(end-1:-1:1)];
Ch=Ch/norm(Ch);
 

% signal filtered by channel
x= filter(Ch,1,d); 


% Noise
v=randn(1,T);  
v=v/norm(v)*10^(-SNRdB/20)*norm(x);

% Received signal
u=x+v;                          

%%%%%%%%%%%%%% WRITE LMS ALGORITHM FOR EQUALIZATION HERE %%%%%%%%%%%%%%%%%%
EqD= round((L+ChL)/2); %delay of input desired response

% Initialize filter coefficients
w_hat = zeros(L, 1);

% Initialize error signal
e = zeros(T, 1);
u  = [zeros(1,EqD) u zeros(1,L) ];

for n=1:T
    un= transpose(u(n+L-1:-1:n));
    e(n) = d(n)-w_hat'*un;

    w_hat = w_hat + mu*un*e(n)';
end

tapInputs = un;

%% Transmission phase
d=round(rand(1,N))*2-1;  
d=d+1i*(round(rand(1,N))*2-1);

% signal filtered by FIR delay line channel
x= filter(Ch,1,d);

% Noise
v=randn(1,N);  
v=v/norm(v)*10^(-SNRdB/20)*norm(x);

% Received signal
u=x+v;  


%%%% Add proper padding (zeros) to u %%%%%%%%



Un=zeros(L,N);% Un is the matrix such that n-th column 
               % corresponds to vector:
               % u(n) = [u(n) u(n-1) ... u(n-L+1)] 
for n=1:N-L
    Un(:,n) = flip(transpose(u(n:n+L-1))); 
end

sb=w_hat'*Un;   % recieved symbol estimation where w_hat is LMS filter in training

%SER(decision part)
sb1=sb/norm(w_hat);  % normalize the output
sb1=sign(real(sb1))+1i*sign(imag(sb1));  %symbol detection
sb2=sb1-d(1:length(sb1));  % error detection
SER=length(find(sb2~=0))/length(sb2); %  SER calculation
disp(SER);
disp(ath);
SERlist(ath) = SER;
ath = ath +1 ;
%end

% plot of transmitted symbols
    subplot(2,2,1), 
    plot(d,'*');   
    grid,title('Input symbols');  xlabel('real part'),ylabel('imaginary part')
    axis([-2 2 -2 2])
    
% plot of received symbols
    subplot(2,2,2),
    plot(u,'o');
    grid, title('Received samples');  xlabel('real part'), ylabel('imaginary part')

% plots of the equalized symbols    
    subplot(2,2,3),
    plot(sb,'o');   
    grid, title('Equalized symbols'), xlabel('real part'), ylabel('imaginary part')

% convergence of LMS algorithm
    subplot(2,2,4),
    plot(abs(e));   % e is estimation error in LMS algorithm
    grid, title('Convergence'), xlabel('n'), ylabel('error signal')
   