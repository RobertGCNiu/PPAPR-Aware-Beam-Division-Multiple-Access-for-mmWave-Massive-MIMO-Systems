%
% Matlab program to reLow Complexity Hybrid Beamforming Based on Orthogonal Constraint and Phase Extraction
% By Wenjuan Pu, Xiaohui Li, Yingchao Lin and Ruiyang Yuan
%
% Written By Niu Guanchong
% Last update: Jun 8, 2017
%

clear all;clc;
%%
%
% System Parameters
%
Num_users=1; % Number of users

% Transmitter 
TX_ant=64; %Number of UPA TX antennas
TX_ant_w=sqrt(TX_ant); % width
TX_ant_h=sqrt(TX_ant); % hieght 
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
N_RF = 2; %Number of RF chians
D_lamda = zeros(N_RF, N_RF);

% Receiver
RX_ant=16; %Number of UPA RX antennas
RX_ant_w=sqrt(RX_ant); % width 
RX_ant_h=sqrt(RX_ant); % hieght

ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);

% Channel model
Num_paths = 2;

% Simulation parameters
ITER = 100;
SNR_dB_range=-20:5:10;  % SNR in dB
Rate_HP=zeros(1,length(SNR_dB_range));
Rate_BS=zeros(1,length(SNR_dB_range));
Rate_optimal=zeros(1,length(SNR_dB_range));

%%
for iter = 1: ITER
    Frf = zeros(TX_ant, N_RF);
    
    % generate random channel matrix
    [H,a_TX,a_RX]=generate_channels(Num_users,TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths);

    %convert 3D Channel Model H to 2D H. SVD can only be used in 2D case.
    H_2D=squeeze(H(1,:,:));
    % H_2D = zeros(RX_ant,TX_ant);
    %     for rx = 1:RX_ant
    %         for tx = 1:TX_ant
    %             H_2D(rx,tx) = H(1,rx,tx);
    %         end
    %     end
    
    [Us, Ss, Vs]=svd(H_2D);
    Frf_Opt_tilde=Vs(:,1:N_RF);
    % Algorithm 1: Line 1-3 For-loop
    for jj = 1:N_RF
            Frf(:, jj) = exp(1j*angle(Frf_Opt_tilde(:, jj)));      
            D_lamda(jj, jj) = Frf_Opt_tilde(:, jj)'*Frf(:,jj)/TX_ant;
    end
    % Algorithm 1: Line 4
    Frf_tilde = Frf * D_lamda;
    % Algorithm 1: Line 4
    Heff = H_2D * Frf_tilde;
    Q = Frf' * Frf;
    [~,Uee,Vee] = svd(Heff*Q^(-1/2));
    for rf = 1:N_RF
        Ue(:,rf) = Vee(:, rf);
    end
    % Algorithm 1: Line 5
    Fbb_tilde = Q^(-1/2)*Ue;% F_tilde_BB = Q^(-1/2) * Ue. 
    % Algorithm 1: Line 6
    Fbb = D_lamda*Fbb_tilde;
    % Algorithm 1: Line 7
    Fbb = sqrt(N_RF/(norm(Frf*Fbb,'fro'))^2)*Fbb;

    for ii=1:length(SNR_dB_range)
        SNR=10^(.1*SNR_dB_range(ii))/Num_users; % SNR value    

        % Channel capacity with TX CSI
%         Rate_optimal(ii)=Rate_optimal(ii)+log2(det(eye(RX_ant)+SNR*(H_2D*(Vs(:,1:N_RF)*Vs(:,1:N_RF)')*H_2D')))/(Num_users*ITER);           
        
        % Channel capacity with full CSI at both TX and RX without power
        % allocation
        Rate_optimal(ii)=Rate_optimal(ii)+log2(det(eye(RX_ant)+SNR*(Ss*Ss')))/(Num_users*ITER);
    
        % Analog only with TX CSI
     %    Rate_BS(ii)=Rate_BS(ii)+log2(det(eye(RX_ant)+SNR*(H_2D*(Frf/TX_ant*Frf')*H_2D')))/(Num_users*ITER);
        Rate_BS(ii)=Rate_BS(ii)+log2(det(eye(RX_ant)+SNR*(H_2D*(Frf*Frf'/(norm(Frf*Frf','fro')^2))*H_2D')))/(Num_users*ITER);

        % Hybrid Precoding with TX CSI
        Rate_HP(ii)=Rate_HP(ii)+log2(det(eye(RX_ant)+SNR*(H_2D*Frf_tilde*(Fbb_tilde*Fbb_tilde')*Frf_tilde'*H_2D')))/(Num_users*ITER); 
    end % End of SNR loop
 end%End of iteration
%%
% Plot results
figure(1);
plot(SNR_dB_range,abs(Rate_optimal),'-s','LineWidth',1.5);
hold on; 
plot(SNR_dB_range,abs(Rate_BS),'-kd','LineWidth',1.5);
plot(SNR_dB_range,abs(Rate_HP),'-r*','LineWidth',1.5);
hold off; grid;
legend('Optimal rate','Analog only','Hybrid Precoding','Location','NorthWest');
 
 
 
 
 
 