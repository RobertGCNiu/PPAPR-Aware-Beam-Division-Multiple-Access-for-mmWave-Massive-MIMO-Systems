%%%%%%%%%%----- Performance of Adaptive Channel Estimation of MmWave Channels-----%%%%%%%
%%% Niu Guanchong
%%%Chinese University of HongKong, Shenzhen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
% 
%--------------------------------------------------------------------------
clear;clc;

% ----------------------------- System Parameters -------------------------
%这个是为生成用户和channel model做准备
Num_users=4; %RF chains of users
%RF_Num_users = 4;
%res = Num_users-RF_Num_users;
    Rate_SU_group = [];
Rate_HP_group = [];
Rate_BS_group = [];
Rate_HB_group = [];
Rate_BS_selected_group = [];
 Rate_HB_greedy_selected_group = [];
 Rate_BS__greedy_selected_group = [];
 Rate_HB_cvx_group = [];
 TX_range = [64  81 100 144  196  256  400] ;
%for TX_ant= TX_range; %Number of UPA TX antennas
TX_ant = 64;
TX_ant_w=sqrt(TX_ant); % width
TX_ant_h=sqrt(TX_ant); % hieght 
ind_TX_w=reshape(repmat([0:1:TX_ant_w-1],TX_ant_h,1),1,TX_ant_w*TX_ant_h);
ind_TX_h=repmat([0:1:TX_ant_h-1],1,TX_ant_w);
a_TX = zeros(TX_ant, Num_users);
RX_ant=16; %Number of UPA RX antennas
RX_ant_w=sqrt(RX_ant); % width 
RX_ant_h=sqrt(RX_ant); % hieght
ind_RX_w=reshape(repmat([0:1:RX_ant_w-1],RX_ant_h,1),1,RX_ant_w*RX_ant_h);
ind_RX_h=repmat([0:1:RX_ant_h-1],1,RX_ant_w);

% ----------------------------- Channel Parameters ------------------------
%这里必须是1，不然最优的Frf就不等于a_TX了
Num_paths=1; %Number of channel paths
%%%%------------------%%amplification constraint of RF chian-----
    lamda = 2;
% ----------------------------- Simulation Parameters ---------------------
SNR_dB_range=0:2:20;  % SNR in dB
Rate_SU=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
Rate_LB=zeros(1,length(SNR_dB_range));% Will carry the lower bound values
Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
Rate_HP=zeros(1,length(SNR_dB_range)); % Will carry the rate of the proposed algorithm (with analog 
Rate_BS_selected = zeros(1,length(SNR_dB_range));
Rate_HB_cvx = zeros(1, length(SNR_dB_range));
Rate_HB_cut = zeros(1, length(SNR_dB_range));
% and zero-forcing digital precoding)
Rate_HB=zeros(1,length(SNR_dB_range));
Rate_HB_2 = zeros(1,length(SNR_dB_range));
Rate_WF=zeros(1,length(SNR_dB_range));
Rate_HB_selected=zeros(1, length(SNR_dB_range));
Rate_BS__greedy_selected=zeros(1, length(SNR_dB_range));
Rate_HB_greedy_selected = zeros(1, length(SNR_dB_range));
%%%这就是realization的次数
ITER=500; % Number of iterations
%总过有多少个用户可以供选择    
users_group =  8;%[6, 8, 12, 15, 20, 30];
Fbb_all = zeros(Num_users, Num_users, ITER);
Fbb_cvx_all = zeros(Num_users, Num_users, ITER);
Fbb__cut_all = zeros(Num_users, Num_users, ITER);
Fbb_selected_all = zeros(Num_users, Num_users, ITER);

% --------------- Simulation starts ---------------------------------------


    
  for users = users_group
       
      Rate_SU=zeros(1,length(SNR_dB_range)); % Will carry the single-user MIMO rate (without interference)
Rate_LB=zeros(1,length(SNR_dB_range));% Will carry the lower bound values
Rate_BS=zeros(1,length(SNR_dB_range));% Will carry the rate with analog-only beamsteering
Rate_HP=zeros(1,length(SNR_dB_range)); % Will carry the rate of the proposed algorithm (with analog 
Rate_BS_selected = zeros(1,length(SNR_dB_range));
Rate_HB_cvx = zeros(1, length(SNR_dB_range));
Rate_HB_cut = zeros(1, length(SNR_dB_range));
% and zero-forcing digital precoding)
Rate_HB=zeros(1,length(SNR_dB_range));
Rate_HB_2 = zeros(1,length(SNR_dB_range));
Rate_WF=zeros(1,length(SNR_dB_range));
Rate_HB_selected=zeros(1, length(SNR_dB_range));
Rate_BS__greedy_selected=zeros(1, length(SNR_dB_range));
Rate_HB_greedy_selected = zeros(1, length(SNR_dB_range));
      for iter = 1:ITER
    % Generate user channels and to make the comparasion生成group的所有用户
    [H_user_select,a_TX_user_select,a_RX_user_select]=generate_channels(users,TX_ant_w,TX_ant_h,RX_ant_w,RX_ant_h,Num_paths); 
     
    % H is a 3-dimensional matrix, with Num_users,RX_ant,TX_ant dimensions
    %%%%%%%%%%%%从这里面挑   Frf 让Fbb尽可能控制在小于 lamda 的范围内
  %  [H_fbb_select, a_TX_fbb_select, a_RX_fbb_select, Fbb_select] = select_Fbb(H_user_select, a_TX_user_select, a_RX_user_select,users,Num_users);
    
  %%%%%%%%%%%%%--------greedy selection
   [Frf_greedy_selected,  H_greedy_selected, Wrf_greedy_selected] =  greedyAlgorithm(H_user_select, a_TX_user_select, a_RX_user_select, users, Num_users);
       for u=1:1:Num_users
        Channel_greedy_selected=zeros(RX_ant,TX_ant);
        Channel_greedy_selected(:,:)= H_greedy_selected(u,:,:);
        He_greedy_selected(u,:)=Wrf_greedy_selected(:,u)'*Channel_greedy_selected*Frf_greedy_selected ;    % Effective channels
       end
    
    % Baseband zero-forcing precoding
    Fbb_greedy_selected=He_greedy_selected'*(He_greedy_selected*He_greedy_selected')^(-1);  %inverse(He) 
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb_greedy_selected(:,u)=Fbb_greedy_selected(:,u)/sqrt((Frf_greedy_selected*Fbb_greedy_selected(:,u))'*(Frf_greedy_selected*Fbb_greedy_selected(:,u)));
    end

    %Fbb_cut= zeros(Num_users, Num_users);
%     for u = 1:Num_users
%         for u_c = 1:Num_users
%           if abs(Fbb_selected(u_c,u))>lamda
%               Fbb_selected(u_c,u)=lamda;
%           else
%               Fbb_selected(u_c,u) = Fbb_selected(u_c, u);
%           end
%         end
%     end
    for u=1:1:Num_users %Constraint of RF chain
       if sum(abs(Fbb_greedy_selected(:,u))>lamda)>=1
          Fbb_greedy_selected(:,u) = cvxToOptimize( He_greedy_selected, Fbb_greedy_selected(:,u), lamda );
       else Fbb_greedy_selected(:,u) = Fbb_greedy_selected(:,u);
       end
     end
   
   
   
    %Select the best users to support
    [Frf_selected,  H_selected, Wrf_selected]= users_selection_based_on_AoD(H_user_select, a_TX_user_select, a_RX_user_select, users, Num_users);
    for u=1:1:Num_users
        Channel_selected=zeros(RX_ant,TX_ant);
        Channel_selected(:,:)= H_selected(u,:,:);
        He_selected(u,:)=Wrf_selected(:,u)'*Channel_selected*Frf_selected ;    % Effective channels
    end
    
    % Baseband zero-forcing precoding
    Fbb_selected=He_selected'*(He_selected*He_selected')^(-1);  %inverse(He) 
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb_selected(:,u)=Fbb_selected(:,u)/sqrt((Frf_selected*Fbb_selected(:,u))'*(Frf_selected*Fbb_selected(:,u)));
    end

    %Fbb_cut= zeros(Num_users, Num_users);
%     for u = 1:Num_users
%         for u_c = 1:Num_users
%           if abs(Fbb_selected(u_c,u))>lamda
%               Fbb_selected(u_c,u)=lamda;
%           else
%               Fbb_selected(u_c,u) = Fbb_selected(u_c, u);
%           end
%         end
%     end
    for u=1:1:Num_users %Constraint of RF chain
       if sum(abs(Fbb_selected(:,u))>lamda)>=1
          Fbb_selected(:,u) = cvxToOptimize( He_selected, Fbb_selected(:,u), lamda );
       else Fbb_selected(:,u) = Fbb_selected(:,u);
       end
     end
    
    
    %Previous Work by Robert Heath
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %=======  Add by Robert, here I want to randomply select Num_users
    %users to make zero forcing.
    %%%从这些用户中“随机选择”一些to support
    H=zeros(Num_users,RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h); 
    H_res=zeros(Num_users,RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);
    rank_user = randperm(users);
     for u = 1:1:Num_users
     H(u,:,:) = H_user_select(rank_user(u),:,:);
     a_TX(:, u) = a_TX_user_select(:, rank_user(u));
     a_RX(:, u) = a_RX_user_select(:, rank_user(u));
     end 


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stage 1 of the proposed algorithm (Analog precoding)
    Frf=zeros(TX_ant,Num_users); % BS RF precoders 
    Wrf=zeros(RX_ant,Num_users); % MS RF precoders 
    %%%path=1的时候，最好的Frf和Wrf就等于steering vector
    for u=1:1:Num_users
        Frf(:,u)=a_TX(:,u);
        Wrf(:,u)=a_RX(:,u);
    end      
    
   %%%%%%%%%%%%%%%%%%%%%%
   %%%==half hybrid, half analog only====%%%%
   
   %%%%This part haven't finished. I need to design the algorithm to decide
   %%%%which user applys hybrid beamforming, and which one uses only analog
   %%%%beamforming
%    users_res = users - RF_Num_users;
%    support_user_res = Num_users - RF_Num_users;
%     for u = 1:1:RF_Num_users
%      H_two_way(u,:,:) = H_user_select(rank_user(u),:,:);
%      a_TX_two_way(:, u) = a_TX_user_select(:, rank_user(u));
%      a_RX_two_way(:, u) = a_RX_user_select(:, rank_user(u));
%      Frf_two_way(:, u) = a_TX_two_way(:,u);
%     Wrf_two_way(:,u) = a_RX_two_way(:, u);
%      end 
%    H_de = zeros(users - RF_Num_users, RX_ant_w*RX_ant_h,TX_ant_w*TX_ant_h);
%    u_index = 0;
%    a_TX_de = zeros(TX_ant_w*TX_ant_h, users-RF_Num_users);
%    a_RX_de = zeros(RX_ant_w*RX_ant_h, users-RF_Num_users);
%     for u = 1:1:users           
%          if ~ismember(u, rank_user(1:RF_Num_users))
%              u_index = u_index + 1;
%              H_de(u_index,:,:) = H_user_select(u,:,:);
%              a_TX_de(:, u_index) = a_TX_user_select(:, u);
%              a_RX_de(:, u_index) = a_RX_user_select(:, u);
%          end
%      end
% 
%   [Frf_selected_res,  H_selected_res, Wrf_selected_res]= withSpace(H_de, a_TX_de, a_RX_de, users_res, support_user_res, Frf_two_way);
% [Frf_selected_res,  H_selected_res, Wrf_selected_res]= users_selection_based_on_AoD(H_de, a_TX_de, a_RX_de, (users-RF_Num_users), abs(Num_users-RF_Num_users), Frf_two_way);
% 
% H_res(1:RF_Num_users, :, :) = H_two_way;
% H_res((RF_Num_users+1:Num_users), :, :) = H_selected_res;
%  Frf_res(:, 1:RF_Num_users) =  Frf_two_way;
%  Frf_res(:, (RF_Num_users+1:Num_users)) =  Frf_selected_res(:, (RF_Num_users+1:Num_users)) ;
% Wrf_res(:,1:RF_Num_users) = a_RX_two_way;
% Wrf_res(:,(RF_Num_users+1:Num_users)) = Wrf_selected_res(:, (RF_Num_users+1:Num_users)) ;
%    
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       for u=1:1:RF_Num_users
%         Channel_res=zeros(RX_ant,TX_ant);
%         Channel_res(:,:)= H_res(u,:,:);
%         He_res(u,:)=Wrf_res(:,u)'*Channel_res*Frf_res ;    % Effective channels
%     end
%     
%     Baseband zero-forcing precoding
%     Fbb_res=He_res'*(He_res*He_res')^(-1);   
%     for u=1:1:RF_Num_users % Normalization of the hybrid precoders
%         Fbb_res(:,u)=Fbb_res(:,u)/sqrt((Frf_res*Fbb_res(:,u))'*(Frf_res*Fbb_res(:,u)));
%     end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%
%这个是做zero-forcing来消除不同用户的interference
    % Construct the effective channels
    for u=1:1:Num_users
        Channel=zeros(RX_ant,TX_ant);
        Channel(:,:)= H(u,:,:);
        He(u,:)=Wrf(:,u)'*Channel*Frf ;    % Effective channels
    end
    
    % Baseband zero-forcing precoding
    Fbb=He'*(He*He')^(-1);  %inverse(He) 
    for u=1:1:Num_users % Normalization of the hybrid precoders
        Fbb(:,u)=Fbb(:,u)/sqrt((Frf*Fbb(:,u))'*(Frf*Fbb(:,u)));
    end

    Fbb_cvx = zeros(Num_users, Num_users);
    for u=1:1:Num_users %Constraint of RF chain
       if sum(abs(Fbb(:,u))>lamda)>=1
          Fbb_cvx(:,u) = cvxToOptimize( He, Fbb(:,u), lamda );
       else Fbb_cvx(:,u) = Fbb(:,u);
       end
    end
%     cut Fbb from threshold
Fbb_cut= zeros(Num_users, Num_users);
    for u = 1:Num_users
        for u_c = 1:Num_users
          if abs(Fbb(u_c,u))>lamda
              Fbb_cut(u_c,u)=lamda;
          else
              Fbb_cut(u_c,u) = Fbb(u_c, u);
          end
        end
    end
    
%     %UE zero-forcing decoding
%     Wbb = (He'*He)^(-1)*He';
%     for u=1:1:Num_users % Normalization of the hybrid precoders
%         Wbb(:,u)=Wbb(:,u)/sqrt((Wrf*Wbb(:,u))'*(Wrf*Wbb(:,u)));
%     end
%%%%%%这个可以先不管
     % For the lower bound 
        [~, Ss, ~]=svd(Frf);
        s_min=(min(diag(Ss)))^2;
        s_max=(max(diag(Ss)))^2;
        G_factor=4/(s_max/s_min+s_min/s_max+2);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%下面是开始计算spectral rate
    % Spectral efficiency calculations
    count=0;
    for SNR_dB=SNR_dB_range
        count=count+1;
        SNR=10^(.1*SNR_dB)/Num_users; % SNR value
              Lamda=[];
        for u=1:1:Num_users
            Int_set=[]; % interference index
            
            for i=1:1:Num_users
                if(i~=u)
                    Int_set=[Int_set i]; 
                end
            end
            Channel=zeros(RX_ant,TX_ant);
            Channel(:,:)= H(u,:,:);
            [U_channel, S_channel, ~]=svd(Channel);
            
            Channel_selected=zeros(RX_ant,TX_ant);
            Channel_selected(:,:)= H_selected(u,:,:);
            [U_channel_selected, S_channel_selected, ~]=svd(Channel_selected);
            
            Channel_greedy_selected = zeros(RX_ant, TX_ant);
            Channel_greedy_selected(:,:) = H_greedy_selected(u,:,:);
            
%             Channel_res(:,:) = H_res(u, :, :);
            
            % Single-user rate
            Rate_SU(count)=Rate_SU(count)+log2(1+SNR*(S_channel_selected(1,1))^2)/(Num_users*ITER);

            
            % Analog-only beamforming
            SINR_BS=(SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf(:,Int_set)).^2))+1);
            Rate_BS(count)=Rate_BS(count)+log2(1+SINR_BS)/(Num_users*ITER);
            
            %%%%%%%%%%%%%%%
            %%%===Analog-only by user selection
            
            SINR_BS_selected=(SNR*(abs(Wrf_selected(:,u)'*Channel_selected*Frf_selected(:,u)).^2))/(SNR*sum((abs(Wrf_selected(:,u)'*Channel_selected*Frf_selected(:,Int_set)).^2))+1);
            Rate_BS_selected(count)=Rate_BS_selected(count)+log2(1+SINR_BS_selected)/(Num_users*ITER);
            
            %%%------------Analog-only by greedy algorithm
            SINR_BS_greedy_selected=(SNR*(abs(Wrf_greedy_selected(:,u)'*Channel_greedy_selected*Frf_greedy_selected(:,u)).^2))/(SNR*sum((abs(Wrf_greedy_selected(:,u)'*Channel_greedy_selected*Frf_greedy_selected(:,Int_set)).^2))+1);
            Rate_BS__greedy_selected(count)=Rate_BS__greedy_selected(count)+log2(1+SINR_BS_greedy_selected)/(Num_users*ITER);
            
            
            % Derived lower bound
            Rate_LB(count)=Rate_LB(count)+log2(1+SNR*S_channel(1,1)^2*G_factor)/(Num_users*ITER);
            
            %Proposed Hybrid precoding for 3 RF
            if u ~=Num_users+1        %The first tree users will always be interfered by fourth channel
         SINR_HB = (SNR*(abs(He(u,:)*Fbb(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf*Fbb(:, Int_set)).^2))+1);            
            else               %The forth user is interfered by other three channel
         SINR_HB = (SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf(:,Int_set)).^2))+1);
%          SINR_HB = (SNR*(abs(Wbb(:,u)'*Wrf'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wbb(:,u)'*Wrf'*Channel*Frf(:,Int_set)).^2))+1);
         %  SINR_HB = (SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2));  
            end   %end for Propsed hybrid precoding
            Rate_HB(count)=Rate_HB(count)+log2(1+SINR_HB)/(Num_users*ITER);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  cvx
             SINR_HB_cvx = (SNR*(abs(He(u,:)*Fbb_cvx(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf*Fbb_cvx(:, Int_set)).^2))+1);  
             Rate_HB_cvx(count)=Rate_HB_cvx(count)+log2(1+SINR_HB_cvx)/(Num_users*ITER);
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             SINR_HB_cut = (SNR*(abs(He(u,:)*Fbb_cut(:,u)).^2))/(SNR*sum((abs(Wrf(:,u)'*Channel*Frf*Fbb_cut(:, Int_set)).^2))+1);  
             Rate_HB_cut(count)=Rate_HB_cut(count)+log2(1+SINR_HB_cut)/(Num_users*ITER);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
           %Proposed Hybrid precoding for 2 RF
             SINR_HB_selected = (SNR*(abs(He_selected(u,:)*Fbb_selected(:,u)).^2))/(SNR*sum((abs(Wrf_selected(:,u)'*Channel_selected*Frf_selected*Fbb_selected(:, Int_set)).^2))+1);  
             Rate_HB_selected(count)=Rate_HB_selected(count)+log2(1+SINR_HB_selected)/(Num_users*ITER);
           %  -------------------    greedy   Fbb_greedy_selected
           SINR_HB_greedy_selected = (SNR*(abs(He_greedy_selected(u,:)*Fbb_greedy_selected(:,u)).^2))/(SNR*sum((abs(Wrf_greedy_selected(:,u)'*Channel_greedy_selected*Frf_greedy_selected*Fbb_greedy_selected(:, Int_set)).^2))+1);  
             Rate_HB_greedy_selected(count)=Rate_HB_greedy_selected(count)+log2(1+SINR_HB_greedy_selected)/(Num_users*ITER);
           
           
           % if u<=RF_Num_users      %The first two users will always be interfered by third and fourth channel
            %      SINR_HB_2 = (SNR*(abs(He_res(u,:)*Fbb_res(:,u)).^2))/(SNR*sum((abs(Wrf_res(:,u)'*Channel_res*Frf_res).^2)+1));            
          %  else                         
             %     SINR_HB_2 = (SNR*(abs(Wrf_res(:,u)'*Channel_res*Frf_res(:,u)).^2))/(SNR*sum(abs(Wrf_res(:,u)'*Channel_res*Frf_res(:,Int_set)).^2)+1);
         %SINR_HB_2 = (SNR*(abs(Wbb(:,u)'*Wrf'*Channel*Frf(:,u)).^2))/(SNR*sum((abs(Wbb(:,u)'*Wrf'*Channel*Frf(:,Int_set)).^2))+1);
        %  SINR_HB_2 = (SNR*(abs(Wrf(:,u)'*Channel*Frf(:,u)).^2));  
            %end   %end for Propsed hybrid precoding
%             Rate_HB_2(count)=Rate_HB_2(count)+log2(1+SINR_HB_2)/(Num_users*ITER);         
        end       %end for different users
        
       %water-filling algorithm
        %    Lamda_r = Lamda';
          %  Gamma=Water_Pouring(Lamda_r,SNR,Num_users);
         %   Rate_WF(count) = Rate_WF(count)+log2(det(eye(Num_users)+SNR*diag(Lamda_r)*diag(Gamma)))/(Num_users*ITER);
        
        % Hybrid Precoding
        Rate_HP(count)=Rate_HP(count)+log2(det(eye(Num_users)+SNR*(He*(Fbb*Fbb')*He')))/(Num_users*ITER);

    end% End of SNR loop
%end

Fbb_all(:,:,iter) = Fbb;
Fbb_cut_all(:,:,iter) = Fbb_cut;
Fbb_cvx_all(:,:,iter) = Fbb_cvx;
Fbb_selected_all(:,:,iter) = Fbb_selected;
   end %End the group
Rate_SU_group = [Rate_SU_group Rate_SU];
Rate_HP_group = [Rate_HP_group Rate_HP];
Rate_BS_group = [Rate_BS_group Rate_BS];
Rate_HB_group = [Rate_HB_group Rate_HB];
Rate_HB_cvx_group = [Rate_HB_cvx_group Rate_HB_cvx];
Rate_BS_selected_group = [Rate_BS_selected_group Rate_BS_selected];
 Rate_HB_greedy_selected_group = [ Rate_HB_greedy_selected_group  Rate_HB_greedy_selected];
 Rate_BS__greedy_selected_group = [Rate_BS__greedy_selected_group Rate_BS__greedy_selected];
end % End of ITER loop

hold on
cdf_test(Fbb_all);
cdf_test(Fbb_cut_all);
cdf_test(Fbb_cvx_all);
cdf_test(Fbb_selected_all);
legend('zf','cut','cvx','selected');
%Plotting the spectral efficiencies
     plot(SNR_dB_range,Rate_SU,'-ms','LineWidth',1.5);
%   hold on; plot(SNR_dB_range,Rate_HP,'--s','LineWidth',1.5);
%  %  hold on; plot(SNR_dB_range,Rate_WF,'--s','LineWidth',1.5);
% % %if Num_paths==1
% %   % hold on; plot(SNR_dB_range,4*Rate_LB,'--k','LineWidth',1.5);
    hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
    hold on; plot(SNR_dB_range,Rate_HB,'-ks','LineWidth',1.5);  %%zero forcing
                      plot(SNR_dB_range,Rate_HB_cvx,'-gs','LineWidth',1.5);
                      plot(SNR_dB_range,Rate_HB_cut,'-s','LineWidth',1.5);
                      plot(SNR_dB_range, Rate_HB_selected, '-bs', 'LineWidth', 1.5);
                      plot(SNR_dB_range, Rate_HB_greedy_selected, '--cs', 'LineWidth',1.5);
     hold on; plot(SNR_dB_range,Rate_BS_selected,'--ms','LineWidth',1.5);
     hold on; plot(SNR_dB_range,Rate_BS__greedy_selected,'-ys','LineWidth',1.5);
%    hold on; plot(SNR_dB_range,Rate_HB_2,'-ks','LineWidth',1.5); 
%     legend('Digital Precoding for single Users (No Interference SVD)','Zero-Forcing Precoding N_{RF} = N_s=4','Analog-only Beamsteering','Hybrid Beamforming with MDP N_{RF}=3 N_s=4','Analog only by selection');
%  legend('Digital Precoding for single Users (No Interference SVD)','Zero-Forcing Precoding N_{RF} = N_s=4','Analog-only Beamsteering','HybridBeamforming N_{RF}=3 N_s=4','HybridBeamforming2RF');
 %  legend('Single-user (No Interference)','Proposed Hybrid Precoding','Lower Bound (Theorem 1)','Analog-only Beamsteering','HybridBeamforming3RF','HybridBeamforming2RF');
% else
%     hold on; plot(SNR_dB_range,Rate_BS,'-ro','LineWidth',1.5);
%     legend('Single-user (No Interference)','Proposed Hybrid Precoding','Analog-only Beamsteering');
% end
%legend('su', 'analog only','zf', 'cvx', 'clipping');
legend('su', 'analog only', 'zf', 'cvx', 'clipping','cvx with selected', 'analog selected', 'greedy algorithm');
xlabel('SNR (dB)','FontSize',12);
ylabel('Average of Spectral Efficiency (bps/ Hz)','FontSize',12);
grid;