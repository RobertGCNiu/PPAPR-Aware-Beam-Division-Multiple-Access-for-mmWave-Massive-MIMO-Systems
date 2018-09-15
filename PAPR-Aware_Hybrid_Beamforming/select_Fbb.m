function [H_fbb_select, a_TX_fbb_select, a_RX_fbb_select, Fbb_select] = select_Fbb(H_user_select, a_TX_user_select, a_RX_user_select,users,Num_users);

     for u = 1:1:users
     H(u,:,:) = H_user_select(u,:,:);
     Frf(:, u) = a_TX_user_select(:, u);
     Wrf(:, u) = a_RX_user_select(:, u);
     end 
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


end

