function [p_out, output] = globalFit(n_para,X,Y,tint,koff1_i,koff1_lb) 
%% Known Parameters
ttl = X(:,1);
ttl_temp = ttl;
lower_B = 1/min(Y(:,1));
upper_koff = 1/tint;
i_bt = 1;
%% The number of Parameters and Initialization
p_out = zeros(1,length(ttl)+8);
a_para = 100*ones(size(ttl));
weights = ones(size(X));
N_obs = numel(Y);
    if n_para == 2                    
            % p  = [kb,     koff,           a_para']
            para = [ 1, 	koff1_i,              a_para'];
            lb   = [ 0,     koff1_lb,       zeros(size(ttl_temp))'];
            ub   = [ Inf,   upper_koff,     Inf*ones(size(ttl_temp))'];      
        
            f1 = @(p)(   (model3(n_para,p,X,tint,ttl_temp)-Y).*weights );
            opts = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-9);
            %% Fit model to data
            [p, ~, residuals,~,output,~,~] = lsqnonlin(f1,para,lb,ub,opts);            
            disp('1 koff');
            output = output.message;
            N_para = n_para+length(ttl_temp);
            BIC1(i_bt,1) = N_obs*log(sum(sum(residuals.^2))/N_obs)...
                + N_para*(log(N_obs)-log(2*pi));
            p_out(1,:) = [p(1:2),1,zeros(1,4),p(3:end), BIC1];
    elseif n_para == 4
            para = [ 1,     koff1_i,    0.5,          2,          a_para'];
            lb   = [ 0,     koff1_lb,   lower_B,      koff1_lb,    zeros(size(ttl_temp))'];
            ub   = [ Inf,   upper_koff, 1-lower_B,    upper_koff, Inf*ones(size(ttl_temp))'];
            f1 = @(p)(   (model3(n_para,p,X,tint,ttl_temp)-Y).*weights );
            opts = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-9);
            %% Fit model to data
            [p, ~, residuals,~,output,~,~] = lsqnonlin(f1,para,lb,ub,opts);
            if p(2) > p(4)
                assign = p(2);
                p(2) = p(4);
                p(4) = assign;
                p(3) = 1 - p(3);
            end
            disp('2 koff');
            output = output.message;
            N_para = n_para+length(ttl_temp);
            BIC2(i_bt,1) = N_obs*log(sum(sum(residuals.^2))/N_obs)...
                + N_para*(log(N_obs)-log(2*pi));
            p_out(1,:) = [p(1:4),1-p(3),zeros(1,2),p(5:end), BIC2];    
    end
end

function f = model3(n_para,para,X,tint,ttl)
 
p = tint./ttl; p = p(:);
if n_para == 2
    kb = para(1);
    koff1 = para(2);
    ampl = para(3:end);
    f = (ampl'*ones(1,size(X,2))).*...
        (exp(-((kb.*p + koff1)*ones(1,size(X,2))).*X));   
elseif n_para == 4
    kb = para(1);    
    koff1 = para(2);
    B1 = para(3);
    koff2 = para(4);
    ampl = para(5:end);
    f = (ampl'*ones(1,size(X,2))).*...
        (B1.*exp(-((kb.*p + koff1)*ones(1,size(X,2))).*X)+...
        (1-B1) .* exp( -(kb.*p + koff2)*ones(1,size(X,2)).*X ));
elseif n_para ==6

    kb = para(1);
    koff1 = para(2);    B1 = para(3);
    koff2 = para(4);    B2 = para(5);
    koff3 = para(6);
    ampl = para(7:end);
    f = (ampl'*ones(1,size(X,2))).*...
        (B1.*exp(-((kb.*p + koff1) * ones(1,size(X,2))).*X)...
        + B2.* exp( -(kb.*p + koff2)*ones(1,size(X,2)).*X )+...
        (1-B1-B2).* exp( -(kb.*p + koff3)*ones(1,size(X,2)).*X ));
end   
end

