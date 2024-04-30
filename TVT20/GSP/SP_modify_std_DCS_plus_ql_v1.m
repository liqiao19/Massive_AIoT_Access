function [bbb mbm_det used_iter] = SP_modify_std_DCS_plus_ql_v1(z,Phi,support_BSAMP,opt)
%v1改得简单一点。z的第二个维度是1.
% CoSaMP for sparse signal recovery
% inputs: s-->spasity
%         z-->measurements
%         Phi-->mensurement matrix
%         max_iter --> max iteration

% output: a-->the recovery signal

% CoSaMP: ACHA (2009) 301-321
d = size(Phi,2); % d-dimensinal signal
active_num=d/opt.M;
% Nf = size(Phi,3);
%Begin CoSaMP
err = 10^(-6);  % residual error
a = zeros(d,1); % Initialize the recovery signal
v = z;
it=0;
stop = 0;
T_new =[];
while ~stop
    T_old = T_new;
    %     for pp=1:Nf
    %     y(:,pp) = Phi(:,:,pp)'*v(:,pp);
    %     end
    y = Phi'*v;
    %     y_temp = sum(abs(y).^2,2);
    y_temp = abs(y).^2;
    ix_temp = reshape(y_temp,opt.M,[]);
    [tmp, ix] = sort(abs(ix_temp), 'descend');
    ix_temp2 = ix(1:opt.Na,:);
    temp_ = [0:active_num-1]*opt.M;
    Omega = temp_ + ix_temp2;
    
    T = union(Omega, T_new);
    %Signal Estimation
    b = zeros(d,1);
    b(T) = Phi(:, T) \ z;
    %Prune
    b_temp = abs(b.^2);
    ix_temp2 = reshape(b_temp,opt.M,[]);
    [tmp2, ix2] = sort(abs(ix_temp2), 'descend');
    ix_temp3 = ix2(1:opt.Na,:);
    T_new = temp_ + ix_temp3;
    a = zeros(d, 1);
    a(T_new) = Phi(:, T_new) \ z;
    v = z - Phi(:,:)*a;
    
    it = it + 1;
    if (it >= opt.Max_iter|| isequal(T_old,T_new) ||  norm(v)<=err*norm(z))
        stop = 1;
    end
end %End CoSaMP iteration
used_iter = it;
support = T_new;
temppp_ = (0:length(support_BSAMP)-1)*opt.M;
mbm_det=support-temppp_-1;
bbb=a(support);