function [q_meas,dq_meas,ddq_meas] = get_states_from_fem(N,Ncc,qout_full,dqout_full,soro,k,d,a,ind_meas)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 8
    ind_meas = 1:size(qout_full,2);
    calc_ddq = 0;
else
    calc_ddq = 1;
end

q_meas = zeros(Ncc,length(ind_meas));
dq_meas = zeros(size(q_meas));
ddq_meas = zeros(size(q_meas));

for n = 1:Ncc
    ni = 1 + round((n-1)*N/Ncc);
    nip1 = round(n*N/Ncc);
    q_meas(n,:) = sum(qout_full(ni:nip1,ind_meas));
    dq_meas(n,:) = sum(dqout_full(ni:nip1,ind_meas));
end

if calc_ddq
    for i = 1:length(ind_meas)
        q_i = qout_full(:,ind_meas(i));
        dq_i = dqout_full(:,ind_meas(i));
        ddq_i = FDab_wDK(soro, q_i, dq_i, -k*eye(N)*q_i - d*eye(N)*dq_i + a*ones(N,1));
        for n = 1:Ncc
            ddq_meas(n,i) = sum(ddq_i(1+round((n-1)*N/Ncc):round(n*N/Ncc)));
        end
    end
else
    ddq_meas = [];
end

end