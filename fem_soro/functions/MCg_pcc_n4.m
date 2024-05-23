function [M,C,g] = MCg_pcc_n4(mu,L,Iz,gp,q,dq)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%
sq = 1/2 * sinc( q./(2*pi) );
xi = m_q2xi_n4(L(1),L(2),L(3),L(4),q(1),q(2),q(3),q(4),sq(1),sq(2),sq(3),sq(4));

%
M_xi = M_xi_h4_n4(Iz(1),Iz(2),Iz(3),Iz(4),mu(1),mu(2),mu(3),mu(4),...
    xi(2),xi(3),xi(4),xi(5),xi(6),xi(7),xi(8),xi(9),xi(10),xi(11),xi(12),xi(13),xi(14));

%
Lc = zeros(size(L));
ind_q = (q~=0);
Lc(ind_q) = (q(ind_q).*cos(q(ind_q)./2)-2*sin(q(ind_q)./2))./(2*q(ind_q).^2);
Lc = L .* Lc;
J_m = Jm_n4(Lc(1),Lc(2),Lc(3),Lc(4));

%
M = J_m'*M_xi*J_m;

%
dLc = zeros(size(L));
dLc(ind_q) = -((q(ind_q).^2 - 8) .* sin(q(ind_q)./2) + 4 * q(ind_q) .* cos(q(ind_q)./2)) ./ (4*q(ind_q).^3);
dLc(~ind_q) = -1/24;
dLc = L .* dLc;
dJ_m = dJm_n4(dLc(1),dLc(2),dLc(3),dLc(4),dq(1),dq(2),dq(3),dq(4));

%
dxi = J_m * dq;
C_xi = C_xi_h4_n4(dxi(1),dxi(2),dxi(3),dxi(4),dxi(5),dxi(6),dxi(7),dxi(8),dxi(9),dxi(10),dxi(11),dxi(12),dxi(13),dxi(14),...
    mu(1),mu(2),mu(3),mu(4),...
    xi(2),xi(3),xi(4),xi(5),xi(6),xi(7),xi(8),xi(9),xi(10),xi(11),xi(12),xi(13),xi(14));
C = J_m'*(M_xi*dJ_m + C_xi*J_m);

%
g_xi = g_xi_h4_n4(gp,mu(1),mu(2),mu(3),mu(4),xi(1),xi(2),xi(3),xi(4),xi(5),xi(6),xi(7),xi(8),xi(9),xi(10),xi(11),xi(12),xi(13),xi(14));
g = J_m'*g_xi;

end