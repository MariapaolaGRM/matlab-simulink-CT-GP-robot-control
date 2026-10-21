function f = gp_mean(x, P)
%#codegen
% x = [ddq1; ddq2; dq1; dq2; q1; q2]  (stesso ordine delle colonne di D_train.X)
z  = (x(:).' - P.mu)./P.sg;              % 1x6
d2 = sum((P.Z - z).^2, 2);               % m x 1
f  = zeros(2,1);
for i = 1:2
    k    = P.sF(i)^2 * exp(-0.5*d2/P.sL(i)^2);
    f(i) = k.' * P.alpha(:,i);           % media del GP, eq. (1)
end
end