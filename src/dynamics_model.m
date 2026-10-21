function ddq = dynamics_model(q,dq,tau)

    p = get_params(true);
    fr = p.Fv.*dq; %+ p.Fc.*tanh(dq/eps_f);
    ddq = H(q,p) \ (tau - C(q,dq,p)*dq - g(q,p) - fr);

end