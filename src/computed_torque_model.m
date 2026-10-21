function tau = computed_torque_model(q,dq,ddq_d,dq_d)

    p = get_params(false);
    tau = H(q,p)*ddq_d + C(q,dq,p)*dq_d + g(q,p);
    
end