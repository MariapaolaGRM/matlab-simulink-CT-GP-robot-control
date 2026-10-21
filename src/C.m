function Cmat = C(q,dq,p)

    q2  = q(2);
    dq1 = dq(1);
    dq2 = dq(2);

    h = -p.ml2*p.a1*p.l2*sin(q2);

    Cmat = [ h*dq2,          h*(dq1+dq2);
            -h*dq1,         0 ];

end