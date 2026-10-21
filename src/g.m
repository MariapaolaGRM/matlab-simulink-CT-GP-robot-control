function gvec = g(q,p)

    q1 = q(1);
    q2 = q(2);

    m1  = p.ml1;
    m2  = p.ml2;
    mm2 = p.mm2;

    l1 = p.l1;
    l2 = p.l2;
    a1 = p.a1;
    g0 = p.g;

    gvec = [ ...
        m1*g0*l1*cos(q1) ...
        + m2*g0*(a1*cos(q1) + l2*cos(q1+q2)) ...
        + mm2*g0*a1*cos(q1);

        m2*g0*l2*cos(q1+q2)
    ];

end