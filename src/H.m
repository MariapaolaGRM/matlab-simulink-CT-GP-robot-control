function Hmat = H(q,p)

    q2 = q(2);

    m1  = p.ml1;
    m2  = p.ml2;
    mm1 = p.mm1;
    mm2 = p.mm2;

    l1 = p.l1;
    l2 = p.l2;

    I1  = p.Il1;
    I2  = p.Il2;
    Im1 = p.Im1;
    Im2 = p.Im2;

    a1 = p.a1;
    kr1 = p.kr1;
    kr2 = p.kr2;

    JP_l1 = [-l1*sin(q(1)), 0;
               l1*cos(q(1)), 0];

    JO_l1 = [0, 0;
             0, 0;
             1, 0];

    JP_l2 = [-a1*sin(q(1))-l2*sin(q(1)+q2), -l2*sin(q(1)+q2);
              a1*cos(q(1))+l2*cos(q(1)+q2),  l2*cos(q(1)+q2)];

    JO_l2 = [0, 0;
             0, 0;
             1, 1];

    JP_m1 = [0, 0;
             0, 0];

    JO_m1 = [0, 0;
             0, 0;
             kr1, 0];

    JP_m2 = [-a1*sin(q(1)), 0;
              a1*cos(q(1)), 0];

    JO_m2 = [0, 0;
             0, 0;
             1, kr2];

    B1 = m1  * JP_l1.' * JP_l1 ...
       + I1  * JO_l1.' * JO_l1 ...
       + mm1 * JP_m1.' * JP_m1 ...
       + Im1 * JO_m1.' * JO_m1;

    B2 = m2  * JP_l2.' * JP_l2 ...
       + I2  * JO_l2.' * JO_l2 ...
       + mm2 * JP_m2.' * JP_m2 ...
       + Im2 * JO_m2.' * JO_m2;

    Hmat = B1 + B2;

end