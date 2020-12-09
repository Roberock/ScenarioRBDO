function compensator=build_comp(gains)

[A,B,C,D]=tf2ss(gains.num,gains.den);
compensator.A=A;
compensator.B=B;
compensator.C=C;
compensator.D=D;