function [f g gv gpb] = CWRDistance(A, B, Va, Vb, Pa, Pb, scale, lambda)
[f1 g1 gv1 gpb1] = GaborTransform(A, A, Va, Va, Pa, Pa, lambda, scale);
[f2 g2 gv2 gpb2] = GaborTransform(A, B, Va, Vb, Pa, Pb, lambda, scale);
[f3] = GaborTransform(B, B, Vb, Vb, Pb, Pb, lambda, scale);
f = f1 - 2*f2 + f3;
g = 2*g1 - 2*g2;
gv = 2*gv1 - 2*gv2;
gpb = 2*gpb1 - 2*gpb2;