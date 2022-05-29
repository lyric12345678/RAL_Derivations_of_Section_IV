clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Objective: Proof of the derivation of our paper: "Semidefinite Optimization
%for 3-D Relative Pose Estimation". 
%Author: Ming Li
%Date: Feb./09/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proof of Appendix A
%% From equation (3) and equation (6) to equation (16)
syms w_l p C s_1l s_2l d_underline_l varepsilon_l
%We provide the result of d_underline_l-w_l.'*w_l
Formula_16_original=d_underline_l-(p+ C*s_2l-s_1l).'*(p+ C*s_2l-s_1l);
%Let's define
varepsilon_l=d_underline_l-p^2-s_1l^2-s_2l^2;                              %Note that d_underline_l=d_l^2-v_bar_l
%Our derivation result (16) is
Formula_16_Derivation=varepsilon_l+2*(s_1l-p).'*C*s_2l+2*s_1l.'*p;
%Then, we know that 
Result_1=simplify(Formula_16_Derivation-Formula_16_original);
% since C^2=1, we know that Formula_15_Derivation= Formula_15_original

%% From equation (16) to equation (17)
syms x y z q_4 q_1 q_2 q_3 x_1l y_1l z_1l x_2l y_2l z_2l d_meas r varepsilon_l
p = [x y z].';
% C = [1-2*q_2^2-2*q_3^2,      2*q_1*q_2-2*q_4*q_3,     2*q_1*q_3+2*q_4*q_2;
%      2*q_1*q_2+2*q_4*q_3,    1-2*q_1^2-2*q_3^2,       2*q_2*q_3-2*q_4*q_1;
%      2*q_1*q_3-2*q_4*q_2,    2*q_2*q_3+2*q_4*q_1,     1-2*q_1^2-2*q_2^2];
C = [q_1^2-q_2^2-q_3^2+q_4^2,      2*q_1*q_2-2*q_4*q_3,           2*q_1*q_3+2*q_4*q_2;
     2*q_1*q_2+2*q_4*q_3,          q_2^2-q_1^2-q_3^2+q_4^2,       2*q_2*q_3-2*q_4*q_1;
     2*q_1*q_3-2*q_4*q_2,          2*q_2*q_3+2*q_4*q_1,           q_3^2+q_4^2-q_1^2-q_2^2];
s_1l =[x_1l y_1l z_1l].';
s_2l =[x_2l y_2l z_2l].';

Formula_16_Derivation=varepsilon_l+2*(s_1l-p).'*C*s_2l+2*s_1l.'*p;
r=C.'*p;
e_l=vec(s_1l*s_2l.');
C_mathcal=vec(C);
Formula_17_Derivation=varepsilon_l+2*e_l.'*C_mathcal+2*[-s_2l.',s_1l.']*[r.',p.'].';
Result_2=simplify(Formula_17_Derivation-Formula_16_Derivation);
%Since Result_2=0, we know that Formula_16_Derivation=Formula_16_Derivation

%% Equation (18), determining matrix G
syms q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44
%We write G as follows.
q_tilde=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44].';
i=[1,5,9,2,4,3,7,6,8,1,5,9,6,8,3,7,1,5,9,2,4,1,5,9];
j=[1,1,1,2,2,3,3,4,4,5,5,5,6,6,7,7,8,8,8,9,9,10,10,10];
s=[1,-1,-1,2,2,2,2,2,-2,-1,1,-1,2,2,-2,2,-1,-1,1,2,-2,1,1,1];
G_sparse=full(sparse(i,j,s,9,10));
% G=[       1    0    0    0    -1   0    0    -1    0    1 ;
%           0    2    0    0     0    0    0    0    2    0;
%           0    0    2    0     0    0    -2   0    0    0;
%           0    2    0    0     0    0     0   0   -2    0;
%          -1    0    0    0     1    0     0  -1    0    1;
%           0    0    0    2     0    2     0   0    0    0;
%           0    0    2    0     0    0     2   0    0    0;
%           0    0    0   -2     0    2     0   0    0    0;
%          -1    0    0    0    -1    0     0   1    0    1];
%  G_sparse-G
C_vec=G_sparse*q_tilde;
C_Derivation=reshape(C_vec,3,3);
C_Origion=[q_11-q_22-q_33+q_44,      2*q_12-2*q_34,         2*q_13+2*q_23;
           2*q_12+2*q_34,        q_22+q_44-q_11-q_33,       2*q_23-2*q_14;
           2*q_13-2*q_24,        2*q_23+2*q_14,         q_33+q_44-q_11-q_22];
       

Result_3=simplify(C_Derivation-C_Origion);
%Since Result_3=0, we know that the matrix G_sparse is correct.

%% From (17) to (19) 
syms a_l
a_l=[2*e_l.'*G_sparse,-2*s_2l.',2*s_1l.',varepsilon_l].';
x_opt=[q_tilde.',r.',p.',1].';
Formula_18_Derivation=a_l.'*x_opt;
%We have used the relationship q_i*q_j=q_ij after (17). Therefore, we
%should make a simple change here. That is we replace C_vec with C_mathcal
%for the implementation of matlab
Formula_16_Derivation_Change=varepsilon_l+2*e_l.'*C_vec+2*[-s_2l.',s_1l.']*[r.',p.'].';
Result_4=simplify(Formula_18_Derivation-Formula_16_Derivation_Change);
%Since Result_3=0, we know that Formula_18_Derivation=Formula_16_Derivation


%% 28 Constraints in Appendix B
%% Constraint 1
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%First Constraint q_11*q_44=q_14^2;
M_1=sparse([1 4],[10 4],[1 -1],17,17);
Constraint_1_Original=q_11*q_44-q_14^2;
Constraint_1_Derivation=x_opt.'*M_1*x_opt;
Result_Constraint_1=simplify(Constraint_1_Original-Constraint_1_Derivation);

%% Constraint 2
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Second Constraint q_22*q_44=q_24^2;
M_2=sparse([5 7],[10 7],[1 -1],17,17);
Constraint_2_Original=q_22*q_44-q_24^2;
Constraint_2_Derivation=x_opt.'*M_2*x_opt;
Result_Constraint_2=simplify(Constraint_2_Original-Constraint_2_Derivation);

%% Constraint 3
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Third Constraint q_33*q_44=q_34^2;
M_3=sparse([8 9],[10 9],[1 -1],17,17);
Constraint_3_Original=q_33*q_44-q_34^2;
Constraint_3_Derivation=x_opt.'*M_3*x_opt;
Result_Constraint_3=simplify(Constraint_3_Original-Constraint_3_Derivation);


%% Constraint 4
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Fourth Constraint q_11*q_22=q_12^2;
M_4=sparse([1 2],[5 2],[1 -1],17,17);
Constraint_4_Original=q_11*q_22-q_12^2;
Constraint_4_Derivation=x_opt.'*M_4*x_opt;
Result_Constraint_4=simplify(Constraint_4_Original-Constraint_4_Derivation);

%% Constraint 5
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Fifth Constraint q_22*q_33=q_22^2;
M_5=sparse([5 6],[8 6],[1 -1],17,17);
Constraint_5_Original=q_22*q_33-q_23^2;
Constraint_5_Derivation=x_opt.'*M_5*x_opt;
Result_Constraint_5=simplify(Constraint_5_Original-Constraint_5_Derivation);

%% Constraint 6
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Seventh Constraint q_33*q_11=q_13^2;
M_6=sparse([8 3],[1 3],[1 -1],17,17);
Constraint_6_Original=q_33*q_11-q_13^2;
Constraint_6_Derivation=x_opt.'*M_6*x_opt;
Result_Constraint_6=simplify(Constraint_6_Original-Constraint_6_Derivation);

%% Constraint 7
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Eigth Constraint q_14*q_24=q_12*q_44;
M_7=sparse([4 2],[7 10],[1 -1],17,17);
Constraint_7_Original= q_14*q_24-q_12*q_44;
Constraint_7_Derivation=x_opt.'*M_7*x_opt;
Result_Constraint_7=simplify(Constraint_7_Original-Constraint_7_Derivation);

%% Constraint 8
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Eighth Constraint q_24*q_34=q_23*q_44;
M_8=sparse([7 6],[9 10],[1 -1],17,17);
Constraint_8_Original= q_24*q_34-q_23*q_44;
Constraint_8_Derivation=x_opt.'*M_8*x_opt;
Result_Constraint_8=simplify(Constraint_8_Original-Constraint_8_Derivation);

%% Constraint 9
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Tenth Constraint q_34*q_14=q_13*q_44;
M_9=sparse([9 3],[4 10],[1 -1],17,17);
Constraint_9_Original= q_34*q_14-q_13*q_44;
Constraint_9_Derivation=x_opt.'*M_9*x_opt;
Result_Constraint_9=simplify(Constraint_9_Original-Constraint_9_Derivation);

%% Constraint 10
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Tenth Constraint q_12*q_13=q_11*q_23;
M_10=sparse([2 1],[3 6],[1 -1],17,17);
Constraint_10_Original= q_12*q_13-q_11*q_23;
Constraint_10_Derivation=x_opt.'*M_10*x_opt;
Result_Constraint_10=simplify(Constraint_10_Original-Constraint_10_Derivation);


%% Constraint 11
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Eleventh Constraint q_12*q_23=q_22*q_13;
M_11=sparse([2 5],[6 3],[1 -1],17,17);
Constraint_11_Original=q_12*q_23-q_22*q_13;
Constraint_11_Derivation=x_opt.'*M_11*x_opt;
Result_Constraint_11=simplify(Constraint_11_Original-Constraint_11_Derivation);

%% Constraint 12
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Twelfth Constraint q_13*q_23=q_33*q_12;
M_12=sparse([3 8],[6 2],[1 -1],17,17);
Constraint_12_Original=q_13*q_23-q_33*q_12;
Constraint_12_Derivation=x_opt.'*M_12*x_opt;
Result_Constraint_12=simplify(Constraint_12_Original-Constraint_12_Derivation);

%% Constraint 13
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Thirteen Constraint q_12*q_34=q_13*q_24;
M_13=sparse([2 3],[9 7],[1 -1],17,17);
Constraint_13_Original=q_12*q_34-q_13*q_24;
Constraint_13_Derivation=x_opt.'*M_13*x_opt;
Result_Constraint_13=simplify(Constraint_13_Original-Constraint_13_Derivation);

%% Constraint 14
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Fifteenth Constraint q_23*q_34=q_33*q_23;
M_14=sparse([6 8],[9 7],[1 -1],17,17);
Constraint_14_Original=q_23*q_34-q_33*q_24;
Constraint_14_Derivation=x_opt.'*M_14*x_opt;
Result_Constraint_14=simplify(Constraint_14_Original-Constraint_14_Derivation);

%% Constraint 15
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Fifteenth Constraint q_13*q_34=q_33*q_14;
M_15=sparse([3 8],[9 4],[1 -1],17,17);
Constraint_15_Original=q_13*q_34-q_33*q_14;
Constraint_15_Derivation=x_opt.'*M_15*x_opt;
Result_Constraint_15=simplify(Constraint_15_Original-Constraint_15_Derivation);

%% Constraint 16
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
% Sixteenth Constraint q_23*q_24=q_22*q_34;
M_16=sparse([6 5],[7 9],[1 -1],17,17);
Constraint_16_Original=q_23*q_24-q_22*q_34;
Constraint_16_Derivation=x_opt.'*M_16*x_opt;
Result_Constraint_16=simplify(Constraint_16_Original-Constraint_16_Derivation);

%% Constraint 17
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
% Seventeenth Constraint q_12*q_24=q_22*q_14;
M_17=sparse([2 5],[7 4],[1 -1],17,17);
Constraint_17_Original=q_12*q_24-q_22*q_14;
Constraint_17_Derivation=x_opt.'*M_17*x_opt;
Result_Constraint_17=simplify(Constraint_17_Original-Constraint_17_Derivation);

%% Constraint 18
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
% Eighteenth Constraint q_13*q_14=q_11*q_34;
M_18=sparse([3 1],[4 9],[1 -1],17,17);
Constraint_18_Original=q_13*q_14-q_11*q_34;
Constraint_18_Derivation=x_opt.'*M_18*x_opt;
Result_Constraint_18=simplify(Constraint_18_Original-Constraint_18_Derivation);

%% Constraint 19
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Nineteenth Constraint q_12*q_14=q_11*q_24;
M_19=sparse([2 1],[4 7],[1 -1],17,17);
Constraint_19_Original=q_12*q_14-q_11*q_24;
Constraint_19_Derivation=x_opt.'*M_19*x_opt;
Result_Constraint_19=simplify(Constraint_19_Original-Constraint_19_Derivation);

%% Constraint 20
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%Twentith Constraint q_12*q_34=q_14*q_23;
M_20=sparse([2 4],[9 6],[1 -1],17,17);
Constraint_20_Original=q_12*q_34-q_14*q_23;
Constraint_20_Derivation=x_opt.'*M_20*x_opt;
Result_Constraint_20=simplify(Constraint_20_Original-Constraint_20_Derivation);

%% Constraint 21
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%21st Constraint ;
syms r_1 r_2 r_3
x_opt_change=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1].';
% r_change=C_Origion.'*p;
% simplify((G_sparse([1,2,3],:)*q_tilde).'*p-r_change(1))
Constraint_21_Original=collect((G_sparse([1,2,3],:)*q_tilde).'*p-r_1,x_opt_change);
M_21=sparse([1 2 3 5 7 8 9 10 11],[14 15 16 14 16 14 15 14 17],[1 2 2 -1 -2 -1 2 1 -1],17,17);
Constraint_21_Derivation=x_opt_change.'*M_21*x_opt_change;
Result_Constraint_21=simplify(Constraint_21_Original-Constraint_21_Derivation);

%% Constraint 22
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%22nd Constraint ;
Constraint_22_Original=collect((G_sparse([4,5,6],:)*q_tilde).'*p-r_2,x_opt_change);
M_22=sparse([2 1 4 5 6 8 9 10 12],[14 15 16 15 16 15 14 15 17],[2 -1 2 1 2 -1 -2 1 -1],17,17);
Constraint_22_Derivation=x_opt_change.'*M_22*x_opt_change;
Result_Constraint_22=simplify(Constraint_22_Original-Constraint_22_Derivation);

%% Constraint 23
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%23rd Constraint ;
Constraint_23_Original=collect((G_sparse([7,8,9],:)*q_tilde).'*p-r_3,x_opt_change);
M_23=sparse([1 3 4 5 6 7 8 10 13],[16 14 15 16 15 14 16 16 17],[-1 2 -2 -1 2 2 1 1 -1],17,17);
Constraint_23_Derivation=x_opt_change.'*M_23*x_opt_change;
Result_Constraint_23=simplify(Constraint_23_Original-Constraint_23_Derivation);


%% Constraint 24
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%24th Constraint p-Cr=0
r_change=[r_1 r_2 r_3].';
Constraint_24_Original=collect((G_sparse([1,4,7],:)*q_tilde).'*r_change-x,x_opt_change);
M_24=sparse([1 2 3 5 7 8 9 10 14],[11 12 13 11 13 11 12 11 17],[1 2 2 -1 2 -1 -2 1 -1],17,17);
Constraint_24_Derivation=x_opt_change.'*M_24*x_opt_change;
Result_Constraint_24=simplify(Constraint_24_Original-Constraint_24_Derivation);

%% Constraint 25
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%25th Constraint
Constraint_25_Original=collect((G_sparse([2,5,8],:)*q_tilde).'*r_change-y,x_opt_change);
M_25=sparse([1 2 4 5 6 8 9 10 15],[12 11 13 12 13 12 11 12 17],[-1 2 -2 1 2 -1 2 1 -1],17,17);
Constraint_25_Derivation=x_opt_change.'*M_25*x_opt_change;
Result_Constraint_25=simplify(Constraint_25_Original-Constraint_25_Derivation);

%% Constraint 26
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%26th Constraint
Constraint_26_Original=collect((G_sparse([3,6,9],:)*q_tilde).'*r_change-z,x_opt_change);
M_26=sparse([1 3 4 5 6 7 8 10 16],[13 11 12 13 12 11 13 13 17],[-1 2 2 -1 2 -2 1 1 -1],17,17);
Constraint_26_Derivation=x_opt_change.'*M_26*x_opt_change;
Result_Constraint_26=simplify(Constraint_26_Original-Constraint_26_Derivation);

%% Constraint 27
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%27th Constraint p^2=d0^2;
syms d0
d0=10;
M_27=sparse([14 15 16 17],[14 15 16 17],[1 1 1 -d0^2],17,17);                 
Constraint_27_Original=p.'*p-d0^2;
Constraint_27_Derivation=x_opt.'*M_27*x_opt;
Result_Constraint_27=simplify(Constraint_27_Original-Constraint_27_Derivation);


%% Constraint 28
%we write out x_opt=[q_11 q_12 q_13 q_14 q_22 q_23 q_24 q_33 q_34 q_44 r_1 r_2 r_3 x y z 1]
%28th Constraint
Constraint_28_Original=collect((q_11+q_22+q_33+q_44-1)^2,x_opt_change);
M_28=sparse([1 1 1 1 1 5 5 5 5 8 8 8 10 10 17],[1 5 8 10 17 5 8 10 17 8 10 17 10 17 17],[1 2 2 2 -2 1 2 2 -2 1 2 -2 1 -2 1],17,17);
Constraint_28_Derivation=x_opt_change.'*M_28*x_opt_change;
Result_Constraint_28=simplify(Constraint_28_Original-Constraint_28_Derivation);






