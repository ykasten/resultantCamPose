function [ inputs ] = extractCoefficients( p1s,qs,Rs,ts )
%EXTRACTCOEFFICIENTS takes a configurations of 6 correspondences  from 6
%cameras to an unknown camera and returns the coefficients for the polynomial system
%p1s - a cell of size 6x1 containing the points in the image of the new unknown camera
%qs - a cell of size 6x1 containing the points in the images of the known
%cameras
%Rs -  a cell of size 6x1 containing the rotations of the known
%cameras
%ts -  a cell of size 6x1 containing the positions of the known
%cameras
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019

q=Rs{1}*qs{1};
p1_x=q(1);
p1_y=q(2);
p1_z=q(3);

q=Rs{2}*qs{2};
p2_x=q(1);
p2_y=q(2);
p2_z=q(3);

q=Rs{3}*qs{3};
p3_x=q(1);
p3_y=q(2);
p3_z=q(3);


q=Rs{4}*qs{4};
p4_x=q(1);
p4_y=q(2);
p4_z=q(3);

q=Rs{5}*qs{5};
p5_x=q(1);
p5_y=q(2);
p5_z=q(3);

q=Rs{6}*qs{6};
p6_x=q(1);
p6_y=q(2);
p6_z=q(3);


p1_x_t=p1s{1}(1);
p1_y_t=p1s{1}(2);

p2_x_t=p1s{2}(1);
p2_y_t=p1s{2}(2);

p3_x_t=p1s{3}(1);
p3_y_t=p1s{3}(2);

p4_x_t=p1s{4}(1);
p4_y_t=p1s{4}(2);

p5_x_t=p1s{5}(1);
p5_y_t=p1s{5}(2);

p6_x_t=p1s{6}(1);
p6_y_t=p1s{6}(2);

t1_1=ts{1}(1) ;t1_2=ts{1}(2); t1_3=ts{1}(3);
t2_1=ts{2}(1)  ;t2_2=ts{2}(2) ; t2_3=ts{2}(3) ;
t3_1=ts{3}(1) ;t3_2=ts{3}(2); t3_3=ts{3}(3);

t4_1=ts{4}(1) ;t4_2=ts{4}(2); t4_3=ts{4}(3);
t5_1=ts{5}(1) ;t5_2=ts{5}(2); t5_3=ts{5}(3);
t6_1=ts{6}(1) ;t6_2=ts{6}(2); t6_3=ts{6}(3);
inputs=[t1_1, t1_2, t1_3, t2_1, t2_2, t2_3, t3_1, t3_2, t3_3, t4_1, t4_2, t4_3, t5_1, t5_2, t5_3, t6_1, t6_2, t6_3, p1_x, p1_y, p1_z, p1_x_t, p1_y_t, p2_x, p2_y, p2_z, p2_x_t, p2_y_t, p3_x, p3_y, p3_z, p3_x_t, p3_y_t, p4_x, p4_y, p4_z, p4_x_t, p4_y_t, p5_x, p5_y, p5_z, p5_x_t, p5_y_t, p6_x, p6_y, p6_z, p6_x_t, p6_y_t];


end

