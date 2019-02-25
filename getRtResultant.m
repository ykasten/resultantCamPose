function [Rs1,ts1,ok ,solsres ] = getRtResultant( p1s,qs,Rs,ts )
%GETRT  This is an implementation for the method from the paper "Resultant
% Based Incremental Recovery of Camera Pose from Pairwise Matches"
% (WACV'19)
%Given 6 correspondences from images with known cameras to an image 
%of unknown camera s.t no more the 3 correspondences come from the same known camera
%return the position and orientation of the unknown camera
%p1s - a cell of size 6x1 containing the points in the image of the new camera
%qs - a cell of size 6x1 containing the points in the images of the known
%cameras
%Rs -  a cell of size 6x1 containing the rotations of the known
%cameras
%ts -  a cell of size 6x1 containing the positions of the known
%cameras
%outputs:
%Rs1 - a cell of size mx1 of solutions for the rotation of the unknown
%camera
%ts1 - a cell of size mx1 of solutions for the translation of the unknown
%camera
%ok - a binary variable that is 0, if something went wrong
%solsres - mx7 matrix representing m solutions for the positions and
%orientations of the camera with the unknown position and orientation. 
%Each row represents one solution where the first 3 entries are the 
%normalized quaternion, and the last 4 are d.  
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019

[ inputs ] = extractCoefficients( p1s,qs,Rs,ts );
[solsres,ok]=solveEquations(inputs);

if ok==0
    Rs1=0;ts1=0; solsres=0;
    return;
end
Rs1=cell(size(solsres,1),1);
ts1=cell(size(solsres,1),1);
%Build the rotations and the translations from the solutions for  d and q
for i=1:size(solsres,1)
    q=[1 solsres(i,1:3)];
    q=q/norm(q);
    d=solsres(i,4:7);
    qstar=[q(1) -q(2:4)];
    t=hamiltonProduct(d,qstar)/norm([1 solsres(i,1:3)]);
    xx=solsres(i,1:3)';
    curR= 2*(xx*xx'-[0 -xx(3) xx(2);xx(3) 0 -xx(1);-xx(2) xx(1) 0])+(1-xx'*xx)*eye(3);
    ss=curR*curR';
    R=curR'/sqrt(ss(1,1));
   
    Rs1{i}=R;
    ts1{i}=t(2:4)';
end
end


function [sols,ok]=solveEquations(inputs)

%Extract coefficient matrices of the Dixon resultant (Eq. 16 in the paper)
[m1,m2,m3,m4,m5,m6,m7,m8,m9,~]=ResultantExtractor(2,inputs);

if rank(m1)<27
    sols=0;
    ok=0;
    return;
end
ok=1;

%Building the generalized Eigen values system (Eq. 18  in the paper)
C1=[zeros(27*7,27) eye(27*7);-m1 -m2 -m3 -m4 -m5 -m6 -m7 -m8];
C2=eye(27*8);
C2(7*27+1:8*27,7*27+1:8*27)=m9;
[V2,D2] = eig(C1,C2);

solutions=diag(D2);

V2=V2(:,(abs(solutions)<10^3));
solutions=solutions((abs(solutions)<10^3));




%Taking the real solutions
V2real=real(V2(:,abs(imag(solutions))<0.000001));
solutionsReal=real(solutions(abs(imag(solutions))<0.000001));

q2s=solutionsReal.';
q3s=V2real(3,:)./V2real(1,:);
q4s=V2real(2,:)./V2real(1,:);

%Given q2,q3,q4 get d linearly
[A,b]=extractLinearEquationsT(inputs,q2s,q3s,q4s);
s=A\b;
sols=[q2s' q3s' q4s' reshape(s,4,size(solutionsReal,1))'] ;


end
