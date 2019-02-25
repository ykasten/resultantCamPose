function [Rs1,ts1,ok ,solsres ,quans] = getRtResultant_4_2( p1s,qs,Rs,ts )
%GETRT  This is an implementation for the method from the paper "Resultant
% Based Incremental Recovery of Camera Pose from Pairwise Matches"
% (WACV'19)
%Given 6 correspondences from images with known cameras to an image 
%of unknown camera s.t 4 correspondences come from the same known camera
%return the position and orientation of the unknown camera
%p1s - a cell of size 6x1 containing the points in the image of the new camera
%qs - a cell of size 6x1 containing the points in the images of the known
%cameras
%Rs -  a cell of size 6x1 containing the rotations of the known
%cameras
%ts -  a cell of size 6x1 containing the positions of the known
%cameras
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019

[ inputs ] = extractCoefficients( p1s,qs,Rs,ts );
[solsres,ok]=solveEquations_4_2(inputs);
if ok==0
    Rs1=0;ts1=0; solsres=0;
    return;
end
Rs1=cell(size(solsres,1),1);
ts1=cell(size(solsres,1),1);
quans=cell(size(solsres,1),1);

%For each solution, return t and R
for i=1:size(solsres,1)
    q=[1 solsres(i,1:3)];
    quans{i}=solsres(i,1:3)';
    q=q/norm(q);
    xx=solsres(i,1:3)';
    curR= 2*(xx*xx'-[0 -xx(3) xx(2);xx(3) 0 -xx(1);-xx(2) xx(1) 0])+(1-xx'*xx)*eye(3);
    ss=curR*curR';
    R=curR'/sqrt(ss(1,1));
    d=solsres(i,4:7);
    qstar=[q(1) -q(2:4)];
    t=hamiltonProduct(d,qstar)/norm([1 solsres(i,1:3)]);
    Rs1{i}=R;
    ts1{i}=t(2:4)';
end



end

function [sols,ok]=solveEquations_4_2(inputs)

%Extract coefficient matrices of the Dixon resultant (Eq. 16 in the paper)
[m1,m2,m3,m4,m5,m6,m7,m8,m9,~]=ResultantExtractor(2,inputs);


m1t=m1(1:23,1:23);
m2t=m2(1:23,1:23);
m3t=m3(1:23,1:23);
m4t=m4(1:23,1:23);
m5t=m5(1:23,1:23);
m6t=m6(1:23,1:23);
m7t=m7(1:23,1:23);
m8t=m8(1:23,1:23);
m9t=m9(1:23,1:23);

%Building the generalized Eigen values system (Eq. 18  in the paper)
C1=[zeros(23*7,23) eye(23*7);-m1t -m2t -m3t -m4t -m5t -m6t -m7t -m8t];
C2=eye(23*8);
C2(7*23+1:8*23,7*23+1:8*23)=m9t;

[V2,D2] = eig(C1,C2);

solutions=diag(D2);

V2=V2(:,(abs(solutions)<10^3));
solutions=solutions((abs(solutions)<10^3));



%Taking the real solutions
solutionsReal=real(solutions(abs(imag(solutions))<0.000001));
q4sa=[];
q3sa=[];
q2sa=[];
dsa=[];

for ii=1:length(solutionsReal)
    q2=solutionsReal(ii);
    
    %For each real solution of q2, find the solutions for q3,q4
    ress=m1+m2*q2+m3*q2^2+m4*q2^3+m5*q2^4+m6*q2^5+m7*q2^6+m8*q2^7+m9*q2^8;
    
    %Since the matrix is singular, use LU decomposition to extract the
    %solutions of q3,q4 (see "More than 3 correspondences from one camera"
    %in the paper)
    
    [~,u]=lu(ress(:,[3,5,6,8:10,12:15, 17:27 1,2,4,7,11,16]));
    q4s=roots(u(22,end:-1:end-5));
    
    q4s=real(q4s(abs(imag(q4s))<10^-7));
    [~,u]=lu(ress(:,[5,6,8:10,12:16, 17:27 1,2,3,4,7,11]));
    
    q3s=u(22,[[end-5 end-4] [end-2: end]])*[ones(size(q4s,1),1) q4s q4s.^2 q4s.^3 q4s.^4]';
    q3s=q3s'/(-u(22,end-3));
    
    %Given q2,q3,q4 get d linearly
    [Att,btt]=extractLinearEquationsT(inputs,repmat(q2,size(q3s,1),1)',q3s',q4s');
    s=Att\btt;
    
    ds=reshape(s,4,[])';
    
    
    
    q4sa=[q4sa;q4s];
    q3sa=[q3sa;q3s];
    q2sa=[q2sa;repmat(q2,length(q4s),1)];
    dsa=[dsa;ds];
    
    
    
end

sols=[q2sa q3sa q4sa dsa];
ok=1;


end

