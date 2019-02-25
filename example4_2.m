% This is an example of using the code from the paper "Resultant
% Based Incremental Recovery of Camera Pose from Pairwise Matches"
% (WACV'19)
% for the case 4+2: 4 point correspondences come from the image of same known
% camera
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019

clear
close all;
%Loading a set of 10 configurations of the problem:
load configurations4_2.mat



errorst=ones(10,1);
errorsR=ones(10,1);


for j=1:10
    configuration=configurations(j);
 
    %Solve the problem
    [ Rs1,ts1 ,ok,sols] = getRtResultant_4_2( configuration.p1s,configuration.qs,configuration.Rs,configuration.ts );
    
    
    %In the case of 4+2 there are redundant solutions that should be
    %dismissed.
    %Assuming a calibrated setup, we have the calibration matrix of the cameras.
    %In the case of the example all the cameras has the same calibration matrix.
    K=[2759.48 0 1520.69
        0 2764.16 1006.81
        0 0 1];
    
    %Filter solutions using the Sampson distance metric
    [ sols,ts1,Rs1 ] = filterSampson( configuration.p1s,configuration.qs,configuration.Rs,configuration.ts,{K,K,K,K,K,K} ,ts1,Rs1,K,sols);
    
    
    %Choose the best solution comparing to the ground truth:
    [~,indddd]= min(sum((sols(:,1:3)-repmat(configuration.qg',size(sols,1),1)).^2,2));
    
    ourquant=sols(indddd,1:3)';
    OurR=Rs1{indddd};
    Ourt=ts1{indddd};
    
    errorst(j)=norm(OurR-configuration.Rg,'fro');
    errorsR(j)=norm(Ourt-configuration.tg,'fro');
end
mean(errorst)
mean(errorsR)

