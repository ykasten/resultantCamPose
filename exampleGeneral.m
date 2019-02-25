% This is an example of using the code from the paper "Resultant
% Based Incremental Recovery of Camera Pose from Pairwise Matches"
% (WACV'19)
% for the general case: up to 3 point correspondences come from the same image of known
% camera
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019

clear
close all
%Loading a set of 10 configurations of the problem:
load configurations

errorst=ones(10,1);
errorsR=ones(10,1);
for j=1:10
    configuration=configurations(j);

    %Solve the problem 
    [ Rs1,ts1 ,ok,sols] = getRtResultant( configuration.p1s,configuration.qs,configuration.Rs,configuration.ts );
    
    %Choose the best solution comparing to the ground truth:
    [~,indddd]= min(sum((sols(:,1:3)-repmat(configuration.qg',size(sols,1),1)).^2,2));
    OurR=Rs1{indddd};
    Ourt=ts1{indddd};
    
    errorst(j)=norm(OurR-configuration.Rg,'fro');
    errorsR(j)=norm(Ourt-configuration.tg,'fro');
end

mean(errorst)
mean(errorsR)

