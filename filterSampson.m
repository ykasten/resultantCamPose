function [ solsres,ts1,Rs1 ] = filterSampson( p1s,qs,Rs,ts,Ks ,ts1,Rs1,K1,solsres)
%FILTERSAMPSON gets solutions for the polynomial system and filter wrong
%solutions using Sampson distance between the corresponding points  
%
%p1s - a cell of size 6x1 containing the points in the image of the new camera
%qs - a cell of size 6x1 containing the points in the images of the known
%cameras
%Rs -  a cell of size 6x1 containing the rotations of the known
%cameras
%ts -  a cell of size 6x1 containing the positions of the known
%cameras
%Ks - a cell of size 6x1 containing the calibration matrices of the known
%cameras
%ts1 - a cell of size mx1 of solutions for the translations of the unknown
%camera
%Rs1 - a cell of size mx1 of solutions for the rotation of the unknown
%camera
%K1 - the calibration of the camera with the unknown position and orientation
%(the calibration assumed to be known)
%solsres - mx7 matrix representing m solutions for the positions and
%orientations of the camera with the unknown position and orientation. 
%Each row represents one solution where the first 3 entries are the 
%normalized quaternion, and the last 4 are d.  
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019

for i=1:length(solsres)
    
    %For each solution for the unknown camera, compute the fundamental 
    %matrix relative the each known camera, and calculate the Sampson distance for the correspondence  
    
    F=inv(K1')*Rs1{i}'*(getCrossM(ts1{i})-getCrossM(ts{1}))*Rs{1}*inv(Ks{1});
    [ er1 ] = sampsonError( F,Ks{1}*qs{1},K1*p1s{1} );
    
    
    F=inv(K1')*Rs1{i}'*(getCrossM(ts1{i})-getCrossM(ts{2}))*Rs{2}*inv(Ks{2});
    [ er2 ] = sampsonError( F,Ks{2}*qs{2},K1*p1s{2} );
    
    F=inv(K1')*Rs1{i}'*(getCrossM(ts1{i})-getCrossM(ts{3}))*Rs{3}*inv(Ks{3});
    [ er3 ] = sampsonError( F,Ks{3}*qs{3},K1*p1s{3} );
    F=inv(K1')*Rs1{i}'*(getCrossM(ts1{i})-getCrossM(ts{4}))*Rs{4}*inv(Ks{4});
    [ er4 ] = sampsonError( F,Ks{4}*qs{4},K1*p1s{4} );
    
    F=inv(K1')*Rs1{i}'*(getCrossM(ts1{i})-getCrossM(ts{5}))*Rs{5}*inv(Ks{5});
    [ er5 ] = sampsonError( F,Ks{5}*qs{5},K1*p1s{5} );
    F=inv(K1')*Rs1{i}'*(getCrossM(ts1{i})-getCrossM(ts{6}))*Rs{6}*inv(Ks{6});
    [ er6 ] = sampsonError( F,Ks{6}*qs{6},K1*p1s{6} );
    
    
    
    errorss(i)=norm([er1;er2;er3;er4;er5;er6]);

end
indsKeep=errorss<10^-4;
solsres=solsres(indsKeep,:);
ts1=ts1(indsKeep);
Rs1=Rs1(indsKeep);
end

