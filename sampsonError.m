function [ errors ] = sampsonError( F,p,q )
%SAMPSONERROR Calculate the Sampson error associated with the 
%correspondences p<->q relative to the fundamental matrix F (q^TFp~0) 
% F - the fundamental matrix
% p - the points in the first image
% q - the points in the second image
% errors - the Sampson errors 
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019
first=F*p;
second=F'*q;

error=((sum(q.*first,1)).^2)./(first(1,:).^2+first(2,:).^2 +second(1,:).^2+second(2,:).^2);
errors=error';

end

