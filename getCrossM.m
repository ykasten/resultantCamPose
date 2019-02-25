function [ res ] = getCrossM( v )
%GETCROSSM get cross product matrix from a vector
% Copyright (C) Yoni Kasten, Weizmann Institute, 2019
res=[0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];


end

