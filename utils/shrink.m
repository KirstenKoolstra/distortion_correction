function [ xs ] = shrink( x,parameter )
%   Detailed explanation goes here

s = abs(x);
ss = s-parameter;
ss = ss.*(ss>0); %take max

s = s+(s<parameter);
ss = ss./s;

xs = ss.*x;

end

