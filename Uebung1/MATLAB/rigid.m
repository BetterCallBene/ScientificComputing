function [ dy ] = rigid( t, y )
%RIGID Summary of this function goes here
%   Detailed explanation goes here
dy = zeros(3,1);    % a column vector
dy(1) = y(2) * y(3);
dy(2) = -y(1) * y(3);
dy(3) = -0.51 * y(1) * y(2);

end

