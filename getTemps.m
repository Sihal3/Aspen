function [chamber_temp, throat_temp] = getTemps(pressure, temp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

start = ((pressure - 160) / 10) * 2 + 1;
index = start;
chamber_temp = zeros(16,1);
throat_temp = zeros(16,1);
for i = 1:16

    chamber_temp(i,1) = temp(index,1);
    throat_temp(i,1) = temp(index + 1,1);

    index = index + 30;


end