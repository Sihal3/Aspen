% NASA CEA Output Reader
% Pressures: 140:10:300
% O/F: 1:.1:3

largeData = readmatrix('nasa_cea_output1.txt');
aeat = largeData(:,1);
rho = largeData(:,2);
gas_pressure = largeData(:,3);
cp = largeData(:,4);
temp = largeData(:,5);
isp = largeData(:,6);
ae = largeData(:,7);
cf = largeData(:,8);

chamber_pressure = (140:10:300);
of_ratio = (1:.1:3);

