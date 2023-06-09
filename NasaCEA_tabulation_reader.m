% NASA CEA Output Reader
% Pressures: 140:10:300
% O/F: 1:.1:3

largeData = readmatrix('nasa_cea_output2.txt');
gamma = largeData(:,1);
rho = largeData(:,2);
gas_pressure = largeData(:,3);
cp = largeData(:,4);
temp = largeData(:,5);
isp = largeData(:,6);
ae = largeData(:,7);
cf = largeData(:,8);

chamber_pressure = (160:10:300);
of_ratio = (1:.2:4);

%of = (1:)

plot_count = 1;
index = 1;
for i = 1:7
    for j = 1:2
        %index = 1;

        subplot(7, 2, plot_count);
        plot_count = plot_count + 1;

        [chamber_temp, throat_temp] = getTemps(chamber_pressure(1,index), temp);
        index = index + 1;

        of_ratio = of_ratio';

        plot(of_ratio, chamber_temp, "-r");
        hold on
        plot(of_ratio, throat_temp, "-b");
        hold off;
        title("Pressure: " + chamber_pressure(1,index));
        xlabel("OF Ratio");
        ylabel("Temperature [K]");
        legend("Chamber Temp", "Throat Temp", "Location", "southeast");

     
    end
end

plot_count = 1;
index = 1;
for i = 1:7
    for j = 1:2
        %index = 1;

        subplot(7, 2, plot_count);
        plot_count = plot_count + 1;

        [chamber_temp, throat_temp] = getTemps(chamber_pressure(1,index), temp);
        index = index + 1;

        of_ratio = of_ratio';

        plot(of_ratio, chamber_temp, "-r");
        hold on
        plot(of_ratio, throat_temp, "-b");
        hold off;
        title("Pressure: " + chamber_pressure(1,index));
        xlabel("OF Ratio");
        ylabel("Temperature [K]");
        legend("Chamber Temp", "Throat Temp", "Location", "southeast");

     
    end
end




