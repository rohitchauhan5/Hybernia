days = 0:5:35;
conc = [515 420 370 250 135 120 60 20];
temp = [29,23,27,25,20,23,23,27];
[ax,b,p] = plotyy(days,temp,days,conc,'bar','plot');
title('Trend Chart for Concentration')
xlabel('Day')
ylabel(ax(1),'Temperature (^{o}C)')
ylabel(ax(2),'Concentration')
p.LineWidth = 3;
p.Color = [0,0.7,0.7];

A = magic(3); B = pascal(3);
C = cat(4, A, B);
