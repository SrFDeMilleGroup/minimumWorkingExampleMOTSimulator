clear all;
close all;
simTrapDynamics =0; %change to 1 if you want to simulate particle trajectory, with random photon scatter, to get 'true' size and temperature
moleculeName = "SrF";
molecule = def_molecule(moleculeName);

colors = [[0.5;0.2;0.8],[0.4;0.7;0.8],[0.8;0.5;0.2],[0.5;0.7;0.2],...
    [0.6;0.6;0.3],[0.6;0.3;0.6],[0.3;0.6;0.6],[0.4;0.4;0.8],[0.4;0.8;0.4],[0.8;0.4;0.4],[0.5;0.5;0.5],[0.9;0.2;0.1],...
    [0.2;0.4;0.7],[0.4;0.4;0.0],[0.0;0.4;0.2],[0.2;0.0;0.4],[0.0;0.0;0.0],[0.8;0.8;0.8]];

dataFolder = 'saveData/SrFRedMOTNormalValuesbFieldSettingThreeDBGradGPerCM12.5ForceThreeDNumLasers5Date20240510_1636';

% pos = {'0.5','1.0','1.5','2.0','2.5','3.0'};
pos = {'0.5','1.5','3.0','4.5','6.0','7.5'};
for i=1:length(pos)
    posForPlot(i) = str2num(pos{i});
end
posInMM = posForPlot.*1e-3;
testFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{1},'MMSameDir.dat');
testData = readtable(testFile);

accelsInVelDirection = zeros(size(testData,1),length(pos));
accelsInPosDirection = zeros(size(testData,1),length(pos));
excitedPop = zeros(size(testData,1),length(pos));
for i=1:length(pos)
    currFile = strcat(dataFolder,'/forcevsSpeedDisplacement',pos{i},'MMSameDir.dat');
    currData = readtable(currFile);
    vels = currData.Speed;
    accelsInVelDirection(:,i) = currData.av;
    accelsInPosDirection(:,i) = currData.ar;
    excitedPop(:,i) = currData.PExc;
end

%reverse sign of accel for negative velocities
for i=1:length(vels)
    if vels(i)<0
        accelsInVelDirection(i,:) = accelsInVelDirection(i,:).*-1;
    end
end

%add opposite positions
accelsInVelDirectionFull = [-flipud(fliplr(accelsInVelDirection)),accelsInVelDirection];
excitedPopFull = [flipud(fliplr(excitedPop)),excitedPop];

sortedPos = sort([-posForPlot,posForPlot]);
oneAxisAccel = @(d,v) interp2(sortedPos,vels,accelsInVelDirectionFull./1e3,d,v);%in mm,ms units

%simulate capture with linear interpolation (spline acts up here for some
%reason.  Probably not a big deal to use linear)
diffEqVals = @(t,p) [p(2);oneAxisAccel(p(1),p(2))];
vsToTryForCapture = [1:1:30];
for i=1:length(vsToTryForCapture)
    currV = vsToTryForCapture(i);
    [ts2,ps2] = ode23(diffEqVals,[0 20],[min(sortedPos);currV]);
%     [ts2,ps2] = ode23(diffEqVals,[0 1e-1],[-0*1e-3;currV]);
    if isnan(ps2(end,1))
        break;
    end
end

%plot the highest value of v for which we have capture (TO DO)
[ts2,ps2] = ode23(diffEqVals,[0 20],[min(sortedPos);currV-1]);
figure(1);
plot(ps2(:,1),ps2(:,2));
title(strcat('particle trajectory for v_{Cap}=',num2str(currV-1),' m/s'))

%now plot a(x,v) heat map
oneAxisAccel = @(d,v) interp2(sortedPos,vels,accelsInVelDirectionFull./1e3,d,v,'spline');%in mm,ms units


%make heat map
vsForHeatMap = [min(vels):.1:max(vels)];
xsForHeatMap = [min(sortedPos):.1:max(sortedPos)];
for i=1:length(vsForHeatMap)
    for j=1:length(xsForHeatMap)
        heatMap(i,j) = oneAxisAccel(xsForHeatMap(j),vsForHeatMap(i));
    end
end
figure(2);
imagesc(xsForHeatMap,vsForHeatMap,heatMap)
colorbar
xlabel('x(mm)');
ylabel('v (m/s)')
xlim([min(sortedPos) max(sortedPos)])
h=colorbar;
h.Title.String = "a (mm/ms^{2})"

%make x,v plots from matrix
vMaxToInt = 4;
xMaxToInt = 7;
[~,minCol] = min(abs(xsForHeatMap+xMaxToInt));
[~,maxCol] = min(abs(xsForHeatMap-xMaxToInt));
[~,minRow] = min(abs(vsForHeatMap+vMaxToInt));
[~,maxRow] = min(abs(vsForHeatMap-vMaxToInt));
accVsPosForPlot = trapz(vsForHeatMap(minRow:maxRow),heatMap(minRow:maxRow,:),1)./(2*vMaxToInt);
accVsVelForPlot = trapz(xsForHeatMap(minCol:maxCol),heatMap(:,minCol:maxCol),2)./(2*xMaxToInt);
figure(3);
hold all;
plot(xsForHeatMap,accVsPosForPlot,'Linewidth',2);
xlabel('x(mm)');
ylabel('a(mm/ms^{2})')
title('Acceleration vs Position')
figure(4);
hold all;
plot(vsForHeatMap,accVsVelForPlot,'LineWidth',2);
xlabel('v(m/s)');
ylabel('a(mm/ms^{2})')
title('Acceleration vs Velocity')

[~,minCol] = min(abs(sortedPos+xMaxToInt));
[~,maxCol] = min(abs(sortedPos-xMaxToInt));
[~,minRow] = min(abs(vels+vMaxToInt));
[~,maxRow] = min(abs(vels-vMaxToInt));
meanExcPop = mean(mean(excitedPopFull(minRow:maxRow,minCol:maxCol)));

%sim trap dynamics
maxTime=40;
gam = molecule.gam;
kA = molecule.kA;
mass = molecule.mass;
hbar = 1.05e-34;
kb = 1.38e-23;
if simTrapDynamics==1
    %     scatterRateData = dlmread(strcat(folder,'forceVsVCenter.dat'));
    %     scatterRate = scatterRateData(11,3).*gamSrF;
    scatterRate = meanExcPop*gam*1e-3;%in 1/ms
    tKick = 1/scatterRate;
    v(1) = currV-3;%mm/ms
    r(1) = -7;%mm
    vKick = hbar*kA/mass;
    for i=1:round(maxTime/tKick)
        randPhi1 = 2*pi*rand;
        randPhi2 = 2*pi*rand;
        v(i+1) = v(i)+vKick*(cos(randPhi1)+cos(randPhi2))+oneAxisAccel(r(i),v(i))*tKick;
        r(i+1) = r(i)+v(i)*tKick;
    end
    simTimes = 0:(tKick):maxTime;
    startTime = maxTime/2;
    startInd = maxTime/2/tKick;
    endInd= i;
    vSq=mean(v(startInd:endInd).^2);
    sigma=sqrt(mean(r(startInd:endInd).^2));
    temp = vSq*mass/kb;
end

    