clc;
clear all;
clear classes;
tic
matFilename = 'OutPut_MAT';
xlsFilename = 'OutPut_XlSX';
mkdir(matFilename);
mkdir(xlsFilename)
tMax = 60*5;                             % time in sec
iterations = 1; 
saveInterval = 1;                       % interval used to save data in time (min)
tCurrent = 0;
aD = AllDrops();
aD.Initialization();
dP = DropParameters();
saveData = 0;
leadingTime = 0;
laggingTime = 0;
while (tCurrent < tMax)
    u = rand;
    nDrops = aD.getNoOfDrops();
    %-----------------Random selection of drops for colasence-------------%
    if (nDrops > dP.minDrops)
        dropOne = randi([1 nDrops]);
        dropTwo = randi([1 nDrops]);
        while (dropOne == dropTwo) 
            dropTwo = randi([1 nDrops]);
        end
    end
    %-----------------------------end-------------------------------------%
    
    %---------------------coagulation ------------------------------------%
    dropCoagulation = randi([1 nDrops]);
    dC = aD.getDrop(dropCoagulation);
    nP = dC.getNoOfParticles();
    p1 = 1;
    p2 = 1;
    if (nP > 10)
        nP = floor(nP);
        p1 = randi(nP);
        p2 = randi(nP);
        while (p1 == p2) 
            p2 = randi(nP);
        end
    end
    
    totalFre = aD.getTotalFre(p1,p2,dropCoagulation);
    %-------------------------end-----------------------------------------%
    
    aD.doIvent(dropOne, dropTwo, p1, p2, dropCoagulation);  
    
    %-------------Interval quisence --------------------------------------%
    
    delT = (-log(1-u)/totalFre);    
    tCurrent = tCurrent + delT;
    %-----------------end-------------------------------------------------%
    
    aD.solveAllDGrowth(laggingTime, tCurrent);
    iterations = iterations + 1;
    laggingTime = tCurrent;
    if (rem(iterations,100) == 0) 
        disp(['mean particle size is   = ' , num2str(aD.getMeanParticleSize())]);
        disp(['No of iterations is     = ' , num2str(iterations)]);
        disp(['Current time            = ' , num2str(tCurrent)]);
        disp(['Current time step       = ' , num2str(delT)]);        
    end
    if ((tCurrent/60) >= saveData)
%       aD.saveAllDropsMat(matFilename, tCurrent);
        aD.saveAllDrops(xlsFilename, tCurrent);
        fileName = strcat(matFilename, '\Main_', num2str(tCurrent), '.mat');
        save(fileName);
        saveData=saveData+saveInterval;       
    end
end
%save('main.mat')
toc
