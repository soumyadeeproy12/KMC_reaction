classdef AllDrops < handle
        
    properties
        allDrops;
    end
    
    methods
        function thisObj = AllDrops()
           thisObj.allDrops = Drop.empty;
        end
        
        function thisObj = Initialization (thisObj)
            iP = InitializingParameters ();         
            for dCount = 1:1:iP.reactingDrops
                thisObj.creatReactingDrop();                
            end
            
            for dCount=1:1:iP.reducingDrops                
                thisObj.creatReducingDrop();
            end
        end
        function coalescenceFre = getCoalescenceFre (thisObj)
            coalescenceFre = 0; 
            dP = DropParameters();
            if (thisObj.getNoOfDrops() > dP.minDrops)
                dP = DropParameters();
                coalescenceFre = 0.5*dP.coalescenceEff*thisObj.getNoOfDrops()^2;
            end
        end
        function breakageFre = getBreakageFre (thisObj)
            dP = DropParameters();
            breakageFre = dP.breakageEff*thisObj.getNoOfDrops();
        end
        function totalFre = getTotalFre (thisObj, p1, p2, dCount)
            totalFre = thisObj.getCoalescenceFre() + thisObj.getBreakageFre() + thisObj.getTotalNucleationFre() + thisObj.getTotalCoagulationFre(p1, p2, dCount);            
        end
        function doIvent(thisObj, dropOne, dropTwo, p1, p2, cCount)
            
            u = rand;
            status = 0;
            %-------p-------All frequency----------------------------------%
            coalescenceFre = thisObj.getCoalescenceFre(); 
            breakageFre = thisObj.getBreakageFre ();
            coaFre = thisObj.getTotalCoagulationFre(p1, p2, cCount);
            nueFre = thisObj.getTotalNucleationFre();            
            totalFre = thisObj.getTotalFre(p1, p2, cCount);
            %--------------------end--------------------------------------%
            
            %--------------All Probabilities------------------------------%
            cProb = coalescenceFre/totalFre;
            bProb = (coalescenceFre+breakageFre)/totalFre;
            coaProb = (coalescenceFre+breakageFre+coaFre)/totalFre;
            nueProb = (coalescenceFre+breakageFre+coaFre+nueFre)/totalFre;
            %--------------------end--------------------------------------%
            if (u > 0 && u <= cProb)
                thisObj.mergeDrop(dropOne, dropTwo);
                status = 1;
            elseif (u > cProb && u <= bProb)
                index = randi([1 thisObj.getNoOfDrops]);
                thisObj.breakDrop(index);
                status = 1;
            elseif (u > bProb && u <= coaProb)
                thisObj.doCoagulation (p1, p2, cCount)
                status = 1;
            elseif (u > coaProb && u <= nueProb)
                thisObj.doNucleation (p1, p2, cCount, u)
                status = 1;
            end
            if (status == 0)
                disp(['Random number is ', num2str(u)]);
                disp(['Nucleation prob ', num2str(nueProb)]);
                disp(['Coagulation prob ', num2str(coaProb)]);
                disp(['Breakage prob ', num2str(bProb)]);
                disp(['Colesence prob ', num2str(cProb)]);
                disp(['number of particles in system', thisObj.getTotalNoOfParticles()]);
                disp('no Event occured');
            end
        end
        function thisObj = creatReactingDrop (thisObj)
            cF = CommonFunction();
            rP = ReactantParameters();
            dP = DropParameters();
            vol  = cF.getVolFromDia(dP.initialDropSize);
            conc = rP.concAuCl;
            moles = conc*vol;
            d = Drop(moles, 0, 0, dP.initialDropSize);  
            d.reactantDrop= 1;
            thisObj.allDrops = [thisObj.allDrops, d];
        end
        function thisObj = creatReducingDrop (thisObj)
            cF = CommonFunction();
            rP = ReactantParameters();
            dP = DropParameters();
            vol  = cF.getVolFromDia(dP.initialDropSize);
            conc = rP.concReducingAgent;
            moles = conc*vol;
            d = Drop(0, moles, 0, dP.initialDropSize); 
            d.reducingDrop = 1;
            thisObj.allDrops = [thisObj.allDrops, d];
        end
        
        %----------------Diffusional growth for all Drops ----------------%
        function solveAllDGrowth (thisObj, tInitial, tFinal)
            nDrops = thisObj.getNoOfDrops();
            for dCount=1:1:nDrops
                d = thisObj.getDrop(dCount);
                if (d.productDrop == 1)
                    d.solveDiffusionalGrowth(tInitial, tFinal);
                end
            end        
        end
        %------------------------end--------------------------------------%
        
        function thisObj = breakDrop (thisObj, bCount)
            bins1 = [];
            cF = CommonFunction();
            d = thisObj.getDrop(bCount);
            volume = d.getVolume();
            dia = cF.getDiaFromVol(volume/2);
            d1 = Drop(d.molesOfReactant/2 , d.molesOfReducingAgent/2, d.molesOfProduct/2, dia);
            d1.reactantDrop = d.reactantDrop;
            d1.reducingDrop = d.reducingDrop;
            d1.productDrop  = d.productDrop;
            if (d1.productDrop == 1)
%                 %-------------uniform distribution in terms of mass-------%
%                 d.devideBins();
%                 bins1 = d.getAllBins() ;
%                 d1.setAllBins(bins1);
%                 %--------------------------end----------------------------%
%                 

            end
            d2 = Drop(d.molesOfReactant/2 , d.molesOfReducingAgent/2, d.molesOfProduct/2, dia);
            d2.reactantDrop = d.reactantDrop;
            d2.reducingDrop = d.reducingDrop;
            d2.productDrop  = d.productDrop;
            if (d1.productDrop == 1 && d2.productDrop == 1)
                %----------non uniform distribution in terms of mass---%
                d.distributionByMass (d1, d2);                
                %-----------------------end----------------------------%                
                
                
%                 %-------------uniform distribution in terms of mass-------%
%                 d2.setAllBins(bins1);
%                 %--------------------------end----------------------------%

            end
            thisObj.removeDrop(d);
            thisObj.addDrop(d1);
            thisObj.addDrop(d2);
        end 
        
        function thisObj = mergeDrop (thisObj, m1Count, m2Count)
            d1 = thisObj.getDrop(m1Count);
            d2 = thisObj.getDrop(m2Count);
            cF = CommonFunction();		  
		    vol = d1.volume + d2.volume;            
		    dia = cF.getDiaFromVol(vol);
            
            if (d1.reactantDrop == 1 && d2.reducingDrop == 1)
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                if (molesReactant <= molesOfReducingAgent)
                    molesOfReducingAgent = molesOfReducingAgent - molesReactant;
                    molesProduct = molesProduct + molesReactant;
                    molesReactant = 0.0;
                end
                if (molesReactant > molesOfReducingAgent)
                    molesReactant = molesReactant - molesOfReducingAgent;
                    molesProduct = molesProduct + molesOfReducingAgent;
                    molesOfReducingAgent = 0.0;
                end                
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = 0;	
                d.reducingDrop = 0;
                d.productDrop = 1;               
                thisObj.addDrop(d);                
            elseif (d1.reducingDrop== 1 && d2.reactantDrop == 1)
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                if (molesReactant <= molesOfReducingAgent)
                    molesOfReducingAgent = molesOfReducingAgent - molesReactant;
                    molesProduct = molesProduct + molesReactant;
                    molesReactant = 0.0;
                end
                if (molesReactant > molesOfReducingAgent)
                    molesReactant = molesReactant - molesOfReducingAgent;
                    molesProduct = molesProduct + molesOfReducingAgent;
                    molesOfReducingAgent = 0.0;
                end                
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = 0;	
                d.reducingDrop = 0;
                d.productDrop = 1;               
                thisObj.addDrop(d);
            elseif (d1.reactantDrop == 1 && d2.productDrop == 1)
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                if (molesReactant <= molesOfReducingAgent)
                    molesOfReducingAgent = molesOfReducingAgent - molesReactant;
                    molesProduct = molesProduct + molesReactant;
                    molesReactant = 0.0;
                end
                if (molesReactant > molesOfReducingAgent)
                    molesReactant = molesReactant - molesOfReducingAgent;
                    molesProduct = molesProduct + molesOfReducingAgent;
                    molesOfReducingAgent = 0.0;
                end                
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = 0;	
                d.reducingDrop = 0;
                d.productDrop = 1;               
                thisObj.addDrop(d);
            elseif (d1.productDrop == 1 && d2.reactantDrop == 1)
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                if (molesReactant <= molesOfReducingAgent)
                    molesOfReducingAgent = molesOfReducingAgent - molesReactant;
                    molesProduct = molesProduct + molesReactant;
                    molesReactant = 0.0;
                end
                if (molesReactant > molesOfReducingAgent)
                    molesReactant = molesReactant - molesOfReducingAgent;
                    molesProduct = molesProduct + molesOfReducingAgent;
                    molesOfReducingAgent = 0.0;
                end                
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = 0;	
                d.reducingDrop = 0;
                d.productDrop = 1;               
                thisObj.addDrop(d);
            elseif (d1.productDrop == 1 && d2.productDrop == 1)
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                if (molesReactant <= molesOfReducingAgent)
                    molesOfReducingAgent = molesOfReducingAgent - molesReactant;
                    molesProduct = molesProduct + molesReactant;
                    molesReactant = 0.0;
                end
                if (molesReactant > molesOfReducingAgent)
                    molesReactant = molesReactant - molesOfReducingAgent;
                    molesProduct = molesProduct + molesOfReducingAgent;
                    molesOfReducingAgent = 0.0;
                end                
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = 0;	
                d.reducingDrop = 0;
                d.productDrop = 1;               
                thisObj.addDrop(d);
            elseif (d1.productDrop == 1 && d2.reducingDrop == 1)
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                if (molesReactant <= molesOfReducingAgent)
                    molesOfReducingAgent = molesOfReducingAgent - molesReactant;
                    molesProduct = molesProduct + molesReactant;
                    molesReactant = 0.0;
                end
                if (molesReactant > molesOfReducingAgent)
                    molesReactant = molesReactant - molesOfReducingAgent;
                    molesProduct = molesProduct + molesOfReducingAgent;
                    molesOfReducingAgent = 0.0;
                end                
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = 0;	
                d.reducingDrop = 0;
                d.productDrop = 1;               
                thisObj.addDrop(d);
            elseif (d1.reducingDrop == 1 && d2.productDrop == 1)
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                if (molesReactant <= molesOfReducingAgent)
                    molesOfReducingAgent = molesOfReducingAgent - molesReactant;
                    molesProduct = molesProduct + molesReactant;
                    molesReactant = 0.0;
                end
                if (molesReactant > molesOfReducingAgent)
                    molesReactant = molesReactant - molesOfReducingAgent;
                    molesProduct = molesProduct + molesOfReducingAgent;
                    molesOfReducingAgent = 0.0;
                end                
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = 0;	
                d.reducingDrop = 0;
                d.productDrop = 1;               
                thisObj.addDrop(d);
            else 
                molesReactant = d1.molesOfReactant + d2.molesOfReactant;
                molesOfReducingAgent = d1.molesOfReducingAgent + d2.molesOfReducingAgent;
                molesProduct = d1.molesOfProduct + d2.molesOfProduct;
                d = Drop (molesReactant, molesOfReducingAgent, molesProduct, dia);
                d.reactantDrop = d1.reactantDrop;	
                d.reducingDrop = d1.reducingDrop;
                d.productDrop = d1.productDrop;               
                thisObj.addDrop(d);
            end
            bins1 = d1.getAllBins();
            bins2 = d2.getAllBins();
            d.mergeBins(bins1, bins2)
            thisObj.removeDrop(d1);
            thisObj.removeDrop(d2);
        end
        
        function meanParticleSize = getMeanParticleSize (thisObj)
            meanParticleSize = 0;
            totalSize = 0;
            avgDrops = 0;
            noOfDrops = thisObj.getNoOfDrops();
            for dCount = 1:1:noOfDrops
                d = thisObj.getDrop(dCount);
                if (d.productDrop == 1)
                    totalSize = totalSize + d.getMeanSize();
                    avgDrops = avgDrops + 1;
                end
            end
            if (avgDrops ~= 0)
                meanParticleSize = totalSize/avgDrops;
            end                
        end
        
        function obj = getDrop(thisObj, gCount) 
            if (gCount > 0)
                obj = thisObj.allDrops(gCount);
            end
        end
        function index = getIndex(thisObj, d)
            index = find(thisObj.allDrops == d);
        end        
        function thisObj = removeDrop(thisObj, d)           
            thisObj.allDrops(thisObj.allDrops == d) = [];
        end
        function thisObj = addDrop(thisObj, d) 
            thisObj.allDrops = [thisObj.allDrops, d];
        end
        function number = getNoOfDrops(thisObj) 
            number = length(thisObj.allDrops);
        end
        function printallDrop(thisObj)
            numberOfDrops = length(thisObj.allDrops);
            for dCount = 1:numberOfDrops 
                d = thisObj.getDrop(dCount);
                disp(['molesOfReactant of object is       = ' , num2str(d.molesOfReactant)]);
                disp(['molesOfProduct of object is        = ', num2str(d.molesOfProduct)]);
                disp(['molesOfReducingAgent of object is  = ', num2str(d.molesOfReducingAgent)]);
                disp([' ......',num2str(dCount)]);
            end
        end 
        function saveAllDropsMat (thisObj, mFileName, time)
            fileName = strcat(mFileName, '\allDrops_', num2str(time), '.mat');
            save(fileName);
        end
        function saveAllDrops (thisObj, xFileName, time)            
            fileName = strcat(xFileName,'\saveAllDrops_',num2str(time),'.xlsx');       
            header={'MolesOfReactant' 'MolesOfReducing' 'MolesOfProduct' };           
            nDrops = thisObj.getNoOfDrops();         
            for nCount =1:1:nDrops
                d = thisObj.getDrop(nCount);
                saveAllDrops(nCount,1) = d.molesOfReactant;
                saveAllDrops(nCount,2) = d.molesOfReducingAgent;
                saveAllDrops(nCount,3) = d.molesOfProduct;
                saveAllDrops(nCount,4) = d.dia;
                saveAllDrops(nCount,5) = d.volume;                
                saveAllDrops(nCount,10)= d.reactantDrop;
                saveAllDrops(nCount,11)= d.reducingDrop;
                saveAllDrops(nCount,12)= d.productDrop;
                saveAllDrops(nCount,13)= d.getNoOfParticles();
                saveAllDrops(nCount,14)= d.getMeanSize();
            end
             xlswrite(fileName,saveAllDrops)
             clear saveAllDrops;
        end   
        %-------------------------After KMC model-------------------------%
        function totalNucleationFre = getTotalNucleationFre (thisObj)
            totalNucleationFre = 0;
            nDrops = thisObj.getNoOfDrops();
            for nCount = 1:1:nDrops
                d = thisObj.getDrop(nCount);
                totalNucleationFre = totalNucleationFre + d.getNucleationRate();
            end
            %totalNucleationFre = 0; %% delete after debugging
        end
        function totalNoOfParticles = getTotalNoOfParticles(thisObj)
            totalNoOfParticles = 0;
            nDrops = thisObj.getNoOfDrops();
            for nCount = 1:1:nDrops
                d = thisObj.getDrop(nCount);
                totalNoOfParticles = totalNoOfParticles + d.getNoOfParticles();
            end            
        end
        function totalCoagulationFre = getTotalCoagulationFre(thisObj, p1, p2, dCount)
            pP = ParticleParameters();
            totalCoagulationFre = 0;
            d=thisObj.getDrop(dCount);
            N = d.getNoOfParticles();
            if (N > pP.minParticles)
                c = Constants();
                sP = SystemParameters ();
                pP = ParticleParameters();
                b1 = d.getBinOfThisParticle(p1);                
                b2 = d.getBinOfThisParticle(p2);                
                r1 = b1.getPRadious();
                r2 = b2.getPRadious();                
                Qp = (2*c.kb*sP.getTemp()/3*sP.getViscosity())*(2+(r1/r2)+(r2/r1));
                totalCoagulationFre = 0.5*pP.beta*(Qp*N)*(sP.getVR()*N); 
            end
            totalCoagulationFre = 0; %% delete after debugging
        end
        function doCoagulation (thisObj, pOne, pTwo, dCount)
            d=thisObj.getDrop(dCount);
            b1 = d.getBinOfThisParticle(pOne);
            b2 = d.getBinOfThisParticle(pTwo);
            size = b1.getSize()+ b2.getSize();
            d.addThisBin (1, size);
            b1.removeNoOfParticles(1);
            b2.removeNoOfParticles(1);
        end
        function doNucleation (thisObj, p1, p2, cDrop, rNumber)
            %disp('Nucleation called ');
            pP = ParticleParameters ();
            c = Constants();
            changeInMole = pP.nucleiSize/c.NA;
            nCount = thisObj.getDropForNucleation(p1, p2, cDrop, rNumber);
           % disp(['Nucleation drops is  ', num2str(nCount)]);
            d = thisObj.getDrop(nCount);
            d.addThisBin(1, pP.nucleiSize);
            d.reduceMolesOfProduct(changeInMole);            
        end
        function nueDrop = getDropForNucleation (thisObj, p1, p2, cDrop, rNumber)
            p1Prob = 0;
            p2Prob = 0;
            nueDrop = 0;
            totalFre = thisObj.getTotalFre(p1, p2, cDrop);
            freExceptNucleation = totalFre - thisObj.getTotalNucleationFre();
            nDrops = thisObj.getNoOfDrops();
            %disp(['totalnumber of drops  is ', num2str(nDrops)]);
            for  dCount = 1:1:nDrops
                d1 = thisObj.getDrop(dCount);               
                p1Prob =  freExceptNucleation / totalFre;
                p2Prob = (freExceptNucleation + d1.getNucleationRate()) / totalFre;
%                  disp(['p1 probability is  ', num2str(p1Prob)]);
%                   disp(['p2 probability is  ', num2str(p2Prob)]);
                if (rNumber >= p1Prob && rNumber < p2Prob)
                   % disp(['Drops selected for nucleation ', num2str(dCount)]);
                    nueDrop = dCount;
                end
                freExceptNucleation = freExceptNucleation + d1.getNucleationRate();
            end           
        end
    end
end

