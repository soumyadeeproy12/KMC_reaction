classdef Drop <handle
   properties (Access = public)
      molesOfReactant, molesOfProduct, molesOfReducingAgent, dia;
      volume;      
      reactantDrop = 0;  reducingDrop = 0; productDrop = 0;
      allBins;
   end
   
   methods
       function thisObj = Drop(mReactant, mReducingAgent, mProduct,  iDia)
           %if nargin == 8
                thisObj.molesOfReactant = mReactant;
                thisObj.molesOfReducingAgent = mReducingAgent;
                thisObj.molesOfProduct = mProduct;                
                thisObj.dia = iDia;
                thisObj.volume = (pi/6.0)*thisObj.dia^3;                
                thisObj.allBins= Bin.empty(1,0);
           %end
       end        
        function growthRate =  getGrowthRate (thisObj)
            pP = ParticleParameters ();
            growthConstant = 0.01E-8;
            lamda = thisObj.getConc() / pP.solubility;
            growthRate = growthConstant*(lamda^0.5-1.0)*(pP.solubility^0.5);		
        end
        function dDdt = diffusionalGrowth(thisObj)
           dDdt =  thisObj.getGrowthRate();
        end
        function solveDiffusionalGrowth (thisObj, tInitial, tFinal)
            pP = ParticleParameters();
            if(thisObj.getConc() > pP.solubility) 
            cF = CommonFunction();            
            c = Constants(); 
            nBins = thisObj.getNoOfBins();
            for bCount = 1:1:nBins
                b = thisObj.getThisBin(bCount);                
                oldSize = b.getSize();
                oldVol = oldSize*pP.volumeOfMolecule;
                oldDia = cF.getDiaFromVol(oldVol);
                tspan = [tInitial tFinal]; 
                [t,N] = ode45(@(t,N) thisObj.diffusionalGrowth(), tspan, oldDia);
                newDia = N(end);
                newVol = cF.getVolFromDia(newDia);
                newSize = newVol/pP.volumeOfMolecule;
                changeMol = ((newSize-oldSize)*b.getNParticles())/c.NA;
                cMol = thisObj.molesOfProduct - changeMol;
                cConc = thisObj.calculateConc(cMol);
                if (cConc > pP.solubility) 
                    b.setSize(newSize);
                    thisObj.reduceMolesOfProduct(changeMol);
                end
                 if (cConc <= pP.solubility)
                     disp(['Current conentration is ', num2str(thisObj.getConc())]);
                     disp(['After difussion growth conentration is ', num2str(cConc)]);
                     disp('Not sufficient material for gowth');
                    return;
                end
                if (thisObj.getConc() <= pP.solubility)
                    disp('Not sufficient material for gowth');
                    return;
                end
            end
            end
        end
        
        function nucleationRate = getNucleationRate (thisObj)            
            sP = SystemParameters ();
			c =  Constants ();
			pP = ParticleParameters ();
			lamda = thisObj.getConc()/pP.solubility();            
            if (lamda >= 1.0)	                
                expNumeratorTerm = -16.0*pi*(c.sigma^3)*(pP.volumeOfMolecule^2);
                expDenominatorTerm = 3.0*((c.kb*sP.temp)^3)*(log(lamda)^2.0);		
                nucleationRate = thisObj.volume*c.A*exp(expNumeratorTerm/expDenominatorTerm);
            else                
				nucleationRate = 0.0;
            end
        end
        
      function reduceMolesOfReactant(thisObj, molesReactantToBeReduce)
         thisObj.molesOfReactant = thisObj.molesOfReactant - molesReactantToBeReduce;
      end
      function reduceMolesOfProduct(thisObj, molesProductToBeReduce)
         nextMoles = thisObj.molesOfProduct - molesProductToBeReduce;
         if (nextMoles >= 0)
            thisObj.molesOfProduct = thisObj.molesOfProduct - molesProductToBeReduce;
         end
      end
      function reducingMolesOfReducingAgent(thiObj, molesReducingAgentToBeReduce) 
          thiObj.molesOfReducingAgent = thiObj.molesOfReducingAgent - molesReducingAgentToBeReduce;
      end      
      function conc = getConc (thisObj)
          conc = 0;
          if (thisObj.molesOfProduct > 0)
            conc = thisObj.molesOfProduct/thisObj.volume;
          end
      end
      function cConc = calculateConc(thisObj, moles)
          cConc = moles/thisObj.volume;
      end
      function vol = getVolume (thisObj)
          vol = thisObj.volume;
      end
      function diameter = getDia(thisObj)
          diameter = thisObj.dia;
      end
      function dAllBins = getAllBins(thisObj)
          dAllBins = thisObj.allBins;
      end
      function setAllBins(thisObj, sAllBins)
          thisObj.allBins = sAllBins;
      end
      function  obj = getThisBin (thisObj, bCount) 
          obj = [];
          if (bCount>0 && bCount <=thisObj.getNoOfBins())
            obj = thisObj.allBins(bCount);
          end
      end
      function obj = getThisSizeBin (thisObj, size)
          obj = [];
          nBins = thisObj.getNoOfBins();
          for nCount=1:1:nBins
              b = thisObj.getThisBin (nCount);
              if (round(size) == round(b.size))
                  obj = b;
              end
          end
      end
      function number = getNoOfBins (thisObj) 
          number = length(thisObj.allBins);
      end
      function devideBins (thisObj)
          nBins = thisObj.getNoOfBins();
          for nCount=1:1:nBins
              b = thisObj.getThisBin(nCount);
              b.devideParticles();
          end
      end
      function distributionByNumber (thisObj, drop1, drop2)          
          tempCount = 0;
          nP = thisObj.getNoOfParticles(); %Plan is to devide particles in bin non uniformly%
          nP = round(nP/2);         
          while (tempCount < nP)
              index = randi([1 thisObj.getNoOfParticles()]);
              b = thisObj.getBinOfThisParticle(index);
              size = b.getSize();
              b.removeNoOfParticles(1);
              if (b.getNParticles() <= 0)
                  b.removeThisBin(b);
              end
              drop1.addThisBin(1, size);
              tempCount = tempCount + 1;              
          end
          drop2Bins = thisObj.getAllBins();
          drop2.setAllBins(drop2Bins);
      end
      function distributionByMass (thisObj, drop1, drop2)  
          tempMass = 0;
          mP = thisObj.getTotalMass(); %Distribution of uniform mass between drops%
          mP = round(mP/2);
          while (tempMass < mP)
              index = randi([1 thisObj.getNoOfParticles()]);
              b = thisObj.getBinOfThisParticle(index);
              size = b.getSize();
              b.removeNoOfParticles(1);
              if (b.getNParticles() <= 0)
                  thisObj.removeThisBin(b);
              end
              drop1.addThisBin(1, size);
              tempMass = tempMass + size;
          end
          drop2Bins = thisObj.getAllBins();
          drop2.setAllBins(drop2Bins);
      end

      function mergeBins (thisObj, binArray1 ,binArray2)         
          tempArray = [binArray1 binArray2];
          nBins = length(tempArray);
          for bCount=1:1:nBins
              b = tempArray(bCount);
              particles = b.getNParticles();
              size = b.getSize();
              thisObj.addThisBin (particles, size);              
          end
      end
      function number = getNoOfParticles (thisObj)
          number = 0;
          nBins = thisObj.getNoOfBins();
          for nCount=1:1:nBins 
              b = thisObj.getThisBin(nCount);
              number = number + b.noOfParticles;
          end          
      end
      function pVol = getTotalPVol (thisObj)
          pVol = 0;
          nBins = thisObj.getNoOfBins();
          for bCount=1:1:nBins 
              b = thisObj.getThisBin(bCount);
              tempVol = b.getPVolume()*b.getNParticles();
              pVol = pVol + tempVol;
              tempVol = 0;
          end
      end  
      function addThisBin (thisObj, particles, size)
          b = thisObj.getThisSizeBin(round(size));
          if (isempty(b)== 1)
              newBin = Bin(particles, size);
              thisObj.allBins = [thisObj.allBins, newBin];
          else
              b.addNoOfparticles(particles);
          end          
      end
      
      function removeThisBin(thisObj, b)
          thisObj.allBins(thisObj.allBins == b) = [];
      end
      function obj = getBinOfThisParticle (thisObj, particleIndex)     
          obj = [];
          number = 0;          
          nBins = thisObj.getNoOfBins();
          for nCount=1:1:nBins
              oldNumber = number;
              b = thisObj.getThisBin(nCount);
              number = number + b.noOfParticles;
              if (particleIndex > oldNumber && particleIndex <= number)
                  if (isempty(b) == 0)
                    obj = b;
                    return;
                  end
              end
          end
      end
      function meanVol = getMeanVol(thisObj)
          meanVol = 0;
          nP = thisObj.getNoOfParticles();
          if (nP > 0)
            vP = thisObj.getTotalPVol();
            meanVol = vP/nP;
          end
      end
      function meanSize = getMeanSize(thisObj)
          meanSize = 0;
          if (thisObj.productDrop == 1)              
                cF = CommonFunction();
                meanVol = thisObj.getMeanVol();
                meanSize = cF.getDiaFromVol(meanVol);              
          end
      end
      function totalMass = getTotalMass(thisObj)
          totalMass = 0;
          bMax = thisObj.getNoOfBins();
          for bCount = 1:1:bMax
              b = thisObj.getThisBin(bCount);
              totalMass = totalMass + b.getBinMass();
          end
      end
      function disPlayProperties (thisObj) 
          disp(['molesOfReactant of object is       = ' , num2str(thisObj.molesOfReactant)]);
          disp(['molesOfProduct of object is        = ', num2str(thisObj.molesOfProduct)]);
          disp(['molesOfReducingAgent of object is  = ', num2str(thisObj.molesOfReducingAgent)]);
          disp(' ......');
      end
   end
end