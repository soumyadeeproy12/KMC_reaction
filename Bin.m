classdef Bin <handle
    %BIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       noOfParticles, size;
    end
    
    methods
        function obj = Bin(nParticles, bSize)
           obj.noOfParticles = nParticles;
           obj.size = bSize;           
        end
        function setSize(thisObj, newSize)
            thisObj.size = newSize;
        end
        function addNoOfparticles (thisObj, number)
            thisObj.noOfParticles = thisObj.noOfParticles + number;
        end
        function removeNoOfParticles (thisObj, number)
            if ((thisObj.noOfParticles - number) >= 0)
                thisObj.noOfParticles = thisObj.noOfParticles - number;
            end
        end
        function value = getSize(thisObj)
            value = thisObj.size;
        end
        function pVolume = getPVolume(thisObj)
            pP = ParticleParameters;
            pVolume = pP.volumeOfMolecule*thisObj.size;
        end
        function pRadious = getPRadious(thisObj)
            cF = CommonFunction(); 
            volume = thisObj.getPVolume();
            pRadious = cF.getRadiousFromVol(volume);
        end
        function pDia = getPDia (thisObj)
            cF = CommonFunction(); 
            volume = thisObj.getPVolume();
            pDia = cF.getDiaFromVol(volume);
        end
        function nParticles = getNParticles (thisObj)
            nParticles = thisObj.noOfParticles;
        end
        function devideParticles (thisObj)
            thisObj.noOfParticles = thisObj.noOfParticles/2 ;
        end
        function mass = getBinMass(thisObj)
            mass = 0;
            if (thisObj.noOfParticles > 0)
                mass = thisObj.noOfParticles*thisObj.size;
            end
        end
    end
end

