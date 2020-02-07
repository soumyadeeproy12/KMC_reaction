classdef DropParameters
    %DROPPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here		
    properties
        coalescenceEff, breakageEff, initialDropSize, minDrops;
    end    
    methods
        function thisObj = DropParameters()
            iP = InitializingParameters;
            thisObj.coalescenceEff = 0.001;
            thisObj.breakageEff =  0.5*iP.nDrops()*thisObj.coalescenceEff;		
            thisObj.initialDropSize = 2.0E-5;
            thisObj.minDrops = 10;
        end
    end    
end

