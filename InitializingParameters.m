classdef InitializingParameters
    %INITIA Summary of this class goes here
    %   Detailed explanation goes here
    properties
       nDrops, reducingDrops, reactingDrops;
    end
    
    methods
        function thisObj = InitializingParameters()
            thisObj.nDrops = 100;
            thisObj.reducingDrops = thisObj.nDrops * 0.5;
            thisObj.reactingDrops = thisObj.nDrops - thisObj.reducingDrops;
        end        
    end
end

