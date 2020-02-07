classdef CommonFunction 
    % This class contains all the common methods needed for simulation;
    properties (Access = private)
        vol, radious, dia;
    end     
    methods
        function thisObj = CommonFunction 
            clear vol radious dia;
        end
        function vol = getVolFromDia (thisObj, dia)		
            vol = (pi/6)*dia^3;		
        end
        
        function radious = getRadiousFromVol (thisObj, vol)		
            radious =  ((3/(4*pi))*vol)^(1/3);
        end
        
        function dia = getDiaFromVol (thisObj, vol) 		
            dia = ((vol*6)/pi)^(1/3);
        end
        
    end
end

