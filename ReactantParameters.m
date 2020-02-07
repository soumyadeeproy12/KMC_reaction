classdef ReactantParameters
    %REACTANTPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
	
    properties (Access = public)
       concAuCl, concReducingAgent;
    end
    
    methods(Access = public)
        function obj = ReactantParameters()
            obj.concAuCl = 0.254*10E-4;
            obj.concReducingAgent = 0.254*10E-3;
        end 
        function concAuCl = getConcAucl (thisObj)
            concAuCl = thisObj.concAuCl;
        end
        function concReducingAgent = getConcReducingAgent(thisObj)
            concReducingAgent = thisObj.concReducingAgent;
        end
    end
end

