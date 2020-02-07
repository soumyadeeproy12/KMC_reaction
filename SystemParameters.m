classdef SystemParameters
    %SYSTEMPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        volumeOfReactor = 1.0*10E-4;
        viscosity = 0.00089;                      %0.0002816; (100C) 
		temp = 273.15;                              %0.0003142; (90C)
                                        			%0.0003540; (80C)
													%0.0008900; (25C)		
    end 
    methods
%         function value = get.temp(thisObj)
%             value = thisObj.temp;
%         end
          function value=getTemp(thisObj)
              value = thisObj.temp;
          end
          function value = getViscosity(thisObj)
              value = thisObj.viscosity;
          end
          function value = getVR(thisObj)
              value = thisObj.volumeOfReactor;
          end
    end
end

