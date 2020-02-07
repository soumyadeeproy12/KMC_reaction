classdef ParticleParameters
    %PARTICLEPARAMETERS Summary of this class goes here    
    properties (Access = public)
        volumeOfMolecule, beta, solubility, diffusivity, nucleiSize, minParticles;
    end
    methods
        function thisObj = ParticleParameters ()
            sP = SystemParameters();
            thisObj.volumeOfMolecule = 1.695*10^-29;
            thisObj.beta = 2.0*10^-4;
            thisObj.solubility = 10^-6.0;
            thisObj.diffusivity = ( (7.4*10^-8*((2.6*18.01528)^0.5)*sP.temp) / (0.89*(10.20552^0.6)))*10E-4;
            thisObj.nucleiSize = 2;
            thisObj.minParticles = 10;
        end
    end
end

