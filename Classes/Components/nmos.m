% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: NMOS transistor

classdef nmos
    properties (SetAccess = immutable)
        Knc; % Transconductance gain coefficients
        Vtn0c; % Zero bias forward threshold voltage coefficients
        gammac; % Body effect parameter coefficients
        lambdac; % Channel modulation factor coefficients
        phic; % Surface potential coefficients
        L = 4e-6; % Channel length [m]
        W = 5.1e-6; % Channel width [m]
    end

    properties
        n; % Aspect ratio factor
    end

    methods
        function obj = nmos(varargin)
            % NMOS Class constructor.
            % Arg:
            % - n: aspect ratio factor
            p = inputParser;
            p.addParameter('n',1);
            p.parse(varargin{:});
            obj.n = p.Results.n;
            load('nmos2.mat','Knc','Vtn0c','gammac','lambdac','phic');
            obj.Knc = Knc;
            obj.Vtn0c = Vtn0c;
            obj.gammac = gammac;
            obj.lambdac = lambdac;
            obj.phic = phic;
        end

        function iD = current(obj,vGS,vSB,vDS)
            % CURRENT Compute the NMOS current as a function of
            % applied voltages.
            % Arg:
            % - vGS: gate source voltage [V]
            % - vSB: source body voltage [V]
            % - vDS: drain source voltage [V]
            % Return:
            % - iD: drain current [A]
            vGS(vDS < 0) = vGS(vDS < 0)-vDS(vDS < 0); 
            if size(vGS,2) > 1
                vgs = vGS';
            else
                vgs = vGS;
            end
            vgs(vgs < 0.61) = 0.61;
            Kn = obj.n*(vgs.^(0:(length(obj.Knc)-1)))*obj.Knc;
            Vtn0 = (vgs.^(0:(length(obj.Vtn0c)-1)))*obj.Vtn0c;
            gamma = (vgs.^(0:(length(obj.gammac)-1)))*obj.gammac;
            lambda = (vgs.^(0:(length(obj.lambdac)-1)))*obj.lambdac;
            phi = (vgs.^(0:(length(obj.phic)-1)))*obj.phic;
            if size(vGS,2) > 1
                Kn = Kn';
                Vtn0 = Vtn0';
                gamma = gamma';
                lambda = lambda';
                phi = phi';
            end
            vSB(vDS < 0) = vSB(vDS < 0)+vDS(vDS < 0);
            vds = vDS;
            vds(vDS < 0) = -vds(vDS < 0);
            Vtn = Vtn0+gamma.*(sqrt(vSB+phi)-sqrt(phi));
            triode = vGS > Vtn & vds < vGS-Vtn;
            saturation = vGS > Vtn & vds >= vGS-Vtn;
            iD = zeros(size(vDS));
            idtriode = Kn.*(1+lambda.*vds).*((vGS-Vtn).*vds-(vds.^2)/2);
            idsat = (Kn/2).*(1+lambda.*vds).*(vGS-Vtn).^2;
            iD(triode) = idtriode(triode);
            iD(saturation) = idsat(saturation);
            iD(vDS < 0) = -iD(vDS < 0);
        end
    end
end