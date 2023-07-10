% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Pulse voltage source

classdef pulsevoltage
    properties
        Voff; % DC offset [V]
        V0; % Amplitude [V]
        f0; % Frequency [Hz]
        phi0; % Phase [rad]
    end

    methods
        function obj = pulsevoltage(varargin)
            % PULSEVOLTAGE Class constructor.
            % Properties:
            % - Voff: DC offset [V]
            % - V0: amplitude [V]
            % - f0: frequency [Hz]
            % - phi0: phase [rad]
            p = inputParser;
            p.addParameter('Voff',0);
            p.addParameter('V0',1);
            p.addParameter('f0',1e3);
            p.addParameter('phi0',0);
            p.parse(varargin{:});
            obj.Voff = p.Results.Voff;
            obj.V0 = p.Results.V0;
            obj.f0 = p.Results.f0;
            obj.phi0 = p.Results.phi0;
        end

        function v = voltage(obj,t)
            % VOLTAGE Compute the source voltage as a function of time.
            % Arg:
            % - t: time [s]
            % Return:
            % - v: voltage [V]
            v = obj.Voff+obj.V0./(1+exp(50*cos(4*pi*obj.f0*t-2* ...
                obj.phi0))).*(1-2./(1+exp(50*sin(2*pi*obj.f0*t- ...
                obj.phi0))));
        end
    end
end