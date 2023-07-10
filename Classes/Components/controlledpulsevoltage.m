% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Controlled pulse voltage source

classdef controlledpulsevoltage
    properties
        Vhigh; % High level voltage [V]
        Vlow; % Low level voltage [V]
        Vthresh; % Threshold voltage [V]
    end

    methods
        function obj = controlledpulsevoltage(varargin)
            % PULSEVOLTAGE Class constructor.
            % Properties:
            % - Vhigh: high level voltage [V]
            % - Vlow: low level voltage [V]
            % - Vthresh: threshold voltage [V]
            p = inputParser;
            p.addParameter('Vhigh',1);
            p.addParameter('Vlow',-1);
            p.addParameter('Vthresh',0);
            p.parse(varargin{:});
            obj.Vhigh = p.Results.Vhigh;
            obj.Vlow = p.Results.Vlow;
            obj.Vthresh = p.Results.Vthresh;
        end

        function v = voltage(obj,vc)
            % VOLTAGE Compute the source voltage as a function of the
            % control voltage.
            % Arg:
            % - t: time [s]
            % Return:
            % - vc: control voltage [V]
            v = obj.Vlow*ones(size(vc));
            v(vc < obj.Vthresh) = obj.Vhigh;
        end
    end
end