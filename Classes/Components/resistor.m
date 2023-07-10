% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Resistor

classdef resistor
    properties
        R; % Resistance [Ohm]
    end

    methods
        function obj = resistor(varargin)
            % RESISTOR Class constructor.
            % Property:
            % - R: resistance [Ohm]
            p = inputParser;
            p.addParameter('R',1e3);
            p.parse(varargin{:});
            obj.R = p.Results.R;
        end

        function i = current(obj,v)
            % CURRENT Compute the resistor current as a function of
            % voltage.
            % Arg:
            % - v: voltage [V]
            % Return:
            % - i: current [A]
            i = v/obj.R;
        end
    end
end