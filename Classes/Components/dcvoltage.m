% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: DC voltage source

classdef dcvoltage
    properties
        Vdc; % Amplitude [V]
    end

    methods
        function obj = dcvoltage(varargin)
            % SINEVOLTAGE Class constructor.
            % Properties:
            % - Vdc: amplitude [V]
            p = inputParser;
            p.addParameter('Vdc',1);
            p.parse(varargin{:});
            obj.Vdc = p.Results.Vdc;
        end
    end
end