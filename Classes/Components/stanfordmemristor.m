% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Stanford model memristor

classdef stanfordmemristor
    properties (SetAccess = immutable)
        q = 1.6e-19; % Elementary charge [C]
        k = 1.38e-23; % Boltzmann constant [J/K]
        I0 = 61.4e-6; % Quasi static-law current coefficient [A]
        V0 = 0.43; % Quasi static-law voltage coefficient [V]
        g0 = 0.275e-9; % Quasi static-law gap coefficient [m]
        v0 = 150; % Gap velocity coefficient [m/s]
        Eag = 1.501; % Activation energy for vacancy generation [eV]
        Ear = 1.5; % Activation energy for vacancy recombination [eV]
        gamma0 = 16.5; % Field enhancement factor coefficient
        beta = 1.25; % Field enhancement gap coefficient
        g1 = 1e-9; % Dynamic law gap coefficient [m]
        a0 = 0.25e-9; % Atomic hopping distance [m]
        L = 5e-9; % Oxide thickness [m]
        gmin = 0.1e-9; % Minimum gap [m]
        gmax = 1.7e-9; % Maximum gap [m]
        T0 = 298; % Room temperature [K]
        Cth = 0.318e-15; % Thermal capacitance [J/K]
        tauth = 0.23e-9; % Thermal time constant [s]
    end

    properties
        gic; % Gap initial value [m]
        Tic; % Temperature initial value [K]
    end

    methods
        function obj = stanfordmemristor(varargin)
            % STANFORDMEMRISTOR Class constructor.
            % Properties:
            % - gic: gap initial value [m]
            % - Tic: temperature initial value [K]
            p = inputParser;
            p.addParameter('gic',obj.gmax);
            p.addParameter('Tic',obj.T0);
            p.parse(varargin{:});
            obj.gic = p.Results.gic;
            obj.Tic = p.Results.Tic;
        end

        function i = current(obj,v,g)
            % CURRENT Compute the memristor current as a function of
            % voltage and gap.
            % Args:
            % - v: voltage [V]
            % - g: gap [m]
            % Return:
            % - i: current [A]
            i = obj.I0*exp(-g/obj.g0).*sinh(v/obj.V0);
        end

        function dgdt = gaptimederivative(obj,v,g,T,gregion)
            % GAPTIMEDERIVATIVE Compute the memristor gap time derivative
            % as a function of voltage, gap and temperature.
            % Args:
            % - v: voltage [V]
            % - g: gap [m]
            % - T: temperature [K]
            % - gregion: gap region
            % Return:
            % - dgdt: gap time derivative [m/s]
            gamma = obj.gamma0-obj.beta*(g/obj.g1).^3;
            dgdt = -obj.v0*(exp(-obj.q*obj.Eag./(obj.k*T)).*exp(gamma* ...
                   obj.a0*obj.q.*v./(obj.L*obj.k*T))-exp(-obj.q* ...
                   obj.Ear./(obj.k*T)).*exp(-gamma*obj.a0*obj.q.*v./ ...
                   (obj.L*obj.k*T)));
            dgdt(gregion == "gmin") = ...
                dgdt(gregion == "gmin").* ...
                (1+sign(dgdt(gregion == "gmin")))/2;
            dgdt(gregion == "gmax") = ...
                dgdt(gregion == "gmax").* ...
                (1-sign(dgdt(gregion == "gmax")))/2;
        end

        function dTdt = temperaturetimederivative(obj,v,g,T)
            % TEMPERATURETIMEDERIVATIVE Compute the memristor temperature
            % time derivative as a function of voltage, gap and
            % temperature.
            % Args:
            % - v: voltage [V]
            % - i: gap [m]
            % - T: temperature [K]
            % Return:
            % - dTdt: temperature time derivative [K/s]
            i = obj.current(v,g);
            dTdt = abs(v.*i)/obj.Cth-(T-obj.T0)/obj.tauth;
        end

        function maxdgdt = maxgaptimederivative(obj,v,g,Tmax)
            % MAXGAPTIMEDERIVATIVE Compute the memristor maximum gap time 
            % derivative as a function of voltage and gap.
            % Args:
            % - v: voltage [V]
            % - g: gap [m]
            % - Tmax: maximum temperature [K]
            % Return:
            % - maxdgdt: gap time derivative [m/s]
            gamma = obj.gamma0-obj.beta*(g/obj.g1).^3;
            v1 = -obj.L*obj.Ear./(gamma*obj.a0);
            v2 = obj.L*obj.Eag./(gamma*obj.a0);
            T = zeros(size(v));
            T(v <= v1 | v >= v2) = obj.T0;
            T(v > v1 & v < v2) = Tmax;
            maxdgdt = -obj.v0*(exp(-obj.q*obj.Eag./(obj.k*T)).* ...
                      exp(gamma*obj.a0*obj.q.*v./(obj.L*obj.k*T))- ...
                      exp(-obj.q*obj.Ear./(obj.k*T)).* ...
                      exp(-gamma*obj.a0*obj.q.*v./(obj.L*obj.k*T)));
        end

        function staticdgdt = staticgaptimederivative(obj,v,g)
            % STATICGAPTIMEDERIVATIVE Compute the memristor gap time
            % derivative assuming the temperature to always be at the
            % equilibrium.
            % Args:
            % - v: voltage [V]
            % - g: gap [m]
            % Return:
            % - staticdgdt: gap time derivative [m/s]
            i = obj.current(v,g);
            T = obj.T0+(obj.tauth/obj.Cth)*abs(v.*i);
            gamma = obj.gamma0-obj.beta*(g/obj.g1).^3;
            staticdgdt = -obj.v0*(exp(-obj.q*obj.Eag./(obj.k*T)).* ...
                         exp(gamma*obj.a0*obj.q.*v./(obj.L*obj.k*T))- ...
                         exp(-obj.q*obj.Ear./(obj.k*T)).* ...
                         exp(-gamma*obj.a0*obj.q.*v./(obj.L*obj.k*T)));
        end
    end
end