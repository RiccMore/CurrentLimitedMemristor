% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: NMOS limited Stanford memristor excitation

classdef nmoslimitedstanfordmemristorexcitation < handle
        properties (SetAccess = immutable)
        Tmax = 500; % Breakdown temperature [K]
    end

    properties
        V1; % Pulse voltage source
        V2; % Controlled pulse voltage source
        V3; % DC voltage source
        M; % Stanford memristor
        N; % NMOS
    end

    properties (Hidden = true)
        gn; % Memristor gap order of magnitude
        Tn; % Memristor temperature order of magnitude
        gminthreshold; % Memristor gmin gap region threshold [m]
        gmaxthreshold; % Memristor gmax gap region threshold [m]
        gregion; % Memristor gap region
    end

    methods
        function obj = nmoslimitedstanfordmemristorexcitation(varargin)
            % NMOSLIMITEDSTANFORDMEMRISTOREXCITATION Class constructor.
            % Properties:
            % - V1: pulse voltage source
            % - V2: controlled pulse voltage source
            % - V3: DC voltage source
            % - M: Stanford memristor
            % - N: NMOS
            p = inputParser;
            p.addParameter('V1',pulsevoltage);
            p.addParameter('V2',controlledpulsevoltage);
            p.addParameter('V3',dcvoltage);
            p.addParameter('M',stanfordmemristor);
            p.addParameter('N',nmos);
            p.parse(varargin{:});
            obj.V1 = p.Results.V1;
            obj.V2 = p.Results.V2;
            obj.V3 = p.Results.V3;
            obj.M = p.Results.M;
            obj.N = p.Results.N;
            obj.gn = 10^(3*round(floor(log10(abs(obj.M.gmin)))/3));
            obj.Tn = 10^(3*round(floor(log10(abs(obj.M.T0)))/3));
            obj.gminthreshold = obj.M.gmin+13*eps(obj.M.gmin);
            obj.gmaxthreshold = obj.M.gmax-2*eps(obj.M.gmax);
            if obj.M.gic <= obj.gminthreshold
                obj.gregion = "gmin";
            elseif obj.M.gic >= obj.gmaxthreshold
                obj.gregion = "gmax";
            else
                obj.gregion = "";
            end
        end

        function sizer(obj,gtarget,dgdtlim)
            % SIZER Set the NMOS gate voltage value to reach a 
            % certain gap value with a set operation.
            % Args:
            % - gtarget: target memristor gap [m]
            % - dgdtlim: gap speed resolution [m/s]
            v = obj.V1.Voff+obj.V1.V0;
            VGS = 0.61:min(1e-3,1.8-obj.V3.Vdc):(1.8-obj.V3.Vdc);
            vDS = 0:min(1e-3,1.8-obj.V3.Vdc):(1.8-obj.V3.Vdc);
            vSB = obj.V3.Vdc*ones(size(vDS));
            vM = v-vDS;
            for a = 1:length(VGS)
                vGS = VGS(a)*ones(size(vDS));
                iM = obj.M.current(vM,gtarget);
                iD = obj.N.current(vGS,vSB,vDS);
                vds = min(interp1(iM-iD,vDS,0,'linear','extrap'));
                vm = v-vds;
                dgdt = obj.M.staticgaptimederivative(vm,gtarget);
                if dgdt <= dgdtlim
                    break
                end
            end
            obj.V2.Vlow = VGS(a);
        end

        function sol = solver(obj,t,g)
            % SOLVER Compute the circuit nodes potentials and branches
            % currents at a given time instant.
			% Arg:
			% - t: time [s]
			% Return:
			% - sol: solution structure
            %        sol.v: nodes potentials [V]
            %        sol.i: branches currents [A]
            sol.v = zeros(4,length(t));
            sol.v(1,:) = obj.V1.voltage(t);
            sol.v(3,:) = obj.V2.voltage(sol.v(1,:));
            sol.v(4,:) = -obj.V3.Vdc*ones(size(t));
            if length(t) > 1
                strclr = '\b';
                fprintf('Signals calculation: \n');
            end
            for a = 1:length(t)
                if length(t) > 1
                    progress = floor(100*a/length(t));
                    fprintf(strclr);
                    strout = sprintf('%d',progress);
                    strclr = repmat('\b',1,length(strout)+2);
                    fprintf(strcat(strout,'%%\n'));
                end
                if sol.v(1,a) > 0
                    vDS = 0:min(1e-3,sol.v(1,a)):sol.v(1,a);
                    vGS = sol.v(3,a)*ones(size(vDS));
                    vSB = -sol.v(4,a)*ones(size(vDS));
                    vM = sol.v(1,a)-vDS;
                    iM = obj.M.current(vM,g(a));
                    iD = obj.N.current(vGS,vSB,vDS);
                    vDS = interp1(iM-iD,vDS,0,'linear','extrap');
                    sol.v(2,a) = min(vDS);
                elseif sol.v(1,a) < 0
                    vDS = sol.v(1,a):min(1e-3,-sol.v(1,a)):0;
                    vGS = sol.v(3,a)*ones(size(vDS));
                    vSB = -sol.v(4,a)*ones(size(vDS));
                    vM = sol.v(1,a)-vDS;
                    iM = obj.M.current(vM,g(a));
                    iD = obj.N.current(vGS,vSB,vDS);
                    vDS = interp1(iM-iD,vDS,0,'linear','extrap');
                    sol.v(2,a) = max(vDS);
                end
            end
            sol.i = obj.N.current(sol.v(3,:),-sol.v(4,:),sol.v(2,:));
        end

        function dxdt = odemodel(obj,t,x)
            % ODEMODEL Solve the circuit state variables time derivatives
            % at a given time instant.
            % Args:
            % - t: time instant [s]
            % - x: state variables
            %      g: gap (normalized magnitude)
            %      T: temperature (normalized magnitude)
            % Return:
            % - dxdt: state variables time derivatives
            %         dgdt: gap derivative (normalized magnitude)
            %         dTdt: temperature derivative (normalized magnitude)
            g = x(1,:)*obj.gn;
            T = x(2,:)*obj.Tn;
            sol = obj.solver(t,g);
            v = sol.v(1,:)-sol.v(2,:);
            dgdt = obj.M.gaptimederivative(v,g,T,obj.gregion);
            dTdt = obj.M.temperaturetimederivative(v,g,T);
            dxdt = [dgdt;dTdt]./[obj.gn;obj.Tn];
        end

        function [value,isterminal,direction] = odeevents(obj,~,x)
            % ODEEVENTS Check the events on the circuit state variables
            % dynamics.
            % Args:
            % - x: state variables
            %      g: gap (normalized magnitude)
            %      T: temperature (normalized magnitude)
            % Returns:
            % - value: events threshold values
            % - isterminal: events simulation terminate conditions
            % - direction: events threshold crossing directions
            g = x(1)*obj.gn;
            T = x(2)*obj.Tn;
            gvalue = g*ones(4,1)-[obj.M.gmin;
                                  obj.gminthreshold;
                                  obj.gmaxthreshold;
                                  obj.M.gmax];
            gisterminal = ones(4,1);
            gdirection = [-1;1;-1;1];
            Tvalue = T-obj.Tmax;
            Tisterminal = 1;
            Tdirection = 1;
            value = [gvalue;Tvalue];
            isterminal = [gisterminal;Tisterminal];
            direction = [gdirection;Tdirection];
        end

        function [sim,fig] = transient(obj,tspan,varargin)
            % TRANSIENT Simulate the circuit on a given time interval.
            % Arg:
            % - tspan: time interval [s]
            % Property:
            % - PlotOutput: plot results enable
            % Returns:
            % - sim: simulation signals
            %        sim.t: time [s]
            %        sim.v: nodes potentials [V]
            %        sim.i: branches currents [A]
            %        sim.g: memristor gap [m]
            %        sim.T: memristor temperature [K]
            %        sim.dgdt: memristor gap time derivative [m/s]
            %        sim.dTdt: memristor temperature time derivative [K/s]
            % - fig: simulation plot
            p = inputParser;
            p.addParameter('PlotOutput',true);
            p.parse(varargin{:});
            PlotOutput = p.Results.PlotOutput;
            if obj.M.gic <= obj.gminthreshold
                obj.gregion = "gmin";
            elseif obj.M.gic >= obj.gmaxthreshold
                obj.gregion = "gmax";
            else
                obj.gregion = "";
            end
            ic = [obj.M.gic;obj.M.Tic]./[obj.gn;obj.Tn];
            options = odeset('Events',@obj.odeevents,'RelTol',1e-6, ...
	                         'AbsTol',1e-6,'MaxStep',tspan(2)/1000, ...
				             'OutputFcn',@odeprogress);
            burntmemristor = false;
            odesol = ode15s(@obj.odemodel,tspan,ic,options);
            greg = repmat(obj.gregion,[1 length(odesol.x)]);
            while odesol.x(end) < tspan(2)
	            if length(odesol.ie) > 1
                    eventstr = sprintf('%e ',odesol.xe);
                    warning('Multiple events at %s s',eventstr);
	            end
	            for a = 1:length(odesol.ie)
		            switch odesol.ie(a)
			            case 1
				            obj.gregion = "gmin";
				            fprintf('Event: memristor gap reached gmin\n');
			            case 2
				            obj.gregion = "";
				            fprintf(['Event: memristor gap greater ' ...
                                     'than gmin\n']);
			            case 3
				            obj.gregion = "";
				            fprintf(['Event: memristor gap smaller ' ...
                                     'than gmax\n']);
			            case 4
				            obj.gregion = "gmax";
				            fprintf('Event: memristor gap reached gmax\n');
			            case 5
				            burntmemristor = true;
				            fprintf(['Event: memristor temperature ' ...
                                     'reached Tmax\n']);
			            otherwise
				            error('Unexpected event');
		            end
	            end
                odesol.ie = [];
                odesol.xe = [];
                odesol.ye = [];
                if burntmemristor
                    fprintf('Memristor burnt\n');
                    break
                end
                odesol = odextend(odesol,@obj.odemodel,tspan(2), ...
                                  odesol.y(:,end),options);
                greg = [greg repmat(obj.gregion, ...
                                    [1 length(odesol.x)-length(greg)])];
            end
            sim.t = odesol.x;
            sim.g = odesol.y(1,:)*obj.gn;
            sim.T = odesol.y(2,:)*obj.Tn;
            sol = obj.solver(sim.t,sim.g);
            sim.v = sol.v;
            sim.i = sol.i;
            obj.gregion = greg;
            v = sim.v(1,:)-sim.v(2,:);
            sim.dgdt = obj.M.gaptimederivative(v,sim.g,sim.T, ...
                                               obj.gregion);
            sim.dTdt = obj.M.temperaturetimederivative(v,sim.g,sim.T);
            if PlotOutput
                fprintf('Plotting simulation...\n');
                plotfontsize = 14;
                fig = figure('Name','TimeSimulation', ...
                             'Units','normalized', ...
                             'OuterPosition',[0 0 1 1]);
                sgtitle('Simulation vs Time', ...
                        'FontSize',1.5*plotfontsize,'FontName','Arial');
                subplot(3,2,1);
                [x,~,ox] = sinotation(sim.t);
                [y1,ny,oy] = sinotation(sim.v(1,:));
                y2 = sim.v(2,:)/(10^(3*ny));
                y3 = y1-y2;
                [xlim,ylim] = setplotlimits(x,y1,"linear");
                plot(x,y1,'LineWidth',2,'DisplayName','Excitation');
                hold on;
                plot(x,y2,'LineWidth',2,'DisplayName','NMOS');
                plot(x,y3,'LineWidth',2,'DisplayName','Memristor');
                title('Voltage');
                xlabel(sprintf('Time [%ss]',ox));
                ylabel(sprintf('Voltage [%sV]',oy));  
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                legend('Location','best','FontSize',0.75*plotfontsize, ...
                       'FontName','Arial');
                box on; grid on;
                subplot(3,2,2);
                [x,~,ox] = sinotation(sim.t);
                [y,~,oy] = sinotation(sim.i);
                [xlim,ylim] = setplotlimits(x,y,"linear");
                plot(x,y,'LineWidth',2);
                title('Current');
                xlabel(sprintf('Time [%ss]',ox));
                ylabel(sprintf('Current [%sA]',oy));
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                box on; grid on;
                subplot(3,2,3);
                [x,~,ox] = sinotation(sim.t);
                [~,ny,oy] = sinotation(obj.M.gmax);
                y = sim.g/(10^(3*ny));
                [xlim,ylim] = setplotlimits(x, ...
                              [obj.M.gmin obj.M.gmax]/(10^(3*ny)), ...
                              "linear");
                plot(x,y,'LineWidth',2); hold on;
                yline(obj.M.gmin/(10^(3*ny)),'--k','LineWidth',2);
                yline(obj.M.gmax/(10^(3*ny)),'--k','LineWidth',2);
                title('Gap');
                xlabel(sprintf('Time [%ss]',ox));
                ylabel(sprintf('Gap [%sm]',oy));
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                box on; grid on;
                subplot(3,2,4);
                [x,~,ox] = sinotation(sim.t);
                [y,~,oy] = sinotation(sim.T);
                [xlim,ylim] = setplotlimits(x,y,"linear");
                plot(x,y,'LineWidth',2);
                title('Temperature');
                xlabel(sprintf('Time [%ss]',ox));
                ylabel(sprintf('Temperature [%sK]',oy));
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                box on; grid on;
                subplot(3,2,5);
                [x,~,ox] = sinotation(sim.t);
                [y,~,oy] = sinotation(sim.dgdt);
                [xlim,ylim] = setplotlimits(x,y,"linear");
                plot(x,y,'LineWidth',2);
                title('Gap Time Derivative');
                xlabel(sprintf('Time [%ss]',ox));
                ylabel(sprintf('Gap/Time [%sm/ns]',oy));
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                box on; grid on;
                subplot(3,2,6);
                [x,~,ox] = sinotation(sim.t);
                [y,~,oy] = sinotation(sim.dTdt/1e9);
                [xlim,ylim] = setplotlimits(x,y,"linear");
                plot(x,y,'LineWidth',2);
                title('Temperature Time Derivative');
                xlabel(sprintf('Time [%ss]',ox));
                ylabel(sprintf('Temperature/Time [%sK/ns]',oy));
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                box on; grid on;
                fprintf('\b\b\b\b');
                fprintf(' done.\n');
                drawnow;
            else
                fig = [];
            end
        end

        function [postproc,fig] = tranprocessing(obj,sim,varargin)
            % TRANPROCESSING Project the memristor gap time derivative
            % simulation trajectory on its maximum case contour surface,
            % evaluated with respect to the memristor gap and voltage
            % values.
            % Arg:
            % - sim: simulation signals
            %        sim.t: time [s]
            %        sim.v: nodes potentials [V]
            %        sim.i: branches currents [A]
            %        sim.g: memristor gap [m]
            %        sim.T: memristor temperature [K]
            %        sim.dgdt: memristor gap time derivative [m/s]
            %        sim.dTdt: memristor temperature time derivative [K/s]
            % Property:
            % - PlotOutput: plot results enable
            % Returns:
            % - postproc: contour surface data
            %             postproc.G: memristor gap mesh grid [m]
            %             postproc.V: memristor voltage mesh grid [V]
            %             postproc.DGDT: memristor maximum gap time
            %                            derivative surface [m/s]
            % - fig: contour surface projection plot
            p = inputParser;
            p.addParameter('PlotOutput',true);
            p.parse(varargin{:});
            PlotOutput = p.Results.PlotOutput;
            gstep = 10^floor(log10(obj.M.gmax-obj.M.gmin)-3);
            GG = obj.M.gmin:gstep:obj.M.gmax;
            vstep = 10^floor(log10(2*obj.V1.V0)-3);
            VV = obj.V1.Voff+(-obj.V1.V0:vstep:obj.V1.V0);
            [GG,VV] = meshgrid(GG,VV);
            DGDT = obj.M.maxgaptimederivative(VV,GG,obj.Tmax);
            strclr = '\b';
            fprintf('Maximum gap derivative computation progress: \n');
            vlim = obj.V1.Voff+[-obj.V1.V0 obj.V1.V0];
            for a = 1:size(GG,2)
                progress = floor(100*a/size(GG,2));
                fprintf(strclr);
                strout = sprintf('%d',progress);
                strclr = repmat('\b',1,length(strout)+2);
                fprintf(strcat(strout,'%%\n'));
                vDS = vlim(1):min(1e-3,-vlim(1)):0;
                vGS = obj.V2.voltage(vlim(1))*ones(size(vDS));
                vSB = obj.V3.Vdc*ones(size(vDS));
                vM = vlim(1)-vDS;
                iM = obj.M.current(vM,GG(1,a));
                iD = obj.N.current(vGS,vSB,vDS);
                vDS = interp1(iM-iD,vDS,0,'linear','extrap');
                vMlow = vlim(1)-max(vDS);
                vDS = 0:min(1e-3,vlim(2)):vlim(2);
                vGS = obj.V2.voltage(vlim(2))*ones(size(vDS));
                vSB = obj.V3.Vdc*ones(size(vDS));
                vM = vlim(2)-vDS;
                iM = obj.M.current(vM,GG(1,a));
                iD = obj.N.current(vGS,vSB,vDS);
                vDS = interp1(iM-iD,vDS,0,'linear','extrap');
                vMhigh = vlim(2)-min(vDS);
                DGDT(VV(:,a) < vMlow | VV(:,a) > vMhigh,a) = nan;
            end
            postproc.G = GG;
            postproc.V = VV;
            postproc.DGDT = DGDT;
            if PlotOutput
                fprintf('Plotting processed simulation...\n');
                plotfontsize = 14;
                fig = figure('Name','ContourPlot', ...
                             'Units','normalized', ...
                             'OuterPosition',[0 0 1 1]);
                [x,nx,ox] = sinotation(GG);
                [y,ny,oy] = sinotation(VV);
                z = 20*log10(abs(DGDT));
                z(z == -Inf) = min(z(~isinf(z)),[],'all');
                [xlim,~] = setplotlimits(x,[],"linear");
                [ylim,~] = setplotlimits(y,[],"linear");
                [zlim,~] = setplotlimits(z,[],"linear");
                levels = (floor(zlim(1)/20):ceil(zlim(2)/20))*20;
                contourf(x,y,z,levels); hold on;
                x = sim.g/(10^(3*nx));
                y = (sim.v(1,:)-sim.v(2,:))/(10^(3*ny));
                z = zeros(size(sim.dgdt));
                plot3(x,y,z,'r','LineWidth',4);
                ticks = 100*(floor(min(levels)/100):ceil(max(levels)/100));
                cb = colorbar('Ticks',ticks, ...
                              'Limits',[min(ticks) max(ticks)]);
                ylabel(cb,'|Worst-Case Gap Derivative| [dB]', ...
                       'FontSize',plotfontsize,'FontName','Arial', ...
                       'Rotation',270);
                title('Gap Derivative Simulation vs Memristor Domain');
                xlabel(sprintf('Gap [%sm]',ox));
                ylabel(sprintf('Voltage [%sV]',oy));
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                fprintf('\b\b\b\b');
                fprintf(' done.\n');
                drawnow;
            else
                fig = [];
            end
        end

        function sim = steptransient(obj,tspan,varargin)
            % STEPTRANSIENT Simulate the circuit on a given time
            % interval for different values of the step parameter.
            % Arg:
            % - tspan: time interval [s]
            % Properties:
            % - V0: V1 amplitude [V]
            % - f0: V1 frequency [Hz]
            % - Vlow: V2 low level voltage [V]
            % - n: NMOS aspect ratio factor
            % Return:
            % - sim: simulations signals
            %        sim.x: step parameter
            %        sim.transim: transient simulation
            p = inputParser;
            p.addParameter('V0',[]);
            p.addParameter('f0',[]);
            p.addParameter('Vlow',[]);
            p.addParameter('n',[]);
            p.parse(varargin{:});
            if ~isempty(p.Results.V0)
                sim = struct('V0',cell(size(p.Results.V0)), ...
                             'tran',[]);
                for a = 1:length(p.Results.V0)
                    clc;
                    progress = floor(100*a/length(p.Results.V0));
                    fprintf('Step simulation progress: %d%%\n',progress);
                    obj.V1.V0 = p.Results.V0(a);
                    sim(a).V0 = p.Results.V0(a);
                    sim(a).tran = obj.transient(tspan,'PlotOutput',false);
                end
            elseif ~isempty(p.Results.f0)
                sim = struct('T0',cell(size(p.Results.f0)), ...
                             'tran',[]);
                for a = 1:length(p.Results.f0)
                    clc;
                    progress = floor(100*a/length(p.Results.f0));
                    fprintf('Step simulation progress: %d%%\n',progress);
                    obj.V1.f0 = p.Results.f0(a);
                    sim(a).tran = obj.transient([tspan(1) tspan(a+1)], ...
                                                'PlotOutput',false);
                    sim(a).T0 = 1/(4*p.Results.f0(a));
                end
            elseif ~isempty(p.Results.Vlow)
                sim = struct('Vlow',cell(size(p.Results.Vlow)), ...
                             'tran',[]);
                for a = 1:length(p.Results.Vlow)
                    clc;
                    progress = floor(100*a/length(p.Results.Vlow));
                    fprintf('Step simulation progress: %d%%\n',progress);
                    obj.V2.Vlow = p.Results.Vlow(a);
                    sim(a).Vlow = p.Results.Vlow(a);
                    sim(a).tran = obj.transient(tspan,'PlotOutput',false);
                end
            elseif ~isempty(p.Results.n)
                sim = struct('ratio',cell(size(p.Results.n)), ...
                             'tran',[]);
                for a = 1:length(p.Results.n)
                    clc;
                    progress = floor(100*a/length(p.Results.n));
                    fprintf('Step simulation progress: %d%%\n',progress);
                    obj.N.n = p.Results.n(a);
                    sim(a).ratio = p.Results.n(a)*obj.N.W/obj.N.L;
                    sim(a).tran = obj.transient(tspan,'PlotOutput',false);
                end
            end
        end

        function [postproc,fig] = steptranprocessing(obj,sim,varargin)
            % STEPTRANPROCESSING Evaluate the set gap value sensitivity to
            % the step parameter.
            % Arg:
            % - sim: simulations signals
            %        sim.x: step parameter
            %        sim.transim: transient simulation
            % Property:
            % - PlotOutput: plot results enable
            % Return:
            % - postproc: mapping data
            %             postproc.x: step parameter
            %             postproc.g: final gap [m]
            %             postproc.dgdx: final gap sensitivity to step
            %                            parameter
            % - fig: sensitivity plot
            p = inputParser;
            p.addParameter('PlotOutput',true);
            p.parse(varargin{:});
            PlotOutput = p.Results.PlotOutput;
            simfields = fields(sim);
            gfinal = zeros(size(sim));
            steppar = zeros(size(sim));
            dgdx = zeros(size(sim));
            strclr = '\b';
            fprintf('Gap mapping computation progress: \n');
            for a = 1:length(sim)
                progress = floor(100*a/length(sim));
                fprintf(strclr);
                strout = sprintf('%d',progress);
                strclr = repmat('\b',1,length(strout)+2);
                fprintf(strcat(strout,'%%\n'));
                gfinal(a) = sim(a).tran.g(end);
                if isequal(simfields{1},'V0')
                    steppar(a) = sim(a).V0;
                    stepparameter = 'Pulse Amplitude';
                    steplabel = 'Voltage [';
                    stepunit = 'V]';
                elseif isequal(simfields{1},'T0')
                    steppar(a) = sim(a).T0;
                    stepparameter = 'Pulse Width';
                    steplabel = 'Time [';
                    stepunit = 's]';
                elseif isequal(simfields{1},'Vlow')
                    steppar(a) = sim(a).Vlow;
                    stepparameter = 'Gate Voltage';
                    steplabel = 'Voltage [';
                    stepunit = 'V]';
                elseif isequal(simfields{1},'ratio')
                    steppar(a) = sim(a).ratio;
                    stepparameter = 'Aspect Ratio';
                    steplabel = 'Aspect Ratio';
                    stepunit = '';
                end
            end
            [steppar,isteppar] = sort(steppar);
            gfinal = gfinal(isteppar);
            postproc.x = linspace(min(steppar),max(steppar), ...
                                  length(steppar));
            postproc.g = interp1(steppar,gfinal,postproc.x,'spline');
            dgdx = diff(postproc.g)./diff(postproc.x);
            postproc.dgdx = interp1(postproc.x(1:end-1),dgdx, ...
                                    postproc.x,'linear','extrap');
            if PlotOutput
                fprintf('Plotting simulation...\n');
                plotfontsize = 14;
                fig = figure('Name','GapMapping', ...
                            'Units','normalized', ...
                            'OuterPosition',[0 0 1 1]);
                sgtitle(sprintf('Gap Mapping vs %s',stepparameter), ...
                        'FontSize',1.5*plotfontsize,'FontName','Arial');
                subplot(1,2,1);
                if isequal(simfields{1},'ratio')
                    x = postproc.x; ox = '';
                else
                    [x,~,ox] = sinotation(postproc.x);
                end
                [~,ny,oy] = sinotation(obj.M.gmax);
                y = postproc.g/(10^(3*ny));
                [xlim,ylim] = setplotlimits(x, ...
                              [obj.M.gmin obj.M.gmax]/(10^(3*ny)), ...
                              "linear");
                plot(x,y,'LineWidth',2); hold on;
                yline(obj.M.gmin/(10^(3*ny)),'--r','LineWidth',2);
                yline(obj.M.gmax/(10^(3*ny)),'--r','LineWidth',2);
                title('Gap');
                xlabel(strcat(steplabel,sprintf('%s',ox),stepunit));
                ylabel(sprintf('Gap [%sm]',oy));
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                box on; grid on;
                subplot(1,2,2);
                if isequal(simfields{1},'ratio')
                    x = postproc.x; ox = '';
                else
                    [x,~,ox] = sinotation(postproc.x);
                end
                [y,~,oy] = sinotation(postproc.dgdx);
                [xlim,ylim] = setplotlimits(x,y,"linear");
                plot(x,y,'LineWidth',2);
                title('Gap Derivative');
                xlabel(strcat(steplabel,sprintf('%s',ox),stepunit));
                if isequal(simfields{1},'ratio')
                    ylabel(sprintf('Gap Derivative [%sm]',oy));
                else
                    ylabel(sprintf('Gap Derivative [%sm/%s',oy,stepunit));
                end
                set(gca,'XLim',xlim,'YLim',ylim, ...
                    'FontSize',plotfontsize,'FontName','Arial', ...
                    'GridAlpha',0.8,'MinorGridAlpha',0.8);
                box on; grid on;
                drawnow;
            else
                fig = [];
            end
        end
    end
end