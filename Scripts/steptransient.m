% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Parameter step transient simulator

%% INITIALIZE WORKSPACE
close all;
clear; clc;
scriptpath = mfilename('fullpath');
[~,eop] = find(scriptpath == '\',2,'last');
cd(scriptpath(1:eop(1)));
addpath(genpath(pwd));
clear;

%% DEFINE PARAMETERS
savedata = true;
fig = gobjects(1,1);

%% DEFINE CIRCUIT
% circuit = "stanfordmemristorexcitation";
% circuit = "resistorlimitedstanfordmemristorexcitation";
circuit = "nmoslimitedstanfordmemristorexcitation";
if circuit == "stanfordmemristorexcitation"
    V = pulsevoltage('Voff',0, ...
                     'V0',2, ...
                     'f0',500, ...
                     'phi0',0);
    M = stanfordmemristor();
    M.gic = M.gmax;
    M.Tic = M.T0;
    C = stanfordmemristorexcitation('V',V, ...
                                    'M',M);
    stepparameter = 'V0';
    step = linspace(1.7,1.85,1001);
    tspan = [0 1/(2*V.f0)];
    clear V M;
elseif circuit == "resistorlimitedstanfordmemristorexcitation"
    V = pulsevoltage('Voff',0, ...
                     'V0',2, ...
                     'f0',500, ...
                     'phi0',0);
    M = stanfordmemristor();
    M.gic = M.gmax;
    M.Tic = M.T0;
    R = resistor();
    C = resistorlimitedstanfordmemristorexcitation('V',V, ...
                                                   'M',M, ...
                                                   'R',R);
    stepparameter = 'R';
    step = linspace(6e3,106e3,101);
    % stepparameter = 'V0';
    % step = linspace(1.8,2.2,101);
    % stepparameter = 'f0';
    % step = linspace(100,500,101);
    if isequal(stepparameter,'f0')
        tspan = [0 1./(2*step)];
    else
        tspan = [0 1/(2*V.f0)];
    end
    clear V M R;
elseif circuit == "nmoslimitedstanfordmemristorexcitation"
    V1 = pulsevoltage('Voff',0, ...
                      'V0',2, ...
                      'f0',500, ...
                      'phi0',0);
    V2 = controlledpulsevoltage('Vhigh',1.6, ...
                                'Vlow',1, ...
                                'Vthresh',0);
    V3 = dcvoltage('Vdc',0.2);
    M = stanfordmemristor();
    M.gic = M.gmax;
    M.Tic = M.T0;
    N = nmos('n',1);
    C = nmoslimitedstanfordmemristorexcitation('V1',V1, ...
                                               'V2',V2, ...
                                               'V3',V3, ...
                                               'M',M, ...
                                               'N',N);
    % stepparameter = 'Vlow';
    % step = linspace(0.6,1.6,101);
    % stepparameter = 'V0';
    % step = linspace(1.8,2.2,101);
    % stepparameter = 'f0';
    % step = 1./linspace(0.001,0.01,101);
    stepparameter = 'n';
    step = linspace(0.85,1.25,101);
    if isequal(stepparameter,'f0')
        tspan = [0 1./(2*step)];
    else
        tspan = [0 1/(2*V1.f0)];
    end
    clear V1 V2 V3 M N;
end

%% SIMULATE CIRCUIT
sim = C.steptransient(tspan,stepparameter,step);
clear step tspan

%% PROCESS THE SIMULATION
[postproc,fig(1)] = C.steptranprocessing(sim);

%% SAVE DATA
if savedata
    fprintf('Saving...\n')
    dirpath = strcat('Data\',circuit,'\steptransient\', ...
                     char(datetime('now','Format','yyMMdd_HHmmss')),'\');
    mkdir(dirpath);
    for a = 1:length(fig)
	    savefig(fig(a),strcat(dirpath,fig(a).Name,'.fig'));
	    saveas(fig(a),strcat(dirpath,fig(a).Name),'svg');
        saveas(fig(a),strcat(dirpath,fig(a).Name),'png');
    end
	clear circuit savedata fig a
    save(strcat(dirpath,'data.mat'),'-regexp','^(?!dirpath.*$).');
    fprintf('\b\b\b\b');
    fprintf(' done.\n');
end
