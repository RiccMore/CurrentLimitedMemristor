% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Transient simulator

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
fig = gobjects(2,1);

%% DEFINE CIRCUIT
% circuit = "stanfordmemristorexcitation";
circuit = "resistorlimitedstanfordmemristorexcitation";
% circuit = "nmoslimitedstanfordmemristorexcitation";
if circuit == "stanfordmemristorexcitation"
    V = pulsevoltage('Voff',0, ...
                     'V0',1.8336316112501002, ...
                     'f0',1000, ...
                     'phi0',-pi);
    M = stanfordmemristor();
    M.gic = mean([M.gmin M.gmax]);
    M.Tic = M.T0;
    C = stanfordmemristorexcitation('V',V, ...
                                    'M',M);
    tspan = [0 1/(V.f0)];
    clear V M;
elseif circuit == "resistorlimitedstanfordmemristorexcitation"
    V = pulsevoltage('Voff',0, ...
                     'V0',2, ...
                     'f0',500, ...
                     'phi0',-pi);
    M = stanfordmemristor();
    M.gic = mean([M.gmin M.gmax]);
    M.Tic = M.T0;
    R = resistor();
    C = resistorlimitedstanfordmemristorexcitation('V',V, ...
                                                   'M',M, ...
                                                   'R',R);
    gtarget = 0.1e-9;
    dgdtlim = -0.1e-6;
    C.sizer(gtarget,dgdtlim);
    tspan = [0 2/V.f0];
    clear V M R dgdtlim;
elseif circuit == "nmoslimitedstanfordmemristorexcitation"
    V1 = pulsevoltage('Voff',0, ...
                      'V0',2, ...
                      'f0',500, ...
                      'phi0',-pi);
    V2 = controlledpulsevoltage('Vhigh',1.6, ...
                                'Vthresh',0);
    V3 = dcvoltage('Vdc',0.2);
    M = stanfordmemristor();
    M.gic = mean([M.gmin M.gmax]);
    M.Tic = M.T0;
    N = nmos();
    N = nmos('n',(5.3e-6)/N.W);
    C = nmoslimitedstanfordmemristorexcitation('V1',V1, ...
                                               'V2',V2, ...
                                               'V3',V3, ...
                                               'M',M, ...
                                               'N',N);
    gtarget = 0.1e-9;
    dgdtlim = -0.1e-6;
    C.sizer(gtarget,dgdtlim);
    % C.V2.Vlow = 1.6;
    tspan = [0 2/V1.f0];
    clear V1 V2 V3 M N dgdtlim;
end

%% SIMULATE CIRCUIT
[sim,fig(1)] = C.transient(tspan); 

%% PROCESS SIMULATION
[postproc,fig(2)] = C.tranprocessing(sim);

%% SAVE DATA
if savedata
    fprintf('Saving...\n')
    dirpath = strcat('Data\',circuit,'\transient\', ...
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