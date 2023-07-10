% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Print ODE progress
%
% Args:
% - t: simulation time instant [s]
% - flag: ODE finite state machine current state
%
% Return:
% - status: ODE stop condition

function status = odeprogress(t,~,flag,varargin)
    persistent tend strclr;
    if isempty(flag)
        tnow = mean(t);
        progress = floor(100*tnow/tend);
        fprintf(strclr);
        strout = sprintf('%d',progress);
        strclr = repmat('\b',1,length(strout)+2);
        fprintf(strcat(strout,'%%\n'));
        status = 0;
    else
        switch flag
            case 'init'
                tend = max(t);
                strclr = '\b';
                fprintf('ODE progress: \n');
            case 'done'
                tend = [];
                strclr = [];
            otherwise
        end
    end
end