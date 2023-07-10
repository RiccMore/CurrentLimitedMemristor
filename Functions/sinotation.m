% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Get a number SI notation
%
% Arg:
% - x: real number
%
% Returns:
% - y: scaled number
% - n: magnitude scaling factor
% - o: order of magnitude

function [y,n,o] = sinotation(x)
    n = floor(log10(max(abs(x),[],'all'))/3);
	if n < -8
		n = -8;
	elseif n > 8
		n = 8;
	end
    y = x/(10^(3*n));
    switch n
        case -8
            o = 'y';
        case -7
            o = 'z';
        case -6
            o = 'a';
        case -5
            o = 'f';
        case -4
            o = 'p';
        case-3
            o = 'n';
        case -2
            o = sprintf('\x03bc');
        case -1
            o = 'm';
        case 0
            o = '';
        case 1
            o = 'k';
        case 2
            o = 'M';
        case 3
            o = 'G';
        case 4
            o = 'T';
        case 5
            o = 'P';
        case 6
            o = 'E';
        case 7
            o = 'Z';
        case 8
            o = 'Y';
        otherwise
            o = '?';
    end
end