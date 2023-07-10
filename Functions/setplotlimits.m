% Company: University of Siena
% Engineer: Riccardo Moretti
% Project: CurrentLimitedMemristor
%
% Description: Set a plot axes limits
%
% Args:
% - x: x axis data
% - y: y axis data
% - yscale: y axis scale
%
% Returns:
% - xlim: x axis limits
% - ylim: y axis limits

function [xlim,ylim] = setplotlimits(x,y,yscale)
    xlim = zeros(1,2);
	if ~isempty(x)
        xlim(1) = min(x,[],'all');
        xlim(2) = max(x,[],'all'); 
		if xlim(1) == xlim(2)
			xlim(1) = 0.9*xlim(1);
			xlim(2) = 1.1*xlim(2);
            xlim = sort(xlim);
		end
	end
    ylim = zeros(1,2);
    if ~isempty(y)
        if max(y,[],'all')-min(y,[],'all') == 0
            if yscale == "linear"
                ylim(1) = 0.9*mean(y,'all');
                ylim(2) = 1.1*mean(y,'all');
                ylim = sort(ylim);
            elseif yscale == "logarithmic"
                ylim(1) = min(y,[],'all')^0.9;
                ylim(2) = max(y,[],'all')^1.1;
                ylim = sort(ylim);
            end
        else
            if yscale == "linear"
                dy = max(y,[],'all')-min(y,[],'all');
                ylim(1) = min(y,[],'all')-0.1*dy;
                ylim(2) = max(y,[],'all')+0.1*dy;
            elseif yscale == "logarithmic"
                dy = max(y,[],'all')/min(y,[],'all');
                ylim(1) = min(y,[],'all')/(dy^0.1);
                ylim(2) = max(y,[],'all')*(dy^0.1);
            end
        end
    end
end