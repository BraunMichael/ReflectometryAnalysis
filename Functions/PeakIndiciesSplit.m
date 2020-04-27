function [top, middle, bottom] = PeakIndiciesSplit(zeros_xindices, zeros_1stderiv_indices, zeros_2ndderiv_indices, peaks_yvalues)
top = [];
middle = [];
bottom = [];

if zeros_1stderiv_indices(1) > zeros_2ndderiv_indices(1)
    %Then 1st one has to be a middle
    if peaks_yvalues(1) < peaks_yvalues(2)
        %Then its on the downhill
        %for loop here to pick the top, middle, bottom values
        for i = 0:floor(max(size(zeros_xindices))/4)-1
            top = [top, zeros_xindices(4*i+4)];
            middle = [middle, zeros_xindices(4*i+1), zeros_xindices(4*i+3)];
            bottom = [bottom, zeros_xindices(4*i+2)];
        end
        remain = mod(max(size(zeros_xindices)), 4);
        switch remain
            case 1
                middle = [middle, zeros_xindices(4*i+4+1)];
            case 2
                middle = [middle, zeros_xindices(4*i+4+1)];
                bottom = [bottom, zeros_xindices(4*i+4+2)];
            case 3
                middle = [middle, zeros_xindices(4*i+4+1)];
                bottom = [bottom, zeros_xindices(4*i+4+2)];
                middle = [middle, zeros_xindices(4*i+4+3)];
        end
    else
        %Then its on the uphill
        for i = 0:floor(max(size(zeros_xindices))/4)-1
            top = [top, zeros_xindices(4*i+2)];
            middle = [middle, zeros_xindices(4*i+1), zeros_xindices(4*i+3)];
            bottom = [bottom, zeros_xindices(4*i+4)];
        end
        remain = mod(max(size(zeros_xindices)), 4);
        switch remain
            case 1
                middle = [middle, zeros_xindices(4*i+4+1)];
            case 2
                middle = [middle, zeros_xindices(4*i+4+1)];
                top = [top, zeros_xindices(4*i+4+2)];
            case 3
                middle = [middle, zeros_xindices(4*i+4+1)];
                top = [top, zeros_xindices(4*i+4+2)];
                middle = [middle, zeros_xindices(4*i+4+3)];
        end
    end
else
    %1st one must be a top or bottom
    if peaks_yvalues(1) > peaks_yvalues(2)
        %Then it starts on a top
        %for loop here to pick the top, middle, bottom values
        for i = 0:floor(max(size(zeros_xindices))/4)-1
            top = [top, zeros_xindices(4*i+1)];
            middle = [middle, zeros_xindices(4*i+2), zeros_xindices(4*i+4)];
            bottom = [bottom, zeros_xindices(4*i+3)];
        end
        remain = mod(max(size(zeros_xindices)), 4);
        switch remain
            case 1
                top = [top, zeros_xindices(4*i+4+1)];
            case 2
                top = [top, zeros_xindices(4*i+4+1)];
                middle = [middle, zeros_xindices(4*i+4+2)];
            case 3
                top = [top, zeros_xindices(4*i+4+1)];
                middle = [middle, zeros_xindices(4*i+4+2)];
                bottom = [bottom, zeros_xindices(4*i+4+3)];
        end
    else
        %Then it starts on a bottom
        %for loop here to pick the top, middle, bottom values
        for i = 0:floor(max(size(zeros_xindices))/4)-1
            bottom = [bottom, zeros_xindices(4*i+1)];
            middle = [middle, zeros_xindices(4*i+2), zeros_xindices(4*i+4)];
            top = [top, zeros_xindices(4*i+3)];
        end
        remain = mod(max(size(zeros_xindices)), 4);
        switch remain
            case 1
                bottom = [bottom, zeros_xindices(4*i+4+1)];
            case 2
                bottom = [bottom, zeros_xindices(4*i+4+1)];
                middle = [middle, zeros_xindices(4*i+4+2)];
            case 3
                bottom = [bottom, zeros_xindices(4*i+4+1)];
                middle = [middle, zeros_xindices(4*i+4+2)];
                top = [top, zeros_xindices(4*i+4+3)];
        end
    end
end