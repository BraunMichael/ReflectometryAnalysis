function [Root_Indexes, Root_Xvalues] = NumericalRootsFunction(Xdatapoints, Ydatapoints)
%Finds zeroes of a purely numeric function
%Exports 2 lists, one with the cloesest index, and one with the extrapolated actual x value of
%the roots
previousx = Xdatapoints(1);
previousy = Ydatapoints(1);
Root_Indexes = [];
Root_Xvalues = [];
intersection = 0;
for i=2:max(size(Xdatapoints))
    currentx = Xdatapoints(i);
    currenty = Ydatapoints(i);
    %currenty*previousy;
    if currenty * previousy <= 0
        if abs(currenty) > abs(previousy)
            Root_Indexes = [Root_Indexes, i-1];
        else
            Root_Indexes = [Root_Indexes, i];
        end
        root = (((currentx - previousx)/(currenty - previousy)) * (intersection-previousy)) + previousx;
        Root_Xvalues = [Root_Xvalues, root];
    end
    previousx = currentx;
    previousy = currenty;
end



