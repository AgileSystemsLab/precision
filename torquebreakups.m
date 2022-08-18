function [probdist, torquewordcolumn] = torquebreakups(numofbin, torqvec)



        choppedtorquevec1 = zeros(1, length(torqvec(:, 1)));
        choppedtorquevec2 = zeros(1, length(torqvec(:, 2)));
        pc1range =length(torqvec(:, 1));
        pc1bin = floor(pc1range./numofbin);
        pc2range = length(torqvec(:, 2));
        pc2bin = floor(pc2range./numofbin);
pc1range = [];
pc2range =[];
torquewordcolumn = zeros(1, length(torqvec(:, 1)));
sortedpc1 = sort(torqvec(:, 1));
sortedpc2 = sort(torqvec(:, 2));

 for y = 1:numofbin
            pc1range = [pc1range, sortedpc1(pc1bin*y)];
            pc2range = [pc2range, sortedpc2(pc2bin*y)];
 end

 if pc1range(end) ~= max(torqvec(:, 1))
     pc1range(end) = max(torqvec(:, 1));
     pc2range(end) = max(torqvec(:, 2));
     
 end

 %lines 29 through 48 split the torque columns into 2 bins each
        choppedtorquevec1(torqvec(:, 1)<=pc1range(1)) = 1;

        choppedtorquevec2(torqvec(:, 2)<=pc2range(1)) = 1;

        if length(pc1range) > 1
            for y = 2:length(pc1range)
                    
                mask1 = torqvec(:, 1)<=pc1range(y);
                mask2 = torqvec(:, 1)> pc1range(y-1);
                mask3 = (mask1+mask2 == 2);
                choppedtorquevec1(mask3) = y;
                mask1 = torqvec(:, 2)<=pc2range(y);
                mask2 = torqvec(:, 2)> pc2range(y-1);
                mask3 = (mask1+mask2 == 2);
                choppedtorquevec2(mask3) = y;
                

            end
        end
       combinedcolumn = [choppedtorquevec1', choppedtorquevec2'];
probdist = [];
        value = 1;%this loop sets the combination of torque bins into a single column of torque words
        for y = 1:max(choppedtorquevec1)
            for x = 1:max(choppedtorquevec2)
                mask1 = choppedtorquevec1 == y;
                mask2 = choppedtorquevec2 == x;
                mask3 = (mask1+mask2 == 2);
                torquewordcolumn(mask3) = value;
                probdist = [probdist, sum(mask3)/length(torquewordcolumn)];
                value = value+1;
            end
        end


end



