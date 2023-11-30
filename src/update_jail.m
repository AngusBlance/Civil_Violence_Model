

function [jail,grid] = update_jail(jail,grid,grid_size)
    %store the length of our jail in the beginning
    if isempty(jail) == 0
        %jail1 = jail
        %create a list for everyone leaving
        leavers = {};
        %create a list of people leaving
        
        for i = 1:length(jail)
            if jail{i}(4) == 0
               
                leavers{end+1} = jail{i};
                
            end
        end
        
        v = 1;
        %put people back onto the board
        count = length(leavers);
        while count ~= 0 && v ~= 100000
            i = randi(grid_size,1);
            j = randi(grid_size,1);
            v=v+1;
            if grid{i,j}(1) == 0
                %leavers{count}
                grid{i,j}= [leavers{count}];
                count = count - 1;
            end
            %assume grid is full if v = 10000 so empty the jail
            
        end
        %jail3 = jail
        %create count
        count = 1;
        %whilst there are people needing to leave the jail
        while count_leavers(jail) ~=0
            

            %remove them from the jail
          
            if jail{count}(4)==0
                jail(count)=[];
                leavers(1) = [];
                count = count -1;
            end
            count=count+1;
        end
        %jail4 = jail
        %jail = jail(not(cellfun(@(x)isequal(x([1]),[2]),jail)));
        N = length(jail);
        
        
        for i = 1:N
                
                
                jail{i}(4) = jail{i}(4) -1;
        end
    end
    
end






% 
% function [jail,grid] = update_jail1(jail,grid,grid_size)
%     %store the length of our jail in the beginning
%     if isempty(jail) == 0
%         %get the length of the jail
%         N = length(jail);
%         %for every ith element from 1-the length of the jail
%         for i = 1:N
%             %if the sentnce is =0
%             if jail{i}(4)==0
%                 q=1
%                 while q ~= 0
%                     %find a place on the grid
%                     for m =1:grid_size
%                         for n = 1:grid_size
%     
%                                 %if the space is free
%                                 if grid{n,m}(1)==0
%                                     %add the current ith jail member to
%                                     %that space
%                                     grid{n,m} = jail{i};
%                                     %change the first number of the guy to
%                                     %5
%                                     jail(i) = {[5,0.9,0.9,-1]};
%                                     %stop the while loop if we find a space
%                                     q=q-1
%                                     break
%                                     
%                                 end
%                                 break
%                             end
%                             break
%                     end
%                     break
%                     end
%                 end
%             
%             
%             end
%         i=1;
%         while count_leavers(jail) ~=0 && isempty(jail) == 0
%             
%                 if jail{i}(1)==5
%                     jail(i)=[];
%                 end
%            i=i+1;
%         end
% 
%         %jail = jail(not(cellfun(@(x)isequal(x([1]),[2]),jail)));
%         N = length(jail);
%         
%         
%         for i = 1:N
%                 
%                 
%                 jail{i}(4) = jail{i}(4) -1;
%         end
%     end
% end



function count = count_leavers(jail)
    n = length(jail);
    count = 0;
    for i = 1:n
        if jail{i}(4) == 0
            count = count +1;
        end
    end
end


