clc
clear all
%%creating paramaters
%
x=6;

% police_n = 64;

crowd_n = 1120;
police_n = 118;
L = 0.82;
T = 0.1;  


age_max = 120;

grid_size = 40; 
vision_l =6;
P_Vis = 6;
smoke_size = 6;
smoke_L = 40;
no=100  ;

%number of iterations
jail = {[9,9,9,-1,-1]};
jail_L=15 ;
%legitimacy of goverment


% T is a value which if G - N > T a active begins to riot (xpected utility
% of not expressing there private grievance)

JustMakeGraphXD = RiotW_Smoke_VS_Without(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,smoke_L,P_Vis,jail,jail_L);
% % [group_green, group_blue, police] = create_ethnic(crowd_n,police_n,age_max);
% % grid = rand_initiallize_grid_ethnic(crowd_n, group_green, group_blue, police_n, police, grid_size);
% % grid = ethnic_agents(grid,grid_size,L,vision_l,T,P_Vis,age_max);
% % [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis);
% % 
% % [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis);
% % [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis);
% % 
% % jail_keys4 = count_jail_keys(jail,4)
% % jail_keys1 = count_jail_keys(jail,1)
% % one = count_keys(grid,1)
% % four = count_keys(grid,4)


% [group_green, group_blue, police] = create_ethnic(crowd_n,police_n,age_max);
% grid = rand_initiallize_grid_ethnic(crowd_n, group_green, group_blue, police_n, police, grid_size);
% grids = {};
% active_keys = {};
% gridkeys = {};
% purpleno = zeros(1,100);
% greenno = zeros(1,30);
% 
% jailno = zeros(1,100);
% 
% 
% purpleno(1,1) = count_keys(grid,1);
% greenno(1,1) = count_keys(grid,4);
% policeno(1,1) = count_keys(grid,2);
% 
% 
% grid_keys(grid,grid_size);
% grids = {grid};
% for i = 1:300  
%     i
% grid = ethnic_agents(grid,grid_size,L,vision_l,T,P_Vis,age_max);
% [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis);
% purpleno(1,i+1) = count_keys(grid,1);
% greenno(1,i+1) = count_keys(grid,4);
% policeno(1,i+1) = count_keys(grid,2);
% jailno(1,i+1) = length(jail);
% 
% grids{end+1} = grid;
% 
% end
% 
% plot(purpleno,'DisplayName','purple no')
% hold on 
% plot(greenno, 'DisplayName','green no')
% plot(jailno, 'DisplayName', 'jail no')
% xlabel('Itterations')
% ylabel('Number of Agents')
% hold off
% 
% legend

% for i = 1:30
% active_keys{end+1} = activekeys(grids{i},grid_size);
% 
% gridkeys{end+1} = grid_keys(grids{i},grid_size);
% 
% end
% 
% A = {}
% for i = 1:30
%     gridkeys{i}
%  A{end+1} = [gridkeys{i} active_keys{i}];
% end



% filename = 'ethnic_cleanse_keys.xlsx';
% writematrix(gridkeys,filename)
% 
% filename2 = 'ethnic_cleanse_Actives.xlsx';
% writematrix(active_keys,filename2)





% [grids,agentno,jailno,riotno] = run_till_no_nosmoke(grid_size,L,vision_l,T,crowd_n, police_n, no,jail,jail_L,P_Vis,2);
% plot(riotno)
% plot(sriotno,'displayname','Rioters with smoke')
% hold on
 
% hold off
% legend
% sgrid = smoke_grid_init(grid_size);





%  smokeriotno = run_till_no_smoke(grid_size,L,vision_l,T,crowd_n, police_n,no,sgrid,smoke_size);
%  riotno = run_till_no_nosmoke(grid_size,L,vision_l,T,crowd_n, police_n, no);
%  plot(1:length(riotno),riotno,1:length(smokeriotno),smokeriotno)




%%my functions
%{

run through a riot until all rioters have been arested
%}

function [grid,riotno,jailno] = effect_of_Police(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,smoke_L,P_Vis,jail)
    jail_L=10000000;
    %create the grid
    [crowd,police] = create_crowd_police(crowd_n,police_n);
    grid  = rand_initiallize_grid(crowd_n, crowd, police_n, police, grid_size);
       riotno = zeros(no,1);
       jailno = zeros(no,1);
       policeno = zeros(no,1);
       policeno(1,1) = count_keys(grid,2);
       riotno(1,1) = 0;
       jailno(1,1) = 0;
       count = 2;
    for i = 1:no
        grid = agents_nosmoke(grid,grid_size,L,vision_l,T);
        riotno(count,1)=count_keys(grid,3);

        grid = removeC(grid,grid_size);
        [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis); 
        
        

        
        policeno(count,1) = count_keys(grid,2);
        jailno(count,1) = length(jail);
        count = count +1;
    end
    plot(riotno)
    hold on 
    plot(jailno)
    plot(policeno)
    legend({'Number Of Rioters','Number of Prisoners','Number of Police'},'Location','east')
    xlabel('Itteration Number')
    ylabel('Number of Agents')
    hold off
end

 

function [Ls,riotno,jailno] = effect_of_legitamcy(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,smoke_L,P_Vis,jail)
    jail_L = 10000000000;
    Ls = 0.90:-0.01:0.0;
    %create the grid
    [crowd,police] = create_crowd_police(crowd_n,police_n);
    grid  = rand_initiallize_grid(crowd_n, crowd, police_n, police, grid_size);
       riotno = zeros(length(Ls),1);
       jailno = zeros(length(Ls),1);
       riotno(1,1) = 0;
       jailno(1,1) = 0;
       count = 2;
    for L = Ls
        grid = agents_nosmoke(grid,grid_size,L,vision_l,T);

        [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis); 
        
        

        riotno(count,1)=count_keys(grid,3);
        
        jailno(count,1) = length(jail);
        count = count +1;
    end
    plot(riotno)
    hold on 
    plot(jailno)
    plot(Ls*1000)
    legend({'Number Of Rioters','Number of Prisoners','Legitamacy'},'Location','north')
    xlabel('Itteration Number')
    ylabel('Number of Agents/Legitimacy*1000')
    hold off
end

function [Ls,riotno,jailno] = sudden_effect_of_legitamcy(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,smoke_L,P_Vis,jail,jail_L)
    Ls =[0.9*ones(1,77),0.6*ones(1,42)];    
    %create the grid
    [crowd,police] = create_crowd_police(crowd_n,police_n);
    grid  = rand_initiallize_grid(crowd_n, crowd, police_n, police, grid_size);
       riotno = zeros(length(Ls),1);
       jailno = zeros(length(Ls),1);
       riotno(1,1) = 0;
       jailno(1,1) = 0;
       count = 2;
    for i = 1:77
        L = 0.9;
        %grid = agents_nosmoke(grid,grid_size,L,vision_l,T);

        [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis); 
        
        

        riotno(count,1)=count_keys(grid,3);
        
        jailno(count,1) = length(jail);
        count = count +1;
    end
    for i = 78:120
        L = 0.7;
        grid = agents_nosmoke(grid,grid_size,L,vision_l,T);
        riotno(count,1)=count_keys(grid,3);
        [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis); 
        
        

        
        
        jailno(count,1) = length(jail);
        count = count +1;

    end
    plot(riotno)
    hold on 
    plot(jailno)
    plot(Ls*1000)
    legend({'Number Of Rioters','Number of Prisoners','Legitamacy'},'Location','northeast')
    xlabel('Itteration Number')
    ylabel('Number of Agents/Legitimacy*1000')

    hold off


end




function [sgrid,riotno] = run_till_no_smoke(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,smoke_L,P_Vis,jail,jail_L,type)
    %create the grid
    if type ==1
        [crowd,police] = create_crowd_police(crowd_n,police_n);
        grid = place_riot(crowd_n, crowd, police_n, police, grid_size);
    end
    if type ==2 
        [crowd,police] = create_crowd_police(crowd_n,police_n);
        grid  = rand_initiallize_grid(crowd_n, crowd, police_n, police, grid_size);
    end
    %create the smoke grid
    sgrid = smoke_grid_init(grid_size);
    %make agents rioters
    grid = agents_smoke(grid,grid_size,L,vision_l,T,sgrid);
    riotno = zeros(no,1);
    agentno = zeros(no,1);

    count = 1;
    riotno(count,1)=count_keys(grid,3);
    
    for i = 1:no
        
        count = count+1;
        
        
        [jail, sgrid, grid] = policemen_smoke(grid, grid_size,smoke_size,smoke_L,sgrid,jail,jail_L,P_Vis);
        sgrid
        grid = agents_smoke(grid,grid_size,L,vision_l,T,sgrid);
        
        riotno(count,1)=count_keys(grid,3);
        agentno(count-1,1) = count_keys(grid,1);
        
        
        
    end
    plot(riotno)
end



function riotno = run_till_no_nosmoke_smoke(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,jail,jail_L)
    %create the grid
    [crowd,police] = create_crowd_police(crowd_n,police_n);
    grid  = rand_initiallize_grid(crowd_n, crowd, police_n, police, grid_size);
    %create the smoke grid
    sgrid = smoke_grid_init(grid_size);
    %make agents rioters
    grid = agents_smoke(grid,grid_size,L,vision_l,T,sgrid);
    riotno = zeros(no,1);
    agentno = zeros(100,1);

    count = 1;
    riotno(count,1)=0;
    
    for i = 1:no/2
        
        count = count+1;
        

        grid = policemen(grid, grid_size,jail,jail_L);
        grid = agents_nosmoke(grid,grid_size,L,vision_l,T);

        riotno(count,1)=count_keys(grid,3);
        agentno(count-1,1) = count_keys(grid,1);
        
    end
    for i = 1:no/2
        
        count = count+1;
        
        
        [jail, sgrid, grid] = policemen_smoke(grid, grid_size,smoke_size,smoke_L,sgrid,jail,jail_L,P_Vis);
        
        grid = agents_smoke(grid,grid_size,L,vision_l,T,sgrid);
        
        riotno(count,1)=count_keys(grid,3);
        agentno(count-1,1) = count_keys(grid,1);
        
        
        
    end
    
end



function [grids,agentno,jailno,riotno] = run_till_no_nosmoke(grid_size,L,vision_l,T,crowd_n, police_n, no,jail,jail_L,P_Vis,type)
    if type ==1
        [crowd,police] = create_crowd_police(crowd_n,police_n);
        grid = place_riot(crowd_n, crowd, police_n, police, grid_size);
    end
    if type ==2 
        [crowd,police] = create_crowd_police(crowd_n,police_n);
        grid  = rand_initiallize_grid(crowd_n, crowd, police_n, police, grid_size);
    end
    
    %grid = agents_nosmoke(grid,grid_size,L,vision_l,T);
    riotno = zeros(no,1);
    agentno = zeros(no,1);
    jailno = zeros(no,1);
    count = 1;
    grids = {grid};
    for i = 1:no
        
        
        [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis); 
        
        grid = agents_nosmoke(grid,grid_size,L,vision_l,T);
        
        
        riotno(count,1)=count_keys(grid,3);
        agentno(count,1) = count_keys(grid,1);
        jailno(count,1) = length(jail);
        grids{end+1} = grid;
        count = count+1;
        
%         numciv = count_keys(grid,1)
%         numriot = count_keys(grid,3)
%         total = numciv+numriot
        
    end
    plot(jailno)
    plot(riotno)
end
    
function plot = run_ethnic(crowd_n, group_green, group_blue, police_n, police, grid_size,age_max,L,vision_l,T,P_Vis,jail,jail_L)
    [group_green, group_blue, police] = create_ethnic(crowd_n,police_n,age_max);
    grid = rand_initiallize_grid_ethnic(crowd_n, group_green, group_blue, police_n, police, grid_size);
    grids = {};
    active_keys = {};
    gridkeys = {};
    purpleno = zeros(1,100);
    greenno = zeros(1,30);
    
    jailno = zeros(1,100);
    
    
    purpleno(1,1) = count_keys(grid,1);
    greenno(1,1) = count_keys(grid,4);
    policeno(1,1) = count_keys(grid,2);
    
    
    grid_keys(grid,grid_size);
    grids = {grid};
    for i = 1:300  
        i
    grid = ethnic_agents(grid,grid_size,L,vision_l,T,P_Vis,age_max);
    [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis);
    purpleno(1,i+1) = count_keys(grid,1);
    greenno(1,i+1) = count_keys(grid,4);
    policeno(1,i+1) = count_keys(grid,2);
    jailno(1,i+1) = length(jail);
    
    grids{end+1} = grid;
    
    end
    
    plot(purpleno,'DisplayName','purple no')
    hold on 
    plot(greenno, 'DisplayName','green no')
    plot(jailno, 'DisplayName', 'jail no')
    xlabel('Itterations')
    ylabel('Number of Agents')
    hold off
    
    legend
end



function extinction_age = run_until_extinction(crowd_n,police_n,age_max,grid_size,L,vision_l,T,P_Vis,jail,jail_L)
    p=0;
     
    [group_green, group_blue, police] = create_ethnic(crowd_n,police_n,age_max);
    grid = rand_initiallize_grid_ethnic(crowd_n, group_green, group_blue, police_n, police, grid_size);
    for i = 1:15000
        grid = ethnic_agents(grid,grid_size,L,vision_l,T,P_Vis,age_max);
        [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis);
        if or(and(count_keys(grid,1) ==0, count_jail_keys(jail,1)==0), and(count_keys(grid,4)==0,count_jail_keys(jail,4) == 0))
            extinction_age = i;
            p=1;
            break
        end
    end
    if p==0
    extinction_age = 15000
    end
end

function extinction_age = affect_police_ethnic(crowd_n,age_max,grid_size,L,vision_l,T,P_Vis,jail,jail_L)
        police_n = 0:3:150;
        extinction_age = zeros(50,50);
        run_until_extinction(crowd_n,0,age_max,grid_size,L,vision_l,T,P_Vis,jail,jail_L);

        for j = police_n
            for i = 1:50
                j
                i
                extinction_age(i,(j/3)+1) = run_until_extinction(crowd_n,i,age_max,grid_size,L,vision_l,T,P_Vis,jail,jail_L);
                extinction_age(i,(j/3)+1)
            end
        end
end

function jorwi= random_place_riot_place(grid_size,L,vision_l,T,crowd_n, police_n, no,jail,jail_L,P_Vis)
        %riot
        jail_L = 10000;
        [grids,agentno,jailno,riotno]=run_till_no_nosmoke(grid_size,L,vision_l,T,crowd_n, police_n, no,jail,jail_L,P_Vis,1)
        %rand
        [grids,ragentno,rjailno,rriotno]=run_till_no_nosmoke(grid_size,L,vision_l,T,crowd_n, police_n, no,jail,jail_L,P_Vis,2)
        plot(riotno,'DisplayName','Riot Placement Rioter')
        hold on
        plot(rriotno,'DisplayName','Random PLacement Rioter')
        plot(jailno,'DisplayName','Riot Placement Jail no')
        plot(rjailno,'DisplayName','Random Placement Jail no')
        ylabel('Number of Agents')
        xlabel('Itteration Number') 
        hold off
        legend
end

function JustMakeGraphXD = RiotW_Smoke_VS_Without(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,smoke_L,P_Vis,jail,jail_L)
    jail_L = 9999999;
        [grids,agentno,jailno,riotno] = run_till_no_nosmoke(grid_size,L,vision_l,T,crowd_n, police_n, no,jail,jail_L,P_Vis,1);
         [sgrid,Priotno] = run_till_no_smoke(grid_size,L,vision_l,T,crowd_n, police_n, no,smoke_size,smoke_L,P_Vis,jail,jail_L,1);
        plot(riotno,'DisplayName','Riot No Without Smoke')
        hold on 
        plot(Priotno,'DisplayName','Riot No With Smoke')
        hold off
        legend
        JustMakeGraphXD = 0;
end
%% agents
function grid = agents_smoke(grid,grid_size,L,vision_l,T,sgrid)
    grid = move_rand(1, grid_size, grid);
    grid = move_rand(3,grid_size, grid);
    grid = agent_to_riot_smoke(grid,grid_size,L,vision_l,T,sgrid);
end

function grid = agents_nosmoke(grid,grid_size,L,vision_l,T)
    grid = move_rand(1, grid_size, grid);
    grid = move_rand(3,grid_size, grid);
    grid = agent_to_riot_nosmoke(grid,grid_size,L,vision_l,T);
end


function [jail, sgrid, grid] = policemen_smoke(grid, grid_size,smoke_size,smoke_L,sgrid,jail,jail_L,P_Vis)
    grid = move_rand(2, grid_size, grid);
    sgrid = smoke(grid,grid_size,smoke_size,sgrid,smoke_L);
    sgrid = updatesmoke(sgrid);
    [jail,grid] = riot_to_arrest(grid_size, grid,jail,jail_L,P_Vis);
    [jail,grid] = update_jail(jail,grid,grid_size);
end



function [jail,grid] = policemen(grid, grid_size,jail,jail_L,P_Vis)
    grid = move_rand(2, grid_size, grid);
    [jail,grid] = update_jail(jail,grid,grid_size);
    [jail,grid] = riot_to_arrest(grid_size, grid,jail,jail_L,P_Vis);
    
    
end


function grid = ethnic_agents(grid,grid_size,L,vision_l,T,P_Vis,age_max)
    grid = move_rand(1, grid_size, grid);
    grid = move_rand(4,grid_size, grid);
    grid = ethnic_group_to_riot(grid,grid_size,L,vision_l,T);
    grid = kill(grid,grid_size,P_Vis);
    grid = cloneing(grid_size,grid,age_max);
    grid = ageing(grid,grid_size);
end


%% model 1

%{
grievance():  
    -Hardship G (physical or economic deprivation) drawn from U(0,1) ]
    uniform distrobution on 0,1
    
    -Legitamacy L percived legitamacy of a regime. This will be arbritary
    number from 0,1
%}
function G  = grievance(H,L)
    a = 1-L;
    G = H*a;
end

%{
arrest_probability():
    -constant k set such that P = 0.9 when C=A=1
    
    -active rioters A

    -cops in the area C
C/A changes for every rioter depending on how many C's and A's are within
the vision of each crowd member within a set vision
%}
function P = arrest_probability(k,C,A)
    P = 1 - exp(-(k*(C/A)));
end


%{
create_crowd_police()
-crowd_n:   number of people in crowd
-police_n:  number of police

-crowd:     our crowd with h r vals and key 1
-police:    our police with h r vals and key 2

creates the crowd aswell as the police such that:
-in the crowd each member has a number 1 in row1 to represent that they
are non rioters aswell as a hardship (H) row2 and a level of
risk they are willing to take (R) in row3. form: 
[0,h1,r1; 0,h2,r2; ...; 0,hn,rn]

-in the police we have our key as 2 instead of 1 and we have H=R=0 as all
the police do is arrest and have no R H 
%}
function [crowd,police] = create_crowd_police(crowd_n,police_n)
    %take our H and R randomly from the uniform distrobution beween 0-1
    H = unifrnd(0,1,1,crowd_n)';
    R = unifrnd(0,1,1,crowd_n)';
    %make a list of 1's from each person in the crowd
    crowdkey = ones(1,crowd_n)';
    %put crowd in the form 
    crowd = [crowdkey,H, R, zeros(1,crowd_n)'];
    %same but police
    police = [2*ones(1,police_n)',zeros(1,police_n)',zeros(1,police_n)',zeros(1,police_n)'];
end

%{
rand_initiallize_grid()
-crowd_n:   number of people in crowd
-crowd:     from create_crowd_police()
-ploice_n:  number of police
-police:    from create_crowd_police
-grid_size: how big we want our grid

-grid:      our initialized grid
%}

function grid  = rand_initiallize_grid(crowd_n, crowd, police_n, police, grid_size)

    %first start our grid out as being a grid of zeros for our given grid
    %size
    for i = 1:grid_size
        for j = 1:grid_size
            grid{i,j} = [0,0,0,0];
        end
    end
    
    
        
%     now we put in our crowd, distrobuted randomly
%     throughout the grid
    while crowd_n ~= 0

        i = randi(grid_size,1);
        j = randi(grid_size,1);
        if grid{i,j}(1) == 0
            grid{i,j}=crowd(crowd_n,1:4);
            crowd_n = crowd_n - 1;
        end
    end
    
    
    %now we put the desired number of plolice randomly into our grid    
    while police_n ~=0
        i = randi(grid_size,1);
        j = randi(grid_size,1);
        if grid{i,j}(1) == 0
            grid{i,j}=police(police_n,1:4);
            police_n = police_n - 1;
        end
    end
end

%{
grid:       the grid
grid_size:  size of the grid
vision_l:   length of the vision

v:          vision at point i j

takes in the grid and a point of the grid, it then returns a subsection
ofthe grid around this point
%}
function v = vision(i,j,grid, vision_l,grid_size)

    w = i - vision_l;
    e = i + vision_l;
    s = j + vision_l;
    n = j - vision_l;
    if w < 1
        w = 1;
    end

    if e > grid_size
        e = grid_size;
    end

    if s >grid_size
        s = grid_size;
    end

    if n < 1
        n = 1;
    end
    dns = s-n+1;
    dew = e-w+1;
    v = grid(w:e,n:s);

    %v = reshape({grid(w:e,n:s)},[dew,dns]);
    
end



%{


takes in the grid and checks the number of rioters around each agent with
the number of police and uses equations to evaluate if the agent becomes a 
rioter aswell
%}


function grid = agent_to_riot_nosmoke(grid,grid_size,L,vision_l,T)
    %loop through each space in the grid
    for i = 1:grid_size
    for j = 1:grid_size
        %if we land on an agent
        if grid{i,j}(1)==1 || grid{i,j}(1) == 3
            %find its grievance number
            G  = grievance(grid{i,j}(2),L);

            %find its risk probablility
            R  = grid{i,j}(3);

            %find the vision at that point
            v = vision(i,j, grid, vision_l, grid_size);
            
            %count number of police in vision
            C=count_keys(v,2);

            %count number of actives(rioters) in vision
            A=count_keys(v,3)+1;
            
            %with our C and A we can now find the probability for our guy
            %at i j to be arrested
            P = arrest_probability(2.3,C,A);
            
            %our agents net risk
            N = R*P;
            
            %expected utility of publicly expressing ones private grievanc
            gmn = G-N;
            
            if gmn >T
                grid{i,j}(1) = 3;
            else
                grid{i,j}(1) = 1;
            end

        end

    end
    end
    %grid;
end





%{

turns agents into rioters taken into acount the smoke
%}
function grid = agent_to_riot_smoke(grid,grid_size,L,vision_l,T,sgrid)
    %loop through each space in the grid
    for i = 1:grid_size
    for j = 1:grid_size
        %if we land on an agent
        if grid{i,j}(1)==1
            %find its grievance number
            G  = grievance(grid{i,j}(2),L);

            %find its risk probablility
            R  = grid{i,j}(3);

            %find the vision at that point
            v = vision(i,j, grid, vision_l, grid_size);
            
            %count number of police in vision
            C=count_keys(v,2);

            %count number of actives(rioters) in vision excluding those
            %with smoke
            A=count_keys_smoke(v,3,sgrid)+1;
            
            %with our C- and A we can now find the probability for our guy
            %at i j to be arrested
            P = arrest_probability(-log(0.1),C,A);
            
            %our agents net risk
            N = R*P;
            
            %expected utility of publicly expressing ones private grievanc
            gmn = G-N;
            
            if gmn >T
                grid{i,j}(1) = 3;
            else
                grid{i,j}(1) = 1;
            end

        end

    end
    end
    %grid;
end







%{
grid:   the matrix you want to count the keys in
key:    the key you want to count

n:      the count

counts the number of a given key in a given matrix
%}
function n = count_keys(grid,key)
    n=0;
    [r,c] = size(grid);
    for x = 1:r
    for y = 1:c
        grid{x,y};
        if grid{x,y}(1)==key
            n=n+1; 
        end

    end
    end
end

%{

count keys in matrix excluding the ones with smoke on them
%}


function n = count_keys_smoke(grid,key,sgrid)
    n=0;
    [r,c] = size(grid);
    for x = 1:r
    for y = 1:c
            
        if grid{x,y}(1)==key && sgrid{x,y}(1) == 0
            n=n+1; 
        end

    end
    end
end

%{
grid_size:  size of grid
grid:       current grid

grid:       updated grid

looks at every police person and arrestsa rioter if it is near by
%}

function [jail,grid] = riot_to_arrest(grid_size, grid,jail,jail_L,P_Vis)
    %for every (i,j)th place on the grid
    
    for i = 1:grid_size
        for j = 1:grid_size
            %if the (i,j)th place is a policeman
            one = 0 ;
            
                if grid{i,j}(1) == 2 
                %check every square 3*3 around the police man
                    
                
                        for n = -P_Vis:P_Vis   
                            for m = -P_Vis:P_Vis
                                %assighn new (i,j)th position we are considering
                                positioni = i + n;
                                positionj = j + m;
                                %check point is on grid and it is a rioter
                                if positioni<=grid_size && positionj <= grid_size && positioni > 0 && positionj > 0 &&  grid{positioni,positionj}(1) == 3
                                    %arrest that rioter
                                   

                                   % grid{positioni,positionj}(1)=1;
                                     l = grid{positioni,positionj};
                                     sentence = randi(jail_L);
                                    l(4)=sentence;
                                    
                                   
                                    jail{end+1} = l;
                                    
                                    grid{positioni,positionj} = [0,0,0,0];
                                    %return so we only arest 1 person per turn
                                    one =1;
                                    break
                                    
                                end
                                
                            end
                            if one ==1 
                                break
                            end
                        
                        end
                           
            end
        end
    end


    
end








%{
crowd_trype:    the type of crowd you want to move
grid_size:      the size of the grid
grid:           the current grid

grid:           the grid after all the agents of a certain type have moved
%}

function grid = move_rand(crowd_type, grid_size, grid)
    v = 0;
    %for every (i,j)th place on the grid
    for i = 1:grid_size
        for j = 1:grid_size
            %while we reach a place with the correct crowd type, and we
            %have made a certain number of attempts v
            while grid{i,j}(1) == crowd_type && v < 1000
                %create a random number -1 to 1 for n-s plane and w-e plane
                %and take the absolute value to get rid of negatives
                directionns = i + randi([-5,5],1);
                directionwe = j + randi([-5,5],1);
                v = v + 1;
                %check if the new position is inside of the grid and not =0
                if directionns<=grid_size && directionwe <= grid_size && directionns > 0 && directionwe > 0
                    %if that space is empty (=0)
                    if grid{directionns,directionwe}(1) == 0
                        %move that person to that location
                        grid{directionns,directionwe} = grid{i,j};
                        %remove that person from where he was previously
                        grid{i,j} =  [0,0,0,0,0];
                    end
                end
            end
        end
    end

end



%%smoke section

function sgrid = smoke(grid,grid_size,smoke_size,sgrid,smoke_L)
    %go through each square on grid
    for i = 1:grid_size
        for j = 1:grid_size
            
           
           if grid{i,j}(1) == 2
               %throw smoke at fools
               sgrid = throwsmoke(i,j,grid,smoke_size,grid_size,sgrid,smoke_L);
           end
        end
    end
end




function sgrid = throwsmoke(i,j,grid,smoke_size,grid_size,sgrid,smoke_L)
    L = length(smoke_size);
    %search around the police officer for rioters
    for n = -L:L
        for m = -L:L
            %assighn new (i,j)th position we are considering
            positioni = i + n;
            positionj = j + m;
            
     
            %check point is on grid and it is a rioter
            if positioni<=grid_size && positionj <= grid_size && positioni > 0 && positionj > 0 && grid{positioni,positionj}(1) == 3
                %create an square area of length vision_l around our rioter
                w = positioni - smoke_size;
                e = positioni + smoke_size;
                s = positionj + smoke_size;
                north = positionj - smoke_size;
                %make each point inside of the grid
                if w < 1
                    w = 1;
                end
            
                if e > grid_size
                    e = grid_size;
                end
            
                if s >grid_size
                    s = grid_size;
                end
            
                if north < 1
                    north = 1;
                end
                
                for k = w:e
                    for o = north:s
                        
                        sgrid{k,o}(1) = 1;
                        sgrid{k,o}(2) = smoke_L;
                    end
                end
                
                
            end
        end
    end
end




function sgrid = smoke_grid_init(grid_size)
    sgrid = {};
    for i = 1:grid_size
        for j = 1:grid_size
            sgrid{i,j} = [0,0];
        end
    end
end




function sgrid = updatesmoke(sgrid)
    N = length(sgrid);
    for i = 1:N
        for j = 1:N
            if sgrid{i,j}(2)~=0
                sgrid{i,j}(2) = sgrid{i,j}(2) -1;
            end
            if sgrid{j,i}(2) ==0 && sgrid{j,i}(1)==1
                sgrid{j,i}(1)=0;
            end
        end
    end
end


function grid = removeC(grid,grid_size)
    n = 0;
    
        for i = 1:grid_size
            for j = 1:grid_size
                if grid{i,j}(1) == 2
                    grid{i,j} = [0,0,0,0];
                    n = n+1;
                    if n == 11
                        return
                    end
                end
                
            end
            
        end
    
    
end

function grid = place_riot(crowd_n, crowd, police_n, police, grid_size,grid)
    for i = 1:grid_size
        for j = 1:grid_size
            grid{i,j} = [0,0,0,0];
        end
    end
    
        for i = 1:grid_size
            for j = 1:grid_size
                grid{i,j} = crowd(crowd_n,1:4);
                if crowd_n == 1
                    grid = place_riotsquad(police_n, police, grid_size,grid);
                    return
                end
                crowd_n = crowd_n -1;
            end
        end
    
end


function grid = place_riotsquad(police_n, police, grid_size,grid)
    for i = grid_size:-1:1
            for j = grid_size:-1:1
                grid{i,j} = police(police_n,1:4);
                if police_n == 1
                    return
                end
                police_n = police_n -1;
            end
    end
end





function gridkeys = grid_keys(grid,grid_size)
    gridkeys = zeros(grid_size,grid_size);
    for i = 1:grid_size
        for j = 1:grid_size
            gridkeys(i,j) = grid{i,j}(1);
        end
    end
end

function activekeys = activekeys(grid,grid_size)
    activekeys = zeros(grid_size,grid_size);
    for i = 1:grid_size
        for j = 1:grid_size
            activekeys(i,j) = grid{i,j}(4);
        end
    end
end
function jail_keys = count_jail_keys(jail,key)
    L = length(jail);
    jail_keys = 0;
    for i = 1:L
        if jail{i}(1) == key
            jail_keys = jail_keys +1;
        end
    end
end


%%
%{
this is the section for the ethnic model part

I need to:
 - create the groups
 - put groups on board
 - make it so you can distinguish between different groups being active
 - let groups murder each other
 - cloneing
 - old age death
%}


%{
this function takes in the number of police and people in each ethnic group
and gives back them,

%}


function [group_green, group_blue, police] = create_ethnic(crowd_n,police_n,age_max)
    %take our H and R randomly from the uniform distrobution beween 0-1
    Hgreen = unifrnd(0,1,1,crowd_n)';
    Rgreen = unifrnd(0,1,1,crowd_n)';

    Hblue = unifrnd(0,1,1,crowd_n)';
    Rblue = unifrnd(0,1,1,crowd_n)';
    %make a list of 1's from each person in the crowd
    crowdkey_blue = ones(1,crowd_n)';
    crowdkey_green = 4*ones(1,crowd_n)';
    %put crowd in the form 
    group_green = [crowdkey_blue,Hgreen, Rgreen, zeros(1,crowd_n)',randi([0,age_max],1,crowd_n)'];

    group_blue = [crowdkey_green,Hblue,Rblue,zeros(1,crowd_n)',randi([0,age_max],1,crowd_n)'];


    %same but police
    if police_n ==0 
        police = [];

    else
        police = [2*ones(1,police_n)',zeros(1,police_n)',zeros(1,police_n)',zeros(1,police_n)',zeros(1,police_n)'];

    end
end

%{
randomly place all agents on grid for ethnic group simulations
%}

function grid  = rand_initiallize_grid_ethnic(crowd_n, group_green, group_blue, police_n, police, grid_size)

    %first start our grid out as being a grid of zeros for our given grid
    %size
    for i = 1:grid_size
        for j = 1:grid_size
            grid{i,j} = [0,0,0,0,0];
        end
    end
    
    
        
%     now we put in our crowd, distrobuted randomly
%     throughout the grid
    p = crowd_n;
    while p ~= 0

        i = randi(grid_size,1);
        j = randi(grid_size,1);
        if grid{i,j}(1) == 0
            grid{i,j}=group_green(crowd_n,1:5);
            p = p - 1;
        end
    end

    while crowd_n ~= 0

        i = randi(grid_size,1);
        j = randi(grid_size,1);
        if grid{i,j}(1) == 0
            grid{i,j}=group_blue(crowd_n,1:5);
            crowd_n = crowd_n - 1;
        end
    end
    
    
    %now we put the desired number of plolice randomly into our grid    
    while police_n ~=0
        i = randi(grid_size,1);
        j = randi(grid_size,1);
        if grid{i,j}(1) == 0
            grid{i,j}=police(police_n,1:5);
            police_n = police_n - 1;
        end
    end
end


%{
now I code the murder
%}

function grid = kill(grid,grid_size,P_Vis)
    %for every (i,j)th place on the grid
    
    for i = 1:grid_size
        for j = 1:grid_size
            %if the (i,j)th place is 
            two= 0;
            
            if grid{i,j}(1) == 1 && grid{i,j}(4) == 1
            %check every square 3*3 around the police man
                
            
                    for n = -2:2   
                        for one = -2:2
                            %assighn new (i,j)th position we are considering
                            positioni = i + n;
                            positionj = j + one;
                            %check point is on grid and it is a rioter of
                            %the opposite group
                            if positioni<=grid_size && positionj <= grid_size && positioni > 0 && positionj > 0 && grid{positioni,positionj}(1) == 4
                                %remove from grid :)
                                grid{positioni,positionj} = [0,0,0,0,0];
                                %return so we only arest 1 person per turn
                                two=1;
                                break
                                
                            end
                            
                        end
                        if two==1 
                            break
                        end
                    
                    end
                           
            else

                if grid{i,j}(1) == 4 && grid{i,j}(4) == 1
                %check every square 3*3 around the police man
                    
                
                        for n = -2:2   
                            for one = -2:2
                                two = 0;
                                %assighn new (i,j)th position we are considering
                                positioni = i + n;
                                positionj = j + one;
                                %check point is on grid and it is a rioter
                                if positioni<=grid_size && positionj <= grid_size && positioni > 0 && positionj > 0 && grid{positioni,positionj}(1) == 1
                                    %remove from grid :)
                                    grid{positioni,positionj} = [0,0,0,0,0];
                                    %return so we only arest 1 person per turn
                                    two=1;
                                    break
                                    
                                end
                                
                            end
                            if two==1 
                                break
                            end
                        
                        end
                               
                end
            end
        end
    end

end


function grid = ethnic_group_to_riot(grid,grid_size,L,vision_l,T)
    %loop through each space in the grid
    for i = 1:grid_size
    for j = 1:grid_size
        %if we land on an agent
        if  grid{i,j}(1) == 1 || grid{i,j}(1) == 4
            %find its grievance number
            G  = grievance(grid{i,j}(2),L);

            %find its risk probablility
            R  = grid{i,j}(3);

            %find the vision at that point
            v = vision(i,j, grid, vision_l, grid_size);
            
            %count number of police in vision
            C=count_keys(v,2);

            %count number of actives(rioters) in vision
            A=count_keys(v,3)+1;
            
            %with our C and A we can now find the probability for our guy
            %at i j to be arrested
            P = arrest_probability(2.3,C,A);
            
            %our agents net risk
            N = R*P;
            
            %expected utility of publicly expressing ones private grievanc
            gmn = G-N;
            
            if gmn >T
                grid{i,j}(4) = 1;
            else
                grid{i,j}(4) = 0;
            end

        end

    end
    end
    %grid;
end

function grid = cloneing(grid_size,grid,age_max)
    %loop throuh the grid
    for i = 1:grid_size
        for j = 1:grid_size
            %if the grid is one of the ethnic groups
            if grid{i,j}(1) == 1 || grid{i,j}(1) == 4
                %create a 1 in 20 chance to clone
                prob = randi(20);
                %if we get 1/20
                if prob == 1
                    %check spaces next to character being cloned
                    one = 0;
                    for x = -1:1
                        for y = -1:1
                            position_x = i + x;
                            position_y = j + y;
                            %check the position is within the grid and is
                            %empty
                            if position_x > 0 && position_x < grid_size+1 && position_y > 0 && position_y < grid_size + 1 &&  grid{position_x, position_y}(1) == 0
                                grid{position_x, position_y} = [grid{i,j}(1),grid{i,j}(2), unifrnd(0,1),0,randi([0,age_max])];
                                one = 1;
                                break
                            end
                        end
%                         if one ==1 
%                             break
%                         end
                    end
                end
                
            end
        end
    end
end


function grid = ageing(grid,grid_size)   
    %loop through the grid
    for i = 1:grid_size
        for j = 1:grid_size
            %check is someone ahs died of old age
            if grid{i,j}(5) == 0 && or(grid{i,j}(1)==1, grid{i,j}(1)==4)
                grid{i,j} = [0,0,0,0,0];
            %if theyre not dead make them age one time
            else
                if or(grid{i,j}(1)==1, grid{i,j}(1)==4)
                    grid{i,j}(5) = grid{i,j}(5) - 1;
                end
            end
        end
    end
end




function [jail,grid] = riot_to_arrest_ethnic(grid_size, grid,jail,jail_L,P_Vis)
    %for every (i,j)th place on the grid
    
    for i = 1:grid_size
        for j = 1:grid_size
            %if the (i,j)th place is a policeman
            one = 0 ;
            
                if grid{i,j}(1) == 2 
                %check every square 3*3 around the police man
                    
                
                        for n = -P_Vis:P_Vis   
                            for m = -P_Vis:P_Vis
                                %assighn new (i,j)th position we are considering
                                positioni = i + n;
                                positionj = j + m;
                                %check point is on grid and it is a rioter
                                if positioni<=grid_size && positionj <= grid_size && positioni > 0 && positionj > 0 &&  grid{positioni,positionj}(4) == 1
                                    %arrest that rioter
                                   

                                   % grid{positioni,positionj}(1)=1;
                                     l = grid{positioni,positionj};
                                     sentence = randi(jail_L);
                                    l(4)=sentence;
                                    
                                   
                                    jail{end+1} = l;
                                    
                                    grid{positioni,positionj} = [0,0,0,0,0];
                                    %return so we only arest 1 person per turn
                                    one =1;
                                    break
                                    
                                end
                                
                            end
%                             if one ==1 
%                                 break
%                             end
                        
                        end
                           
            end
        end
    end


    
end
