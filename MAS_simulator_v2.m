close all
clear all
clc
rng('default');   

% CNRS-LAAS-SNU MAS simulator
% 19 Jan 2021, 15 March 2021, 8 Jan 2024
% Shim / Tanwani
% for the paper

savefigures = 0;    % If 1, then save figures every time
showanimation = 1;  % If 1, show animation
saveanimation = 0;  % If 1, make a movie file


%% Setting environment

% Emulation of the jump times in the paper
SimTime = 20;   % Simulation Time (e.g., 40 secs)
epsilon = 0.01;    % the parameter epsilon (uncoupled case if epsilon == 0)
%epsilon = 0.1;
%epsilon = 1;
%epsilon = 0;
if epsilon ~= 0
    sigma1 = 0.5;
    sigma2 = 1.5;
    N0 = 2;
    tau = 0;        % The state tau, current value
    TS = [0];       % Storage for generated jump times: Time Stamps
    while (TS(end) < SimTime)
        r_min = max([1-tau,0])*epsilon/sigma2;
        r_max = (N0-tau)*epsilon/sigma1;
        NTI = r_min+(r_max-r_min)*rand;  % Next Time Interval
        TS = [TS, TS(end)+NTI];
    
        r_min = tau + NTI*sigma1/epsilon;
        r_max = min([N0,tau + NTI*sigma2/epsilon]);
        tau = r_min+(r_max-r_min)*rand - 1;
    end
    Njump = size(TS,2); % number of jumps
else    % uncoupled case
    TS = [0,SimTime];
    Njump = 2;
end

% Total number of agents
Nnode = 4; 

% Make the matrix W from the edges (1,2), (2,3), (3,4), (4,1)
E(1).W = Make_W(Nnode,1,2);
E(2).W = Make_W(Nnode,2,3);
E(3).W = Make_W(Nnode,3,4);
E(4).W = Make_W(Nnode,4,1);
Nedge = 4;  % number of edge

% Switching sequence lister
% for the ring graph used in the example
SS = [];    % Storage for switching sequence
Stmp = [];  % temporary storage for switching sequence
MaxStmp = 0;
while size(SS,2) < Njump
    Stmp = [];
    while sum(ismember([1:Nedge],Stmp)) ~= Nedge
        Stmp = [Stmp, randi(Nedge)];
    end     % completion of a connected graph by the edge sequence 
    MaxStmp = max([MaxStmp,size(Stmp,2)]);
    SS = [SS,Stmp];
end
disp(sprintf('Maximum length of the strings = %d',MaxStmp))


%% Setting the dynamic system & Simulate

% initial_condition: each row is the initial condition for each agent
initial_condition = [
     1, -2;
     2, 1;
     -1,2;
     -2,-1
     ];
%initial_condition = rand(Nnode,2);
X0 = initial_condition(:);

% parameter for the dynamics
mu = [-.1; 1; 2; 3];

% Simulate coupled MAS
Trecord = [];
Xrecord = [];
for i = 1:Njump-1
  %disp(sprintf('Simulating for time = %g', Tjump(i,1)))
  % run flow dynamics
  [t,x] = ode45( @(t,x) flow_map( t, x, Nnode, mu ), ...
      [TS(i), TS(i+1)], X0);    
  Trecord = [Trecord; t(1:end-1)];
  Xrecord = [Xrecord; x(1:end-1,:)];
  % run jump dynamics (uncoupled case when epsilon == 0) 
  if epsilon > 0
    X0 = jump_map( t(end), x(end,:)', Nnode, E(SS(i)).W )';
  end
end


%% Plot the simulation result

% Color Map
CM = {'k','b','r',[.2 .8 .6],[.5 .6 .7],'g','y'};

plot_time_span = 5;   % Second interval for increasing the number in figure

% Positioning figures
numFigures = 4;
for i = 1:numFigures
    h(i) = figure;
    pos = get(gcf,'Position');
    figWidth = pos(3);
    pos(1) = (i-1)*figWidth;
    set(gcf, 'Position', pos);
end

% Draw X
figure(h(1))
for j=1:Nnode
    plot(Trecord,Xrecord(:,j),'Color',CM{j})
    hold on
end
set(gca,'FontSize',15)
if savefigures == 1
    saveas(gcf,sprintf('plotx_%g.png',epsilon))
end

% Draw Y
figure(h(2))
for j=1:Nnode
    plot(Trecord,Xrecord(:,Nnode+j),'Color',CM{j})
    hold on
end
set(gca,'FontSize',15)
if savefigures == 1
    saveas(gcf,sprintf('ploty_%g.png',epsilon))
end

% Draw Phase Portrait
figure(h(3))
for j = 1:Nnode
    plot( Xrecord(:,j), Xrecord(:,Nnode+j), 'Color', CM{j} ) % Draw each sol
    hold on
end
count = 0;
for tim = 0:plot_time_span:Trecord(end)
  indexvector = find( Trecord > tim );
  indexfirst = indexvector(1);
  for j=1:Nnode
    text(Xrecord(indexfirst,j),Xrecord(indexfirst,Nnode+j),int2str(count),'Color',CM{j},'FontSize',14);
  end
  count = count+1;
end
if savefigures == 1
    saveas(gcf,sprintf('plotcoupled_%g.png',epsilon))
end

%% Animation
if showanimation == 1
    Animate(Xrecord,h(4),saveanimation)
end

return

%% Functions

%% Flow dynamics
% Note that the state x is a column vector
function dx = flow_map(t,x,Nnode,mu)

N = Nnode;

z = x(1:N);
y = x(N+1:2*N);

fx = mu.*(z.*z-1);
gx = y;
dx = [ -z+y; (1-fx).*(-z+y) - gx ];

end

%% Jump dynamics
% param is the structure containing some parameters coming from the workspace.
function xp = jump_map(t,x,Nnode,W)

N = Nnode;

z = x(1:N);
y = x(N+1:2*N);

xp = [ z; W*y ];

end

%% Making W 
function W = Make_W(Nnode,i,j)

W = diag(ones(Nnode,1));
W(i,i) = 0.5;
W(i,j) = 0.5;
W(j,i) = 0.5;
W(j,j) = 0.5;

end

%% Animation
function Animate(xx,fighandle,saveflag)

[a,b] = size(xx);

x1 = xx(:,1);
v1 = xx(:,5);
x2 = xx(:,2);
v2 = xx(:,6);
x3 = xx(:,3);
v3 = xx(:,7);
x4 = xx(:,4);
v4 = xx(:,8);

figure(fighandle)
plot(x1,v1,'r:')
hold on
plot(x2,v2,'b:')
plot(x3,v3,'g:')
plot(x4,v4,'m:')
axis([-3, 3, -4, 4])

if saveflag == 1
    v = VideoWriter('sync1.avi');
    %v.FileFormat = 'mp4';
    %v.FrameRate = 10;
    open(v);
end

speedfactor = round(a/30);
m = a;
for i=1:speedfactor:m,
    p1 = plot(x1(i),v1(i),'ro','MarkerSize',10,'MarkerFaceColor','r');
    p2 = plot(x2(i),v2(i),'bo','MarkerSize',10,'MarkerFaceColor','b');
    p3 = plot(x3(i),v3(i),'go','MarkerSize',10,'MarkerFaceColor','g');
    p4 = plot(x4(i),v4(i),'mo','MarkerSize',10,'MarkerFaceColor','m');
    drawnow;
    pause(0.1);

    if saveflag == 1
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    delete(p1);
    delete(p2); 
    delete(p3);
    delete(p4);
end

if saveflag == 1
    close(v)
end

end


