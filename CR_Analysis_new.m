% post velocity analysis 

nm='CR-2024-07-19-Exp-14';  % Change this 

% Array of experiment names
experiment_array = {'CR-2024-07-19-Exp-14'};

% Initialize a results table to store the output for each experiment
results = table('Size', [length(experiment_array), 9], ...
                'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
                'VariableNames', {'Experiment', 'cb_Pmax', 'cb_tmax', 'cb_Pmin', 'cb_tmin', 'cr_Pmax', 'cr_tmax', 'cr_Pmin', 'cr_tmin'});

% Loop through each experiment in the array
for idx = 1:length(experiment_array)

    nm = experiment_array{idx};  % Experiment name
    Nav=1;
    dir = AIV_names(nm);
    qt  = load(dir.mat.tank);

    % Find when experiment starts t1
    b = AIV_vstats(nm,Nav);
    t1=max(find(b.va.mean>.98*mean(b.va.mean(1:5)))); % Find the time where the table changes speed

    % Estimatge location of the coral
    R = sqrt(qt.Xw.^2+qt.Yw.^2);
    Q = atan2(qt.Yw,qt.Xw);
    r=linspace(0,50,100);
    j=0*r;
    for j=1:length(r)-1
      f=R>=r(j) & R<r(j+1);
      mxq(j)=max(b.svq(f));
      mnq(j)=mean(b.svq(f));
    end
    r(end)=[];
    f=find(mxq>max(mxq)/4 & r<45 & r>5);
    r1=min(r(f));  %change this if inner radius not right
    r2=max(r(f));  %change this if outer radius not right

    % Find the region where the azimuthal velocity is largest
    %Find the min and max radius where azimuthal velox is positive

    [maxValue, f] = max(b.svq(:));  % Find the maximum value and its index PARTH
    [i, j] = ind2sub(size(b.svq), f);  % Convert linear index to subscripts

    X=qt.Xw(j);
    Y=qt.Yw(i);
    clf;
    imagesc(max(b.svq,0));set(gca,'ydir','normal');
    hold('on');
    plot(j,i,'o','markerfacecolor',[0 0 0]);
    q = atan2(Y,X); %change if doesn't work
    %q=-2.1;
    dq = pi*(r2-r1)/(r1+r2)/2;
    q1=q-2*dq; %change this if ending angle not right
    q2=q+0*dq; %change this if starting angle not right
    msk = R>=r1 & R<=r2 & mod(Q-q1+pi,2*pi)-pi>=0 & mod(q2-Q+pi,2*pi)-pi>=0 ;

    %Check that red box is in the right location behind the coral
    % particle in a sector behind coral determined automatically
    p.rlim=[r1 r2];   
    p.qlim=[q1 q2];
    if isfield(p,'xlim') p=rmfield(p,{'xlim','ylim'});end
    p.msk=msk;
    p.ax=[-50 50];
    p.ay=[-50 50];
    p.framerate  = 0.99;

    p.t2=50; %maximum 50 frames for quiver
    figure(11)
    AIV_quiver(nm,Nav,p)
    while (true)
    txt = input("Is this region correct (Y/N)?","s")
    switch(txt)
        case('Y')
            disp('Region is correct')
            break;
     case('N')
          disp('Region is incorrect');
          disp('Edit r1, r2 (lines 24 & 25) and q1,q2 (lines 42 & 43)')
          error('Exiting')
        otherwise
            disp('Unknown response')
    end
    end

    % Make different files of trajectory data
    clear('p');
    force=false;
    p.topo=false;
    p.Np   = 1e4;

    p.framerate  = 0.99;
    p.t1=t1;
    Nav=1;
    %p.xlim = [-50.0 +50.0];p.ylim = [-50.0 +50.0];AIV_trajectories(nm,'1',Nav,p,force); %particles in the entire tank
    p.xlim = [ -0.1  +0.1];p.ylim = [-50.0  +0.0];AIV_trajectories(nm,'2',Nav,p,force); %particle on a line before coral
    %p.xlim = [-25.0 -10.0];p.ylim = [-20.0  -5.0];AIV_trajectories(nm,'3',Nav,p,force); %particles in a box behind coral

    % particle in a sector behind coral determined automatically
    p.rlim=[r1 r2];
    p.qlim=[q1 q2];
    if isfield(p,'xlim') p=rmfield(p,{'xlim','ylim'}); end
    p.msk=msk;
    p.ax=[-50 50];
    p.ay=[-50 50];
    p.framerate  = 0.99;
    %AIV_trajectories(nm,'4',Nav,p,force); 

    % particle in a sector 25% larger than the one determined automatically
    p5=p;
    p5.rlim = p.rlim+0.25*diff(p.rlim)*[-1 1];
    p5.qlim = p.qlim+0.25*diff(p.qlim)*[-1 1];

    AIV_trajectories(nm,'5',Nav,p5,force); 

    % b4=AIV_trajectories(nm,'4');
    % b4.linecolor  = 0.8*[1 1 1];
    % b4.markersize = 5;
    % b4.markerfacecolor = 'r'; %[0 0 1];
    % b4.markeredgecolor = 'k'; %[0 0 0];

    % Red particles - warm water

    b5=AIV_trajectories(nm,'5');
    p6=p;
    p6.rlim = p.rlim+0.0*diff(p.rlim)*[-1 1];
    p6.qlim = p.qlim+0.0*diff(p.qlim)*[-1 1];
    b6=AIV_traj_clip_polar(b5,p6.rlim,p6.qlim);

    b6.linecolor  = 0.8*[1 1 1];
    b6.markersize = 3;
    b6.markerfacecolor = 'r'; %[0 0 1];
    b6.markeredgecolor = 'r'; %[0 0 0];

    % Blue particles - fresh water
    b2=AIV_trajectories(nm,'2');

    b2.linecolor  = 0.8*[1 1 1];
    b2.markersize = 3;
    b2.markerfacecolor = 'b'; %[0 0 1];
    b2.markeredgecolor = 'b'; %[0 0 0];

    %CC=zeros(5,9);
    figure(1)
    subplot(2,1,1)
    disp('Blue')
    cb=CR_traj_fig1(b2,p6);
    title('Blue')
    subplot(2,1,2)
    disp('Red')
    cr=CR_traj_fig1(b6,p6);
    title('Red')
    results(idx, :) = {nm, cb.Pmax, cb.tmax, cb.Pmin, cb.tmin, cr.Pmax, cr.tmax, cr.Pmin, cr.tmin};

    %CC(3,:)=[0 cb.Pmax cb.tmax cb.Pmin cb.tmin cr.Pmax cr.tmax cr.Pmin cr.tmin]

    % % figure(4)
    % % subplot(2,1,1)
    % % plot(CC(:,1),CC(:,2)./CC(3,2),'b^')
    % % hold on
    % % plot(CC(:,1),CC(:,4)./CC(3,4),'bv')
    % % plot(CC(:,1),CC(:,8)./CC(3,8),'rv')
    % % xlabel('$\%$','interpreter','latex');
    % % ylabel('$Pmax Pmin$','interpreter','latex');
    % % 
    % % subplot(2,1,2)
    % % plot(CC(:,1),CC(:,3)./CC(3,3),'b^')
    % % hold on
    % % plot(CC(:,1),CC(:,5)./CC(3,5),'bv')
    % % plot(CC(:,1),CC(:,9)./CC(3,9),'rv')
    % % xlabel('$\%$','interpreter','latex');
    % % ylabel('$t_{max} t_{min}$\,(s)','interpreter','latex');

%     figure(2)
%     CR_traj_fig2(b2,p6)
% 
%     figure(3)
%     CR_traj_fig2(b6,p6)
% 
%     figure(4)
%     AIV_quiver(nm,Nav,p6,[b2 b6])
end
% Export results to CSV file
writetable(results, 'experiment_results.csv');

% Display the results
disp(results);
