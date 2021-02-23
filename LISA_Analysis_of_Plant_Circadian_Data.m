%% Title: Moore and Domijan, 2021, "Local Indicators of Spatial Autocorrelation (LISA) and Their Application to Circadian Spatio-temporal Data"
%% Author: Joshua Moore (University of Liverpool, Department of Mathematical Sciences)

close all
clc, clear
warning('off')
%% Set directory to where .mat files are stored
s = dir(pwd);
leaffile = 'CCA1LL_6.mat'


if isfile(fullfile(pwd,'breakxaxis.m')) && isfile(fullfile(pwd,'breakyaxis.m'))
    rng default
    nbins = 20;
    radius=1;
    skipPart = input('Do you want to show the example figures (Fig 1 & S1 in paper)? Y/N: ', 's');
    if (skipPart) == 'y' || (skipPart) == 'Y'
        Figure1=figure('NumberTitle', 'off', 'Name', 'Figure 1');
        set(Figure1,'Position',[0 0 1500 1500]) % Set figure size
        FigureS1=figure('NumberTitle', 'off', 'Name', 'Figure S1');
        set(gcf, 'Position', get(0, 'Screensize'));
        
        %% Plot example figure
        %% Case 2: Two clusters with behaviour that differs from mean.
        phases = rand(50,50).*pi; % Generate random matrix
        phases(10:30,10:20) = pi; % Cluster One: Phases = pi
        phases(40:45,30:45) = 0.*pi; % Cluster Two: Phases = 0
        phases = padarray(phases,[3 3],NaN,'both');
        
        
        %% Figure S1 Case 2
        [indx,indy] = find(~isnan(phases) == 1); % Identify indices of leaf cells
        Irand = zeros(1,1000); % Initialise randomised I
        for t = 1:1000 % 1000 Iterations
            idx = randperm(length(indx)); % Permute locations of cells
            
            indx2 = indx(idx); % Apply permutation to x coord of cells
            indy2 = indy(idx); % Apply permutation to x coord of cells
            random = phases;
            for j = 1:length(indx)
                random(indx(j),indy(j)) = phases(indx2(j),indy2(j)); % Conditional Randomisation
            end
            Irand(t) = nanmean(LocalMoran2(random,radius),'all'); % Compute local Moran for current randomisation
        end
        temp = reshape(squeeze(Irand),1,[]); % Squeeze Crand values into list
        temp(isnan(temp))=[]; % Remove any NaN
        pd = fitdist(temp','Kernel'); % Fit normal distribution to data
        xvals = -3:0.0001:3;  % Squeeze Crand values into list
        p = cdf(pd,xvals); % Compute cdf from fitted kernel distribution
        ll = (xvals(find(p >= 0.05,1))); % Alpha = 0.05
        ul = (xvals(find(p >= 0.95,1))); % Alpha = 0.05
        
        set(0, 'CurrentFigure', FigureS1)
        % Moran Case 2
        subplot(2,3,2)
        y = pdf(pd,xvals);
        d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of I values
        hold on
        plot(xvals,(y/max(y)).*max(d.Values),'b','LineWidth',4) % Plot pdf of I values
        hold on
        yvals = 0:1:650;
        plot(ones(size(yvals)).*nanmean(LocalMoran2(phases,radius),'all'),yvals,'r','LineWidth',4) % Plot actual Local Moran Value as red line
        plot(ones(size(yvals)).*ul,yvals,'--k','LineWidth',4) % Plot upper bound (alpha = 0.05)
        plot(ones(size(yvals)).*ll,yvals,'--k','LineWidth',4) % Plot lower bound (alpha = 0.05)
        xlim([-0.05 0.3]) % set xlim
        ylim([0 150]) % set ylim
        xticks([-0.05 0 0.05 0.25 0.3]) % set xticks
        set(gca,'fontsize',35) % set fontsize
        xlabel('$I^\theta$','fontsize',40,'interpreter','latex') % xlabel
        title('\textrm{\textbf{Case 2}}','fontsize',40,'interpreter','latex'); % Case number above each column
        breakxaxis([0.06,0.24]); % Insert x axis break
        
        
        [indx,indy] = find(~isnan(phases) == 1); % Identify indices of leaf cells
        Crand = zeros(1,1000); % Initialise randomised C
        for t = 1:1000 % 1000 Iterations
            idx = randperm(length(indx)); % Permute locations of cells
            indx2 = indx(idx); % Apply permutation to x coord of cells
            indy2 = indy(idx); % Apply permutation to x coord of cells
            random = phases;
            for j = 1:length(indx)
                random(indx(j),indy(j)) = phases(indx2(j),indy2(j)); % Conditional Randomisation
            end
            Crand(t) = nanmean(localGeary(random,radius),'all'); % Compute local Moran for current randomisation
        end
        
        
        temp = reshape(squeeze(Crand),1,[]); % Squeeze Crand values into list
        temp(isnan(temp))=[]; % Remove any NaN
        pd = fitdist(temp','Kernel'); % Fit normal distribution to data
        xvals = -3:0.0001:3;  % Squeeze Crand values into list
        p = cdf(pd,xvals); % Compute cdf from fitted kernel distribution
        ll = (xvals(find(p >= 0.05,1))); % Alpha = 0.05
        ul = (xvals(find(p >= 0.95,1))); % Alpha = 0.05
        
        subplot(2,3,5)
        % Geary Case 2
        y = pdf(pd,xvals);
        d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of C values
        hold on
        plot(xvals,(y/max(y)).*max(d.Values),'b','LineWidth',4) % Plot pdf of I values
        hold on
        yvals = 0:1:650;
        plot(ones(size(yvals)).*nanmean(localGeary(phases,radius),'all'),yvals,'r','LineWidth',4) % Plot actual Local Geary Value as red line
        plot(ones(size(yvals)).*ul,yvals,'--k','LineWidth',4) % Plot upper bound (alpha = 0.05)
        plot(ones(size(yvals)).*ll,yvals,'--k','LineWidth',4) % Plot lower bound (alpha = 0.05)
        xlim([0.7 1.05]) % Set xlim
        ylim([0 150]) % Set ylim
        xticks([0.7 0.75 0.95 1 1.05]) % Xticks
        set(gca,'fontsize',35) % Set fontsize
        xlabel('$C^\theta$','fontsize',40,'interpreter','latex') % Insert xlabel
        breakxaxis([0.76,0.9]); % Insert axis break
        
        
        
        
        %% Figure 1 Case 2
        set(0, 'CurrentFigure', Figure1)
        R = KuramotoOP((phases)); % Calculate r from phase data
        % SUBPLOT(3,3,2)
        ax(4) =subplot(3,3,2);
        imagesc(phases) % Plot image of phases
        xlabel(['$r$ = ',num2str(round(nanmean(R,'all'),3))],'fontsize',40,'interpreter','latex') % Report r value
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        hold on
        title('\textrm{\textbf{Case 2}}','fontsize',40,'interpreter','latex') % All upper row cells have case number title
        caxis([0 pi]) % Colorbar limits
        
        
        % SUBPLOT(3,3,5)
        ax(5) =subplot(3,3,5);
        [LI] = LocalMoran2((phases),radius); % Compute Local Moran
        imagesc(LI); % Plot image of local Moran
        caxis([0 1]) % Set colourbar limits
        xlabel(['$I^{\theta}$ = ',num2str(round(nanmean(LI,'all'),3))],'fontsize',40,'interpreter','latex') % Report global Moran value
        hold on
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        
        
        % SUBPLOT(3,3,6)
        ax(6) =subplot(3,3,8);
        % Compute Local Geary
        [LC] = localGeary((phases),radius);
        imagesc(LC); % Plot image of local Geary
        caxis([0 1]) % Set colourbar limits
        xlabel(['$C^{\theta}$ = ',num2str(round(nanmean(LC,'all'),3))],'fontsize',40,'interpreter','latex') % Report global Geary value
        hold on
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        
        
        
        %% Case 1: Randomly distributed phases.
        phases2 = rand(50,50).*pi;
        phases2 = padarray(phases2,[3 3],NaN,'both'); % Pad array with NaN
        
        [indx,indy] = find(~isnan(phases) == 1); % Identify indices of leaf cells
        Irand = zeros(1,1000);
        for t = 1:1000 % 1000 Iterations
            idx = randperm(length(indx)); % Permute locations of cells
            
            indx2 = indx(idx); % Apply permutation to x coord of cells
            indy2 = indy(idx); % Apply permutation to x coord of cells
            random = phases2;
            for j = 1:length(indx)
                random(indx(j),indy(j)) = phases2(indx2(j),indy2(j)); % Conditional Randomisation
            end
            Irand(t) = nanmean(LocalMoran2(random,radius),'all'); % Compute local Moran for current randomisation
        end
        
        temp = reshape(squeeze(Irand),1,[]); % Squeeze Crand values into list
        temp(isnan(temp))=[]; % Remove any NaN
        pd = fitdist(temp','Kernel'); % Fit normal distribution to data
        xvals = -3:0.0001:3;  % Squeeze Crand values into list
        p = cdf(pd,xvals); % Compute cdf from fitted kernel distribution
        ll = (xvals(find(p >= 0.05,1))); % Alpha = 0.05
        ul = (xvals(find(p >= 0.95,1))); % Alpha = 0.05
        
        
        %% Figure S1 Case 1
        set(0, 'CurrentFigure', FigureS1)
        % Moran Case 1
        subplot(2,3,1)
        y = pdf(pd,xvals);
        d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of I values
        hold on
        plot(xvals,(y/max(y)).*max(d.Values),'b','LineWidth',4) % Plot pdf of I values
        hold on
        yvals = 0:1:650;
        plot(ones(size(yvals)).*nanmean(LocalMoran2(phases,radius),'all'),yvals,'r','LineWidth',4) % Plot actual Local Moran Value as red line
        plot(ones(size(yvals)).*ul,yvals,'--k','LineWidth',4) % Plot upper bound (alpha = 0.05)
        plot(ones(size(yvals)).*ll,yvals,'--k','LineWidth',4) % Plot lower bound (alpha = 0.05)
        xlim([-0.05 0.05]) % Set xlim
        ylim([0 150]) % Set ylim
        xticks([-0.05 0 0.05]) % Set xticks
        set(gca,'fontsize',35) % Set font size for figure
        text(-0.1775,155,'(a)','fontsize',40,'interpreter','latex') % Insert (a) label first column only
        ylabel('Frequency','fontsize',40,'interpreter','latex') % ylabel for first column only
        xlabel('$I^\theta$','fontsize',40,'interpreter','latex') % Insert xlabel
        title('\textrm{\textbf{Case 1}}','fontsize',40,'interpreter','latex')% Case number above each column
        
        
        [indx,indy] = find(~isnan(phases2) == 1); % Identify indices of leaf cells
        radius=1;
        Crand = zeros(1,1000);
        for t = 1:1000 % 1000 Iterations
            idx = randperm(length(indx)); % Permute locations of cells
            
            indx2 = indx(idx); % Apply permutation to x coord of cells
            indy2 = indy(idx); % Apply permutation to x coord of cells
            random = phases2;
            for j = 1:length(indx)
                random(indx(j),indy(j)) = phases2(indx2(j),indy2(j)); % Conditional Randomisation
            end
            Crand(t) = nanmean(localGeary(random,radius),'all'); % Compute local Moran for current randomisation
        end
        
        
        temp = reshape(squeeze(Crand),1,[]); % Squeeze Crand values into list
        temp(isnan(temp))=[]; % Remove any NaN
        pd = fitdist(temp','Kernel'); % Fit normal distribution to data
        xvals = -3:0.0001:3;  % Squeeze Crand values into list
        p = cdf(pd,xvals); % Compute cdf from fitted kernel distribution
        ll = (xvals(find(p >= 0.05,1))); % Alpha = 0.05
        ul = (xvals(find(p >= 0.95,1))); % Alpha = 0.05
        
        subplot(2,3,4)
        % Geary Case 1
        y = pdf(pd,xvals);
        d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of C values
        hold on
        plot(xvals,(y/max(y)).*max(d.Values),'b','LineWidth',4) % Plot pdf of I values
        hold on
        yvals = 0:1:650;
        plot(ones(size(yvals)).*nanmean(localGeary(phases,radius),'all'),yvals,'r','LineWidth',4) % Plot actual Local Geary Value as red line
        plot(ones(size(yvals)).*ul,yvals,'--k','LineWidth',4) % Plot upper bound (alpha = 0.05)
        plot(ones(size(yvals)).*ll,yvals,'--k','LineWidth',4) % Plot lower bound (alpha = 0.05)
        xlim([0.935 1.02]) % Set xlim
        xticks([0.9 0.95 0.975 1 1.05]) % Assign xticks
        ylim([0 150]) % Set ylim
        set(gca,'fontsize',35) % Set font size
        text(0.55,155,'(b)','fontsize',40,'interpreter','latex') % Insert (b) label first column only
        ylabel('Frequency','fontsize',40,'interpreter','latex') % ylabel for first column only
        xlabel('$C^\theta$','fontsize',40,'interpreter','latex') % Insert xlabel
        
        
        
        %% Figure 1 Case 1
        set(0, 'CurrentFigure', Figure1)
        R = KuramotoOP((phases2)); % Calculate r from phase data
        % SUBPLOT(3,3,1)
        ax(1) = subplot(3,3,1);
        imagesc(phases2) % Plots phase as image
        caxis([0 pi]) % Set colourlimits
        xlabel(['$r$ = ',num2str(round(nanmean(R,'all'),3))],'fontsize',40,'interpreter','latex') % Report r value
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        alabel=ylabel('\textrm{\textbf{(a)}}','fontsize',40,'interpreter','latex','Rotation',0);% First column add (a) label
        alabel.Position(1) = -2; % Adjust alabel location
        alabel.Position(2) = 4; % Adjust alabel location
        hold on
        title('\textrm{\textbf{Case 1}}','fontsize',40,'interpreter','latex') % Add title with case number at top of column
        
        % SUBPLOT(3,3,2)
        ax(2) =subplot(3,3,4);
        LI = LocalMoran2((phases2),radius);  % Calculate Local Moran from phase data
        imagesc(LI) % Plot local Moran as image
        caxis([0 1]) % Set colourbar limits
        xlabel(['$I^{\theta}$ = ',num2str(round(nanmean(LI,'all'),3))],'fontsize',40,'interpreter','latex') % Report global Moran value
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        blabel=ylabel('\textrm{\textbf{(b)}}','fontsize',40,'interpreter','latex','Rotation',0); % First column add (b) label
        blabel.Position(1) = -2; % Adjust blabel location
        blabel.Position(2) = 4; % Adjust blabel location
        
        
        % SUBPLOT(3,3,7)
        ax(3) = subplot(3,3,7);
        % Compute local Geary
        [LC] = localGeary((phases2),radius); % Compute local Geary
        imagesc(LC); % Plot image of local Geary
        caxis([0 1]) % Set colourbar limits
        xlabel(['$C^{\theta}$ = ',num2str(round(nanmean(LC,'all'),3))],'fontsize',40,'interpreter','latex') % Report global Geary value
        hold on
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        clabel=ylabel('\textrm{\textbf{(c)}}','fontsize',40,'interpreter','latex','Rotation',0); % Add clabel
        clabel.Position(1) = -2; % Adjust location of clabel
        clabel.Position(2) = 4;  % Adjust location of clabel
        
        
        
        %% Case 3: Clusters with behaviour roughly equal to the mean.
        phases2 = rand(50,50).*pi;
        phases2(10:30,10:20) = 0.5.*pi; % Cluster One: Phases = pi/2
        phases2(40:45,30:45) = 0.5.*pi; % Cluster Two: Phases = pi/2
        phases2 = padarray(phases2,[3 3],NaN,'both');
        
        %% Figure S1 Case 3
        [indx,indy] = find(~isnan(phases) == 1); % Identify indices of leaf cells
        Irand = zeros(1,1000); % Initialise Randomisation I
        for t = 1:1000 % 1000 Iterations
            idx = randperm(length(indx)); % Permute locations of cells
            
            indx2 = indx(idx); % Apply permutation to x coord of cells
            indy2 = indy(idx); % Apply permutation to x coord of cells
            random = phases2;
            for j = 1:length(indx)
                random(indx(j),indy(j)) = phases2(indx2(j),indy2(j)); % Conditional Randomisation
            end
            Irand(t) = nanmean(LocalMoran2(random,radius),'all'); % Compute local Moran for current randomisation
        end
        
        
        temp = reshape(squeeze(Irand),1,[]); % Squeeze Crand values into list
        temp(isnan(temp))=[]; % Remove any NaN
        pd = fitdist(temp','Kernel'); % Fit normal distribution to data
        xvals = -3:0.0001:3;  % Squeeze Crand values into list
        p = cdf(pd,xvals); % Compute cdf from fitted kernel distribution
        ll = (xvals(find(p >= 0.05,1))); % Alpha = 0.05
        ul = (xvals(find(p >= 0.95,1))); % Alpha = 0.05
        
        
        
        set(0, 'CurrentFigure', FigureS1)
        % Moran Case 3
        subplot(2,3,3)
        y = pdf(pd,xvals);
        d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of I values
        hold on
        plot(xvals,(y/max(y)).*max(d.Values),'b','LineWidth',4) % Plot pdf of I values
        hold on
        yvals = 0:1:650;
        plot(ones(size(yvals)).*nanmean(LocalMoran2(phases,radius),'all'),yvals,'r','LineWidth',4) % Plot actual Local Moran Value as red line
        plot(ones(size(yvals)).*ul,yvals,'--k','LineWidth',4) % Plot upper bound (alpha = 0.05)
        plot(ones(size(yvals)).*ll,yvals,'--k','LineWidth',4) % Plot lower bound (alpha = 0.05)
        ylim([0 150])
        xlim([-0.05 0.05])
        xticks([-0.05 0 0.05])
        set(gca,'fontsize',35)
        xlabel('$I^\theta$','fontsize',40,'interpreter','latex')
        title('\textrm{\textbf{Case 3}}','fontsize',40,'interpreter','latex')
        
        [indx,indy] = find(~isnan(phases2) == 1); % Identify indices of leaf cells
        Crand = zeros(1,1000); % Initialise Randomised C
        for t = 1:1000 % 1000 Iterations
            idx = randperm(length(indx)); % Permute locations of cells
            indx2 = indx(idx); % Apply permutation to x coord of cells
            indy2 = indy(idx); % Apply permutation to x coord of cells
            random = phases2;
            for j = 1:length(indx)
                random(indx(j),indy(j)) = phases2(indx2(j),indy2(j)); % Conditional Randomisation
            end
            Crand(t) = nanmean(localGeary(random,radius),'all'); % Compute local Moran for current randomisation
        end
        
        
        temp = reshape(squeeze(Crand),1,[]); % Squeeze Crand values into list
        temp(isnan(temp))=[]; % Remove any NaN
        pd = fitdist(temp','Kernel'); % Fit normal distribution to data
        xvals = -3:0.0001:3;  % Squeeze Crand values into list
        p = cdf(pd,xvals); % Compute cdf from fitted kernel distribution
        ll = (xvals(find(p >= 0.05,1))); % Alpha = 0.05
        ul = (xvals(find(p >= 0.95,1))); % Alpha = 0.05
        
        subplot(2,3,6)
        % Geary Case 3
        y = pdf(pd,xvals);
        d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of C values
        hold on
        plot(xvals,(y/max(y)).*max(d.Values),'b','LineWidth',4) % Plot pdf of I values
        hold on
        yvals = 0:1:650; % yvals for ul,ll and local Geary
        plot(ones(size(yvals)).*nanmean(localGeary(phases,radius),'all'),yvals,'r','LineWidth',4) % Plot actual Local Geary Value as red line
        plot(ones(size(yvals)).*ul,yvals,'--k','LineWidth',4) % Plot upper bound (alpha = 0.05)
        plot(ones(size(yvals)).*ll,yvals,'--k','LineWidth',4) % Plot lower bound (alpha = 0.05)
        xlim([0.935 1.02]) % Set xlim
        xticks([0.9 0.95 0.975 1 1.05]) % Set xticks
        ylim([0 150]) % Set ylim
        set(gca,'fontsize',35) % Set fontsize
        xlabel('$C^\theta$','fontsize',40,'interpreter','latex') % xlabel
        
        
        
        
        
        
        %% Figure 1 Case 3
        set(0, 'CurrentFigure', Figure1)
        R = KuramotoOP((phases2));
        % SUBPLOT(3,3,3)
        ax(7) =subplot(3,3,3);
        imagesc(phases2)
        h=colorbar('XTickLabel',{'0','\pi/2','\pi'},'XTick', 0:pi/2:pi,'FontSize',25); % add colourbar with specified xticks
        ylabel(h, '$\theta \; (rad)$','fontsize',40,'interpreter','latex') % Add label to colourbar
        caxis([0 pi]) % Set colour limits
        xlabel(['$r$ = ',num2str(round(nanmean(R,'all'),3))],'fontsize',40,'interpreter','latex') % Report r value
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        hold on
        title('\textrm{\textbf{Case 3}}','fontsize',40,'interpreter','latex') % Add case number at top of column
        
        
        % SUBPLOT(3,3,6)
        ax(8) =subplot(3,3,6);
        % Compute local Moran
        [LI] = LocalMoran2((phases2),radius);
        imagesc(LI);% Plot image of local Moran
        h=colorbar; % Colourbar
        h.FontSize = 25; % Set font size of colourbar
        ylabel(h, '$I^{\theta}_{i}$','fontsize',40,'interpreter','latex') % Assign label to colourbar
        caxis([0 1]) % Set colourbar limits
        xlabel(['$I^{\theta}$ = ',num2str(round(nanmean(LI,'all'),3))],'fontsize',40,'interpreter','latex') % Report global Moran value
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        
        % SUBPLOT(3,3,9)
        ax(9) =subplot(3,3,9);
        % Compute local Geary
        [LC] = localGeary((phases2),radius);
        imagesc(LC); % Plot image of local Geary
        h=colorbar; % Colourbar
        h.FontSize = 25; % Set font size of colourbar
        ylabel(h, '$C^{\theta}_{i}$','fontsize',40,'interpreter','latex') % Assign label to colourbar
        caxis([0 1]) % Set colourbar limits
        xlabel(['$C^{\theta}$ = ',num2str(round(nanmean(LC,'all'),3))],'fontsize',40,'interpreter','latex') % Report global Geary value
        set(gca,'xtick',[]) % Remove xtick
        set(gca,'ytick',[]) % Remove ytick
        xlim([4 53]) % Set xlim
        ylim([4 53]) % Set ylim
        
        
        
        
        % Adjust size of subplots.
        set(ax(1), 'Position', [0.1, 0.7, 0.24, 0.24]);
        set(ax(2), 'Position', [0.1, 0.4, 0.24, 0.24]);
        set(ax(3), 'Position', [0.1, 0.1, 0.24, 0.24]);
        set(ax(4), 'Position', [0.375, 0.7, 0.24, 0.24]);
        set(ax(5), 'Position', [0.375, 0.4, 0.24, 0.24]);
        set(ax(6), 'Position', [0.375, 0.1, 0.24, 0.24]);
        set(ax(7), 'Position', [0.65, 0.7, 0.24, 0.24]);
        set(ax(8), 'Position', [0.65, 0.4, 0.24, 0.24]);
        set(ax(9), 'Position', [0.65, 0.1, 0.24, 0.24]);
        
        
        % Set Correct colourmaps
        colormap(ax(1),parula);
        colormap(ax(2),parula);
        colormap(ax(4),parula);
        colormap(ax(5),parula);
        colormap(ax(7),parula);
        colormap(ax(8),parula);
        % Geary's C scale is reversed so that high spatial autocorrelation is still
        % yellow
        colormap(ax(3),flipud(parula));
        colormap(ax(6),flipud(parula));
        colormap(ax(9),flipud(parula));
        
        
        
        
    end
    
    
    leafn = 1; % Initialise leaf number counter counter for subplots

    if isfile(leaffile)
        clear x1 x2 y1 y2 % Clear previous xlim and ylim
        data = load(leaffile);
        
        rawexpression = data.timecourse; % raw expression is given by lucpre_comat
        % Uncomment below if you wish to see selected leaf pixels
        %                 figure
        %                 imagesc(nanmean(rawexpression,3).*(nanmean(rawexpression,3)>nanmean(rawexpression,'all')*1.5))
        %
        leaf = (nanmean(rawexpression,3)>nanmean(rawexpression,'all')*1.5); % Remove pixels that are not part of leaf. Cells outside of leaf are well identified as those with expression < 1.5*mean.
        
        
        
        %% Calculate phases using Pikovsky method
        phases = zeros(size(rawexpression)); % Initialise phase array
        for m = 1:length(rawexpression(:,1,1))
            for j = 1:length(rawexpression(1,:,1))
                % For m,j = all spatial indices
                if leaf(m,j) == 1 % If cell m,j is part of leaf
                    [peaks,indices] = findpeaks(lowpass(squeeze(rawexpression(m,j,:)),0.01)); % Identify peaks in raw expression
                    for t = 1:length(rawexpression(1,1,:))
                        index1 = find(indices > t,1); % Find the first peak occuring before time t
                        index2 = find(indices <= t,1,'last'); % Find the first peak occuring after time t
                        if isempty(index2)
                            indices = [1; indices]; % If t is before first peak work introduce peak at t=1 (removed later).
                            index2 = 1; % Index of peak before t is therefore 1.
                        end
                        if isempty(index1)
                            indices = [indices; length(rawexpression(1,1,:))]; % If t is later than final peak, introduce peak at final time point.
                            index1 = length(indices); % Index of final peak is therefore the length of indices.
                        end
                        % Use linear interpolation to identify phases in between peaks
                        phases(m,j,t) = 2*pi*index2 + ((t - indices(index2))/(indices(index1) - indices(index2)))*2*pi;
                    end
                else
                    % Cells outside of leaf: phases = NaN;
                    phases(m,j,:) = NaN;
                end
            end
        end
        
        
        %% Remove first 24 hours of phases due to poor performance of Pikovsky. For consistency consider only 72 hours of data.
        phases= phases(:,:,(24/(2/3)):min((24/(2/3))+72/(2/3),length(phases(1,1,:)))); % If length of data is less than 72 hours, consider as much as possible.
        
        
        
        
        
        radius = 1; % Radius of Influence: This is used in local Moran, local Geary and Moran scatter plot calculations.
        % This is synonymous to the radius at which one would
        % expect cells to interact.
        
        %% Synchrony Measures Calculation
        % Initialise Measures
        LocalMI = zeros(size(phases)); % Local Moran
        devphase = zeros(size(phases(:,:,1)));  % Deviation of phase from average
        LocalGC = zeros(size(phases)); % Local Geary C
        GearyC = zeros(1,length(phases(1,1,:))); % Global Geary C
        KuramotoOrderParameter = zeros(1,length(phases(1,1,:))); % Kuramoto Order Parameter
        MoransI = zeros(1,length(phases(1,1,:))); % Global Morans I
        for t = 1:length(phases(1,1,:))
            [LocalMI(:,:,t),devphase(:,:,t)] = LocalMoran2((phases(:,:,t)),radius); %Compute Local Moran
            [LocalGC(:,:,t),GearyC(t)] = localGeary((phases(:,:,t)),radius); % Compute Local Geary
            KuramotoOrderParameter(t) = KuramotoOP(phases(:,:,t)); % Compute Kuramoto Order Parameter
            MoransI(t) = nanmean(LocalMI(:,:,t),'all'); % Global Moran's I is the mean of local components
        end
        
        
        
        TimeAveragedLocalMoran = nanmean(LocalMI,3);
        %% Find good xlim and ylim for imagesc plots
        n = 0; % Initialise counter
        for m = 1:length(TimeAveragedLocalMoran(1,:))
            
            if any(~isnan(TimeAveragedLocalMoran(:,m))) && n == 0
                % Find first column with leaf cell
                x1 = m; % Column with first leaf cell
                n = 1; % Counter: Indicates whether looking for first or final leaf column. n = 1: found first column, looking for final
            elseif  any(~isnan(TimeAveragedLocalMoran(:,m))) == 0 && n == 1
                % Find final column with leaf cell
                x2 = m; % Column with final leaf cell
                n = 2; % n = 2: Found both final and first column
            end
        end
        if ~exist('x2','var')
            % If x2 can't be identified, use length of data in x direction
            x2 = length(TimeAveragedLocalMoran(1,:));
        end
        if ~exist('x1','var')
            % If x1 can't be identified, use first column
            x1 = 1;
        end
        
        n=0; % Initialise counter
        for m = 1:length(TimeAveragedLocalMoran(:,1))
            if any(~isnan(TimeAveragedLocalMoran(m,:))) && n == 0
                % Find first row with leaf cell
                y1 = m; % Row with first leaf cell
                n = 1; % Counter: Indicates whether looking for first or final leaf row. n = 1: found first row, looking for final
            elseif any(~isnan(TimeAveragedLocalMoran(m,:))) == 0 && n == 1
                % Find final row with leaf cell
                y2 = m;  % Row with final leaf cell
                n = 2; % n = 2: Found both final and first Row
            end
        end
        if ~exist('y2','var')
            % If y2 can't be identified, use length of data in y direction
            y2 = length(TimeAveragedLocalMoran(:,1));
        end
        if ~exist('y1','var')
            % If y1 can't be identified, use first row
            y1 = 1;
        end
        
        
        
        
        
        
        
        
        %% Idenitfy Significant Local Moran Clusters: Uncomment for significance calculations for leaf 20
        Significance = zeros(size(phases)); % Initialise Significance
        commandwindow
        skipPart = input('Do you want to conduct significance calculations? Y/N', 's');
        if (skipPart) == 'y' || (skipPart) == 'Y'
            
            figure('NumberTitle', 'off', 'Name', 'Figure S5');
            [indx,indy] = find(leaf == 1); % Identify indices of leaf cells
            for m = 1:length(LocalMI(1,1,:))
                LISArand = zeros(length(phases(:,1,1)),length(phases(1,:,1)),1000);
                for t = 1:1000 % 1000 Iterations
                    idx = randperm(length(indx)); % Permute locations of cells
                    
                    indx2 = indx(idx); % Apply permutation to x coord of cells
                    indy2 = indy(idx); % Apply permutation to x coord of cells
                    random = phases(:,:,m);
                    for j = 1:length(indx)
                        random(indx(j),indy(j)) = phases(indx2(j),indy2(j),m); % Conditional Randomisation
                    end
                    LISArand(:,:,t) = LocalMoran2(random,radius); % Compute local Moran for current randomisation
                end
                
                % Initialise upper and lower limit for significance calculations
                ul = zeros(size(leaf));
                ll = ul;
                set(gcf,'Position',[0 0 2500 650])
                
                
                for t = 1:length(LISArand(:,1,1))
                    for k = 1:length(LISArand(1,:,1))
                        if leaf(t,k) == 1
                            temp = squeeze(LISArand(t,k,:));
                            temp = reshape(temp,1,[]);
                            temp(isnan(temp))=[];
                            pd = fitdist(temp','Kernel'); % Fit normal distribution to data
                            xvals = -6+0.05:0.05:6;
                            p = cdf(pd,xvals); % Compute cdf from fitted normal distribution
                            %                     ll(t,k) = (xvals(find(p >= 0.05/nansum(leaf,'all'),1))); % Bonferroni Condition
                            %                     ul(t,k) = (xvals(find(p >= 1-(0.05/nansum(leaf,'all')),1)));  % Bonferroni Condition
                            ll(t,k) = (xvals(find(p >= 0.05,1))); % Alpha = 0.05
                            ul(t,k) = (xvals(find(p >= 0.95,1))); % Alpha = 0.05
                            
                            if (k == 18 && t == 12 && m == 1 && sum(leaffile == 'CCA1LL_6.mat')==length(leaffile))
                                
                                subplot(2,12,[4 7])
                                d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of local Moran values
                                hold on
                                y = pdf(pd,xvals);  % y values of pdf at xvals
                                plot(xvals,(y/max(y))*max(d.Values),'b','LineWidth',5);  % Plot pdf of local Moran values
                                xlim([-2.5 2.5]) % Set xlim
                                ylim([0 550])
                                set(gca,'fontsize',16) % Increase size of xticklabels
                                ylabel('\textrm{Frequency}','interpreter','latex','FontSize',25) % Insert y label
                                xlabel('\textrm{Moran''s \(I\)}','interpreter','latex','FontSize',25) % Insert x label
                                ta=title('\textrm{\textbf{(b)}}','interpreter','latex','FontSize',30); % Add b label
                                ta.Position(1) = -3.5; % Adjust x coordinate of b label
                                ta.Position(2) = 550; % Adjust y coordinate of b label
                                yvals = 0:1:650; % Set y values
                                plot(ones(size(yvals)).*ll(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*ul(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*(xvals(find(p >= 0.05/nansum(leaf,'all'),1))),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*xvals(find(p >= 1-(0.05/nansum(leaf,'all')),1)),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*(xvals(find(p >= (1-(0.05)^(1/nansum(leaf,'all'))),1))),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*xvals(find(p >= (1-(1-(0.05)^(1/nansum(leaf,'all')))),1)),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*LocalMI(t,k,m),yvals,'r','LineWidth',4) % Plot actual local Moran value
                                text(-2.4,450,'1','interpreter','latex','FontSize',25) % Insert cell label
                                subplot(2,12,[1 2 13 14])
                                % Add coloured pixel to first subplot
                                im = zeros(length(leaf(:,1)),length(leaf(1,:)),3);
                                im(t,k,1) = 1; % Red pixel at location t,k
                                p=imagesc(im); % plot pixel
                                set(p,'AlphaData',nansum(im,3)>0) % Set all not t,k to invisible
                                hold on
                                text(k+1,t+1,'1','interpreter','latex','FontSize',25) % Label pixel
                            elseif (k == 10 && t == 20 && m == 1 && sum(leaffile == 'CCA1LL_6.mat')
                                subplot(2,12,[9 12])
                                d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of local Moran values
                                hold on
                                y = pdf(pd,xvals);  % y values of pdf at xvals
                                plot(xvals,(y/max(y))*max(d.Values),'b','LineWidth',5);  % Plot pdf of local Moran values
                                
                                xlim([-2.5 2.5]) % Set xlim
                                ylim([0 550])
                                set(gca,'fontsize',16) % Increase size of xticklabels
                                ylabel('\textrm{Frequency}','interpreter','latex','FontSize',25) % Insert y label
                                xlabel('\textrm{Moran''s \(I\)}','interpreter','latex','FontSize',25) % Insert x label
                                yvals = 0:1:650; % Set y values
                                plot(ones(size(yvals)).*ll(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*ul(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*(xvals(find(p >= 0.05/nansum(leaf,'all'),1))),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*xvals(find(p >= 1-(0.05/nansum(leaf,'all')),1)),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*(xvals(find(p >= (1-(0.05)^(1/nansum(leaf,'all'))),1))),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*xvals(find(p >= (1-(1-(0.05)^(1/nansum(leaf,'all')))),1)),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*LocalMI(t,k,m),yvals,'r','LineWidth',4) % Plot actual local Moran value
                                text(-2.4,450,'2','interpreter','latex','FontSize',25)
                                subplot(2,12,[1 2 13 14])
                                % Add coloured pixel to first subplot
                                im = zeros(length(leaf(:,1)),length(leaf(1,:)),3);
                                im(t,k,1) = 1; % Red pixel at location t,k
                                p=imagesc(im); % plot pixel
                                set(p,'AlphaData',nansum(im,3)>0) % Set all not t,k to invisible
                                hold on
                                text(k+1,t+1,'2','interpreter','latex','FontSize',25) % Label pixel
                            elseif (k == 25 && t == 28 && m == 1 && sum(leaffile == 'CCA1LL_6.mat')
                                subplot(2,12,[16 19])
                                d=histogram(temp',nbins,'FaceColor','w','LineWidth',2);
                                hold on
                                y = pdf(pd,xvals);
                                plot(xvals,(y/max(y))*max(d.Values),'b','LineWidth',5);
                                xlim([-2.5 2.5]) % Set xlim
                                ylim([0 550])
                                set(gca,'fontsize',16) % Increase size of xticklabels
                                ylabel('\textrm{Frequency}','interpreter','latex','FontSize',25) % Insert y label
                                xlabel('\textrm{Moran''s \(I\)}','interpreter','latex','FontSize',25) % Insert x label
                                
                                yvals = 0:1:650; % Set y values
                                plot(ones(size(yvals)).*ll(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*ul(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*(xvals(find(p >= 0.05/nansum(leaf,'all'),1))),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*xvals(find(p >= 1-(0.05/nansum(leaf,'all')),1)),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*(xvals(find(p >= (1-(0.05)^(1/nansum(leaf,'all'))),1))),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*xvals(find(p >= (1-(1-(0.05)^(1/nansum(leaf,'all')))),1)),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*LocalMI(t,k,m),yvals,'r','LineWidth',4) % Plot actual local Moran value
                                text(-2.4,450,'3','interpreter','latex','FontSize',25)
                                subplot(2,12,[1 2 13 14])
                                % Add coloured pixel to first subplot
                                im = zeros(length(leaf(:,1)),length(leaf(1,:)),3);
                                im(t,k,1) = 1; % Red pixel at location t,k
                                p=imagesc(im); % plot pixel
                                set(p,'AlphaData',nansum(im,3)>0) % Set all not t,k to invisible
                                hold on
                                text(k+1,t+1,'3','interpreter','latex','FontSize',25) % Label pixel
                            elseif (k == 19 && t == 32 && m == 1 && sum(leaffile == 'CCA1LL_6.mat')
                                subplot(2,12,[21 24])
                                d=histogram(temp',nbins,'FaceColor','w','LineWidth',2); % Plot histogram of local Moran values
                                hold on
                                y = pdf(pd,xvals);  % y values of pdf at xvals
                                plot(xvals,(y/max(y))*max(d.Values),'b','LineWidth',5);  % Plot pdf of local Moran values
                                xlim([-2.5 2.5]) % Set xlim
                                ylim([0 550]) % Set ylim
                                set(gca,'fontsize',16) % Increase size of xticklabels
                                ylabel('\textrm{Frequency}','interpreter','latex','FontSize',25) % Insert y label
                                xlabel('\textrm{Moran''s \(I\)}','interpreter','latex','FontSize',25) % Insert x label
                                yvals = 0:1:650; % Set y values
                                plot(ones(size(yvals)).*ll(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*ul(t,k),yvals,'--k','LineWidth',4) % Alpha = 0.05
                                plot(ones(size(yvals)).*(xvals(find(p >= 0.05/nansum(leaf,'all'),1))),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*xvals(find(p >= 1-(0.05/nansum(leaf,'all')),1)),yvals,'-.m','LineWidth',4) % Bonferroni Condition
                                plot(ones(size(yvals)).*(xvals(find(p >= (1-(0.05)^(1/nansum(leaf,'all'))),1))),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*xvals(find(p >= (1-(1-(0.05)^(1/nansum(leaf,'all')))),1)),yvals,':','LineWidth',4,'Color',[190 190 190]/255) % Sidak Condition
                                plot(ones(size(yvals)).*LocalMI(t,k,m),yvals,'r','LineWidth',4) % Plot actual local Moran value
                                text(-2.4,450,'4','interpreter','latex','FontSize',25)
                                subplot(2,12,[1 2 13 14])
                                % Add coloured pixel to first subplot
                                im = zeros(length(leaf(:,1)),length(leaf(1,:)),3);
                                im(t,k,1) = 1; % Red pixel at location t,k
                                p=imagesc(im); % plot pixel
                                set(p,'AlphaData',nansum(im,3)>0) % Set all not t,k to invisible
                                hold on
                                text(k+1,t+1,'4','interpreter','latex','FontSize',25) % Place label South East of cell
                                subplot(2,12,[1 2 13 14])
                                Clusters = leaf;
                                binaryImage = Clusters; % binary image for cluster location
                                boundaries = bwboundaries(binaryImage); % Identify cluster boundary
                                x = boundaries{1}(:, 2); % x coord of boundary
                                y = boundaries{1}(:, 1); % y coord of boundary
                                hold on
                                plot(x,y,'color','k','LineWidth',4) % plot leaf outline
                                xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
                                ylim([y1-1 y2+1])
                                set(gca,'xtick',[]) % Remove x ticks
                                set(gca,'ytick',[]) % Remove y ticks
                                hold on
                                tb=title('\textrm{\textbf{(a)}}','interpreter','latex','FontSize',30); % Add b label
                                tb.Position(1) = x1-7; % Adjust x coordinate of b label
                                tb.Position(2) = y1-1; % Adjust y coordinate of b label
                                
                            end
                        end
                    end
                end
                
                temp = Significance(:,:,m); % Initialise 2D matrix
                temp(LocalMI(:,:,m)>=ul) = 1; % Cells in temp with local Moran greater than upper limit are significant
                Significance(:,:,m) = temp; % Assign temp to significance matrix
            end
            
        end
        
        
        
        
        
        
        
        
        %% Moran scatter plot computation
        K = ones(2*radius + 1,2*radius + 1); % Determine influence matrix from radius. This is the kernel used in convolutions.
        K(radius + 1,radius + 1) = 0; % Cell do not interact with themselves -> central cell of K = 0
        
        % Initialise Moran scatter plot data
        MoranScatterData = zeros(size(devphase));
        MoranScatterX = zeros(length(devphase(1,1,:)),length(reshape(TimeAveragedLocalMoran,1,[])));
        MoranScatterY = zeros(length(devphase(1,1,:)),length(reshape(TimeAveragedLocalMoran,1,[])));
        
        for t = 1:length(devphase(1,1,:))
            temp = (devphase(:,:,t)); % Assign variable of interest
            temp(isnan(TimeAveragedLocalMoran)) = 0; % Remove NaN for convolution
            spatiallag = conv2(temp,K,'same')./(conv2(leaf,K,'same')); % Determine spatial lagged component of devphase = (sum of neighbours)/(number of neighbours). This is the convolutional equivalent of row standardisation.
            temp(~leaf) = NaN; % temp that are not part of leaf = NaN
            spatiallag(~leaf) = NaN; % spatiallag that are not part of leaf = NaN
            
            % Reshape into one-dimensional vectors
            MoranScatterX(t,:) = reshape(temp,1,[]);
            MoranScatterY(t,:) = reshape(spatiallag,1,[]);
            
            % Center on the origin (subtract mean). For Moran scatte plot,
            % plot these vectors against eachother e.g.
            % plot(MoranScatterX(t,:),MoranScatterY(t,:))
            MoranScatterX(t,:) = MoranScatterX(t,:);%-nanmean(MoranScatterX(t,:));
            MoranScatterY(t,:) = MoranScatterY(t,:);%-nanmean(MoranScatterY(t,:));
            
            % Reshape into shape of leaf
            MoranScatterMatY = reshape(MoranScatterY(t,:),size(leaf));
            MoranScatterMatX = reshape(MoranScatterX(t,:),size(leaf));
            
            % Identify different clsuters
            temp = MoranScatterData(:,:,t);
            temp(MoranScatterMatY > 0 & MoranScatterMatX < 0) = 1; % Low-High
            temp(MoranScatterMatX > 0 & MoranScatterMatY < 0) = 2; % High-Low
            temp((MoranScatterMatX > 0 & MoranScatterMatY > 0)) = 3; % High-High
            temp((MoranScatterMatX < 0 & MoranScatterMatY < 0))= 4; %Low-Low
            MoranScatterData(:,:,t) = temp;
        end
        MoranScatterClusters = zeros(size(MoranScatterData(:,:,1)));
        MoranScatterMean = nanmean(MoranScatterData(:,:,end),3); % behaviour final time
        MoranScatterStandardDev = nanstd(MoranScatterData(:,:,end),[],3); % Standard deviation of behaviour at final time
        MoranScatterClusters(MoranScatterMean == 3 & MoranScatterStandardDev == 0) = 0; % If high-high at final time
        MoranScatterClusters(MoranScatterMean == 4 & MoranScatterStandardDev == 0) = 1; % If low-low at final time
        MoranScatterClusters(MoranScatterMean ~= 4 & MoranScatterMean ~= 3) = 2; % If anything else at final time
        MoranScatterClusters(~leaf) = NaN; % Remove zeros due to NaN values.
        
        % Identify cells that are significant at final time
        Zt = zeros(size(Significance(:,:,1)));
        SignficanceAtT = nanmean(Significance(:,:,end),3); % If cells significant for all time are desired: remove "(:,:,end)" and average through time e.g. nanmean(Significance,3).
        Zt(abs(SignficanceAtT) == 1) = 1; % One if significant
        Zt(abs(SignficanceAtT) == 0) = 0; % Else Zt = 0
        Zt(~leaf) = NaN; %Remove zeros due to NaN values
        
        
        
        %% Moran scatter plot figure: Final Time Point
        [x,y] = prepareCurveData((MoranScatterX(end,:))',(MoranScatterY(end,:))'); % Prepare curve data for linear fit (Mainly just removing NaN & non unique x, y).
        f=fit(x,y,'poly1'); % Linear fit
        figure('NumberTitle', 'off', 'Name', 'Figure 2','Position',[0 0 1000 500]);
        subplot(1,2,1)
        plot((-5:5),zeros(length(-5:5)),'--k') % Mark y-axis
        hold on
        plot(zeros(length(-5:5)),(-5:5),'--k') % Mark x-axis
        scatter((MoranScatterX(end,:)),(MoranScatterY(end,:)),'r','filled') % Overlay scatter points so error bars are behind
        xlim([-1.5 2]) % Set consistent limits for all leaves
        ylim([-1.5 2]) % Set consistent limits for all leaves
        plot((-5:5),nanmean(f.p1)*(-5:5)+f.p2,'k','LineWidth',1.5) % Plot time averaged Moran's I
        text(1.1, 1.75,'\textrm{Q1}','interpreter','latex','FontSize',25) % Quadrant 1 label
        text(1.1,-1.25,'\textrm{Q4}','interpreter','latex','FontSize',25) % Quadrant 2 label
        text(-0.6, -1.25,'\textrm{Q3}','interpreter','latex','FontSize',25) % Quadrant 3 label
        text(-0.6, 1.75,'\textrm{Q2}','interpreter','latex','FontSize',25) % Quadrant 4 label
        
        xticks([-1 0 1 2]) % Set x tick values
        yticks([-1 0 1 2]) % Set y tick values
        d = get(gca,'XTickLabel'); % Find the x tick label
        set(gca,'XTickLabel',d,'fontsize',16) % Increase size of xticklabels
        ta=title('\textrm{\textbf{(a)}}','interpreter','latex','FontSize',30); % Create (a) subplot label
        ta.Position(1) = -2; % Adjust x coordinate of (a)
        ta.Position(2) = 1.9; % Adjust y coordinate of (a)
        xlabel('$d(\theta_i, \bar{\theta})$','interpreter','latex','FontSize',30) % Label x-axis
        ylabel('$\sum_j w_{ij}d(\theta_j, \bar{\theta})$','interpreter','latex','FontSize',30) % Label y-axis
        
        
        d = sign(nanmean(reshape(MoranScatterX(end,:),size(leaf)),3)).*...
            (abs(1-(nanmean(reshape(MoranScatterX(end,:),size(leaf)),3)./nanmean(reshape(MoranScatterY(end,:),size(leaf)),3)))); % Plot how close gradient is to one
        subplot(1,2,2)
        ax1=subplot(1,2,2); % Create axis 1 for one plot
        ax2=subplot(1,2,2); % Create axis 2 for overlay
        set(ax1,'Visible','off'); % Set ax1 visibility
        set(ax2,'Visible','off'); % Set ax1 visibility
        n =80; % Set number of increments in colorbar
        cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]'; % Set blue half of colourbar
        cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]'; % Set red half of colourbar
        cmap = [cmap1; cmap2(2:end, :)]; % Combine for complete colourmap
        
        p1=imagesc(ax2,d); % Plot d variable as image on ax2
        colormap(cmap) % Use custom colourmap
        caxis([-0.1 0.1]) % Colormap has limits [-0.1 0.1]
        h = colorbar('XTick', -0.1:0.05:0.1,'FontSize',25); %Set colorbar increments and fontsize
        slabel=ylabel(h, '$S_{i}$','fontsize',30,'interpreter','latex'); % Label colourbar
        slabel.Position(1) = 2.75; % Shift colourbar label slightly to left
        hold on
        set(p1,'AlphaData',~isnan(MoranScatterClusters)) % Set NaN to white.
        im= []; % initialise image
        for i = 1:3
            temp = MoranScatterClusters==2;
            temp(~(MoranScatterClusters==2)) = 1; % Set to one if cells is in Q1 or Q3
            temp = temp-1;
            im(:,:,i) = temp; % Im is rgb array which is black for all cells in Q2 or Q4
        end
        p2 = imagesc(im); % Plot im
        set(p2,'AlphaData',MoranScatterClusters==2)
        set(gca,'xtick',[]) % Remove x ticks
        set(gca,'ytick',[]) % Remove y ticks
        xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
        ylim([y1-1 y2+1])
        tb=title('\textrm{\textbf{(b)}}','interpreter','latex','FontSize',30); % Add b label
        tb.Position(1) = 3.5; % Adjust x coordinate of b label
        tb.Position(2) = 20; % Adjust y coordinate of b label
        
        
        
        
        
        
        %% Compute condensed Moran scatter plot
        MoranScatterClusters = zeros(size(MoranScatterData(:,:,1)));
        MoranScatterMean = nanmean(MoranScatterData(:,:,:),3); % Average behaviour through time
        MoranScatterStandardDev = nanstd(MoranScatterData(:,:,:),[],3); % Standard deviation of behaviour through time
        MoranScatterClusters(MoranScatterMean == 3 & MoranScatterStandardDev == 0) = 0; % If high-high for all time
        MoranScatterClusters(MoranScatterMean == 4 & MoranScatterStandardDev == 0) = 1; % If low-low for all time
        MoranScatterClusters(MoranScatterMean ~= 4 & MoranScatterMean ~= 3) = 2; % If anything else at any point in time
        MoranScatterClusters(~leaf) = NaN; % Remove zeros due to NaN values.
        
        % Identify cells that are significant at final time
        Zt = zeros(size(Significance(:,:,1)));
        SignficanceAtT =  nanmean(Significance,3); % If cells significant for all time are desired: remove "(:,:,end)" and average through time e.g. nanmean(Significance,3).
        Zt(abs(SignficanceAtT) == 1) = 1; % One  if significant
        Zt(abs(SignficanceAtT) == 0) = 0; % Else Zt = 0
        Zt(~leaf) = NaN; %Remove zeros due to NaN values
        
        
        
        
        
        %% Condense Moran Scatter Plot
        [x,y] = prepareCurveData(nanmean(MoranScatterX,1)',nanmean(MoranScatterY,1)'); % Prepare curve data
        f=fit(x,y,'poly1'); % Fit to linear model
        figure('NumberTitle', 'off', 'Name', 'Figure 3','Position',[0 0 1000 500]);
        subplot(1,2,1)
        plot((-5:5),zeros(length(-5:5)),'--k') % Mark y-axis
        hold on
        plot(zeros(length(-5:5)),(-5:5),'--k') % Mark x-axis
        errorbar(nanmean(MoranScatterX,1),nanmean(MoranScatterY,1),nanmean(MoranScatterY,1)-nanmin(MoranScatterY,[],1),nanmean(MoranScatterY,1)-nanmax(MoranScatterY,[],1),nanmean(MoranScatterX,1)-nanmin(MoranScatterX,[],1),nanmean(MoranScatterX,1)-nanmax(MoranScatterX,[],1),'Marker','none','LineStyle','none','Color',[0 0 0]); % Plot time average scatter plot data with errorbars
        scatter(nanmean(MoranScatterX,1),nanmean(MoranScatterY,1),'r','filled') % Overlay scatter points so error bars are behind
        xlim([-1.5 2]) % Set consistent limits for all leaves
        ylim([-1.5 2]) % Set consistent limits for all leaves
        plot((-5:5),nanmean(f.p1)*(-5:5)+f.p2,'k','LineWidth',1.5) % Plot time averaged Moran's I
        text(1.1, 1.75,'\textrm{Q1}','interpreter','latex','FontSize',25) % Quadrant 1 label
        text(1.1,-1.25,'\textrm{Q4}','interpreter','latex','FontSize',25) % Quadrant 2 label
        text(-0.6, -1.25,'\textrm{Q3}','interpreter','latex','FontSize',25) % Quadrant 3 label
        text(-0.6, 1.75,'\textrm{Q2}','interpreter','latex','FontSize',25) % Quadrant 4 label
        xticks([-1 0 1 2]) % Set x tick values
        yticks([-1 0 1 2]) % Set y tick values
        d = get(gca,'XTickLabel'); % Find the x tick label
        set(gca,'XTickLabel',d,'fontsize',16) % Increase size of xticklabels
        ta=title('\textrm{\textbf{(a)}}','interpreter','latex','FontSize',30); % Create (a) subplot label
        ta.Position(1) = -2; % Adjust x coordinate of (a)
        ta.Position(2) = 1.9; % Adjust y coordinate of (a)
        xlabel('$d(\theta_i, \bar{\theta})$','interpreter','latex','FontSize',30) % Label x-axis
        ylabel('$\sum_j w_{ij}d(\theta_j, \bar{\theta})$','interpreter','latex','FontSize',30) % Label y-axis
        
        d = sign((reshape(nanmean(MoranScatterX,1),size(leaf)))).*...
            (abs(1-((reshape(nanmean(MoranScatterX,1),size(leaf)))./(reshape(nanmean(MoranScatterY,1),size(leaf)))))); % How close is point to line of gradient one?
        subplot(1,2,2)
        ax1=subplot(1,2,2); % Create axis 1 for one plot
        ax2=subplot(1,2,2); % Create axis 2 for overlay
        set(ax1,'Visible','off'); % Set ax1 visibility
        set(ax2,'Visible','off'); % Set ax1 visibility
        n =80; % Set number of increments in colorbar
        cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]'; % Set blue half of colourbar
        cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]'; % Set red half of colourbar
        cmap = [cmap1; cmap2(2:end, :)]; % Combine for complete colourmap
        
        p1=imagesc(ax2,d); % Plot d variable as image on ax2
        colormap(cmap) % Use custom colourmap
        caxis([-0.1 0.1]) % Colormap has limits [-0.1 0.1]
        h = colorbar('XTick', -0.1:0.05:0.1,'FontSize',15); %Set colorbar increments and fontsize
        slabel=ylabel(h, '$\bar{S}_{i}$','fontsize',30,'interpreter','latex'); % Set colourbar label
        slabel.Position(1) = 2.75; % Shift colourbar label to left
        hold on
        set(p1,'AlphaData',~isnan(MoranScatterClusters)) % Set NaN to white.
        im= []; % initialise image
        for i = 1:3
            temp = MoranScatterClusters==2;
            temp(~(MoranScatterClusters==2)) = 1; % Set to one if cells are in Q1 or Q3
            temp = temp-1;
            im(:,:,i) = temp; % Im is rgb array which is black for all cells in Q2 or Q4
        end
        p2 = imagesc(im); % Plot im
        set(p2,'AlphaData',MoranScatterClusters==2)
        set(gca,'xtick',[]) % Remove x ticks
        set(gca,'ytick',[]) % Remove y ticks
        xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
        ylim([y1-1 y2+1])
        tb=title('\textrm{\textbf{(b)}}','interpreter','latex','FontSize',30); % Add b label
        tb.Position(1) = x1-5; % Adjust x coordinate of b label
        tb.Position(2) = y1; % Adjust y coordinate of b label
        
        
        
        
        
        
        
        
        %% Global Measure Plot
        figure
        set(gcf,'Position',[0 0 800 800])
        % Second Row
        plot((1:length(KuramotoOrderParameter))*2/3,KuramotoOrderParameter,'LineWidth',2) % plot Kuramoto order parameter vs. time (increments of 2/3 hours)
        hold on
        plot((1:length(MoransI))*2/3,MoransI,'LineWidth',2)% plot global Moran's I vs. time (increments of 2/3 hours)
        plot((1:length(GearyC))*2/3,GearyC,'LineWidth',2)% plot global Geary's C vs. time (increments of 2/3 hours)
        xticks([0:24:length(MoransI)*2/3]) % Xticks in 24 hour intervals
        xlim([0 72]) % Set x lim to 72 hours
        ylim([0 1]) % Each measure is generally between 0 and 1.
        a = get(gca,'XLabel');
        set(gca,'XLabel',a,'fontsize',20) % Set font size for ticks
        yticks([0 0.1 0.6 0.8 0.9 1])   % Set y ticks
        breakyaxis([0.15,0.7]) % Insert axis break
        text(28,-0.1,'\textrm{Time (h)}','interpreter','latex','FontSize',30) % x-axis is in hours
        text(-7.5,0.25,'$r(t)$, $I^{\theta}(t)$ \textrm{or} $C^{\theta}(t)$','interpreter','latex','FontSize',30,'rotation',90)  % Insert y label as text (problem with ylabel and breakyaxis)
        
        
        %% Time Averaged local Moran, local Geary and dtheta plots
        figure
        set(gcf,'Position',[0 0 2000 800]) % Set size of figure
        subplot(1,3,1)
        h1=imagesc(nanmean(LocalMI,3));
        set(h1,'AlphaData',leaf)
        set(gca,'xtick',[]) % Remove x ticks
        set(gca,'ytick',[]) % Remove y ticks
        xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
        ylim([y1-1 y2+1])
        caxis([0 1])
        colormap((parula)); % Colourmap is parula
        hL=colorbar('fontsize',25,'Location','southoutside');
        ylabel(hL, '${I}^{\theta}_{i}$','fontsize',40,'interpreter','latex') % Label colourbar
        Clusters = (Zt == 1); % Clusters are made up of significant cells
        L = bwlabel(Clusters); % Label clusters
        % Plot clusters
        for j = 1:max(L,[],'all')
            Clusters = (L==j); % Find current cluster in L
            binaryImage = Clusters; % binary image for cluster location
            boundaries = bwboundaries(binaryImage); % Identify cluster boundary
            x = boundaries{1}(:, 2); % x coord of boundary
            y = boundaries{1}(:, 1); % y coord of boundary
            hold on
            plot(x,y,'color','red','LineWidth',4) % plot each of the clusters
        end
        hold on
        % Subplot 2: Local Geary C
        s2 = subplot(1,3,2);
        cla
        h1=imagesc(nanmean(LocalGC,3));
        set(h1,'AlphaData',leaf)
        set(gca,'xtick',[]) % Remove x ticks
        set(gca,'ytick',[]) % Remove y ticks
        xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
        ylim([y1-1 y2+1])
        colormap(s2,flipud(parula)); % Colourmap is opposite for local Geary
        caxis([0 0.2])
        hL=colorbar('fontsize',25,'Location','southoutside');
        ylabel(hL, '${C}^{\theta}_{i}$','fontsize',40,'interpreter','latex') % Label colourbar
        hold on
        % Subplot 3: d(theta_i,\bar(theta))
        subplot(1,3,3)
        cla
        h1=imagesc(nanmean(devphase,3));
        set(h1,'AlphaData',leaf)
        set(gca,'xtick',[]) % Remove x ticks
        set(gca,'ytick',[]) % Remove y ticks
        xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
        ylim([y1-1 y2+1])
        caxis([-1 1])
        colormap((parula));
        hL=colorbar('fontsize',25,'Location','southoutside');
        ylabel(hL, '$Time\, Averaged\, d(\theta_{i},\bar{\theta})$','fontsize',40,'interpreter','latex') % Label colourbar
        
        
        
        %% Video of local Moran's I, local Geary's C and d(theta_i,\bar{\theta}) over the whole time course
        figure(200 + leafn)
        VidName = ['Video_',leaffile]; % Define video name = Video + leaf number
        v = VideoWriter(VidName); % Ddefine file with name = VidName
        v.Quality = 100; % Set quality of video file
        v.FrameRate = 5; % Set frame rate for playback
        open(v) % open v, ready for writing
        set(gcf,'Position',[0 0 2000 800]) % Set size of figure
        clear Movie % Clear previous movie variable
        n = 1; % Set counter
        for i = 1:108
            cla
            % Subplot 2: Local Geary C
            subplot(1,3,1)
            h1=imagesc(LocalMI(:,:,i));
            set(h1,'AlphaData',leaf)
            set(gca,'xtick',[]) % Remove x ticks
            set(gca,'ytick',[]) % Remove y ticks
            xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
            ylim([y1-1 y2+1])
            caxis([0 1])
            colormap((parula)); % Colourmap is parula
            hL=colorbar('fontsize',25,'Location','southoutside');
            ylabel(hL, '${I}^{\theta}_{i}$','fontsize',40,'interpreter','latex') % Label colourbar
            hold on
            % Subplot 2: Local Geary C
            s2 = subplot(1,3,2);
            cla
            h1=imagesc(LocalGC(:,:,i));
            set(h1,'AlphaData',leaf)
            set(gca,'xtick',[]) % Remove x ticks
            set(gca,'ytick',[]) % Remove y ticks
            xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
            ylim([y1-1 y2+1])
            colormap(s2,flipud(parula)); % Colourmap is opposite for local Geary
            caxis([0 0.2])
            hL=colorbar('fontsize',25,'Location','southoutside');
            ylabel(hL, '${C}^{\theta}_{i}$','fontsize',40,'interpreter','latex') % Label colourbar
            hold on
            title(['\textrm{Time }= ',num2str(round(24 + (i)*2/3,2)),' \textrm{hours}'],'Fontsize',30,'interpreter','latex')
            % Subplot 3: d(theta_i,\bar(theta))
            subplot(1,3,3)
            cla
            h1=imagesc(devphase(:,:,i));
            set(h1,'AlphaData',leaf)
            set(gca,'xtick',[]) % Remove x ticks
            set(gca,'ytick',[]) % Remove y ticks
            xlim([x1-1 x2+1]) % Set limits based on x1,x2 and y1,y2 calculated previously
            ylim([y1-1 y2+1])
            caxis([-1 1])
            colormap((parula));
            hL=colorbar('fontsize',25,'Location','southoutside');
            ylabel(hL, '$d\theta_{i}$','fontsize',40,'interpreter','latex') % Label colourbar
            Movie(n) = getframe(gcf); % Current movie frame is the entire figure
            n = n + 1;
        end
        
        writeVideo(v,Movie) % Write video footage for future reference
        
        close(v) % close video file
        
        
        
        
        
        
        
    else
        
        disp([leaffile,' is not located in the current directory. To download this dataset please visit:'])
        disp('https://www.research.ed.ac.uk/portal/en/datasets/data-from-spontaneous-spatiotemporal-waves-of-gene-expression-from-biological-clocks-in-the-leaf-wenden-toner-et-al-pnas-2012-109-67576762(904f8897-1eda-4737-a3fa-e5e08f078a49).html')
        
        
        
    end
    
    
else
    disp('Matlab functions breakxaxis and breakyaxis are not present in the current directory')
    disp('For breakxaxis.m, please visit: https://www.mathworks.com/matlabcentral/fileexchange/42905-break-x-axis')
    disp('For breakyaxis.m, please visit: https://www.mathworks.com/matlabcentral/fileexchange/45760-break-y-axis')
    
end

function [KuramotoOrderParameter] = KuramotoOP(phases)
%% Definition: Kuramoto order parameter for population of phase oscillators. This will give an indication of the synchrony at a certain time.
%% Input:
% Phases: Phase matrix = snapshot of oscillatory phases at a certain
%         time point
%% Output:
% KuramotoOrderParameter: Kuramoto order parameter, KuramotoOrderParameter. KuramotoOrderParameter = 1 is perfect synchrony. KuramotoOrderParameter = 0 is complete asynchrony.


%% Code
% What are the dimensions of the data?
Ny = length(phases(:,1));
Nx = length(phases(1,:));


% Initialise x and y components of phases.
xvals=zeros(1,nansum(~isnan(phases),'all'));
yvals=zeros(1,nansum(~isnan(phases),'all'));
n=1; % Initialise counter
for j = 1:Nx
    for m = 1:Ny
        if isnan(phases(m,j))==0 % If phase is not from outside of leaf
            xvals(n) = cos(phases(m,j)); % x = cos(theta)
            yvals(n) = sin(phases(m,j)); % y = sin(theta)
            n=n+1; % Updata counter
        end
    end
end

% Compute average x and y
meanx = mean(xvals);
meany = mean(yvals);
% KuramotoOrderParameter is the length of the vector v = meanx (i) + meany (j)
KuramotoOrderParameter = (meanx^2 + meany^2)^0.5;
end



function [LocalGC, GearyC] = localGeary(phases,radius)
%% Definition: Function for the calulation of local Geary's GearyC distribution from a snapshot of oscillatory phase at a single time point.
%% Input:
% Phases: This a matrix containing the phase values from a single time point.
% radius: radius used to generate the spatial weight matrix. e.g. if
%         radius if 1, cells within one cell of eachother are said to
%         interact.
%% Output:
% LocalGC: local version of Geary's GearyC. Gives a direct indication of cells
%     that are coherent with respect to their neighbours.
% GearyC: Global Geary's GearyC. For our variant of local Geary's GearyC, this is the
%    average of the local Geary's GearyC distribution.


% Which cells are part of the leaf?
leaf = ~isnan(phases);


x = cos(phases); % x component of phase values
y = sin(phases); % y component of phase values



K = ones(2*radius + 1,2*radius + 1); % Kernel used for spatial weight computation
K(radius+1,radius+1) = 0; % Cells are not coupled to themselves to centre of kernel is 0.


phasehat = angle(nanmean(x,'all')+nanmean(y,'all')*i); % Compute average phase as the argument of the complex number given by the mean x and y coords.


devphase = angle(cos(phases-phasehat)+sin(phases-phasehat)*i).*leaf;  % Compute deviation of phase from average using distance measure. This is the arc length separating the phase and the population average.

phases(isnan(phases))=0; % Set NaN values from outside of leaf to zero: Prevents NaN in convolution.


w = conv2(leaf,K,'same');  % How many cells is each cell coupled to?
W = nansum(w,'all'); % W = sum of the spatially weighted matrix over all indices.


phases(~leaf)=NaN; % As we do not use convolution for this, phases outside leaf can be set to NaN again.
sumtop = zeros(size(phases)); % Initialise sum for top of local Geary calulcation

for k = 1:length(phases(1,:))
    for j = 1:length(phases(:,1))
        % For each cell in phases
        for m = max(k-radius,1):min(k+radius:length(phases(1,:)))
            for n = max(j-radius,1):min(j+radius:length(phases(:,1)))
                % For cells within a radius given by input
                if n~=m && ~isnan(phases(n,m)) && ~isnan(phases(j,k))
                    % n~m: w = 0 for i~=j,
                    % ~isnan(phases(n,m)): Don't take cells from outside leaf
                    % ~isnan(phases(j,k)): Do not conside cells outside of leaf
                    
                    % Sum up values of the arcelngth separating the distance measures between i and j
                    sumtop(j,k) = sumtop(j,k) + (angle(cos(devphase(j,k)-devphase(n,m))+sin(devphase(j,k)-devphase(n,m))*i)^2)/w(j,k);
                end
            end
        end
    end
end


sumtop = sumtop.*leaf; % Remove any values outside of leaf
sumtop(~leaf)=NaN; % Set cells outside of leaf to NaN
Ncur = nansum(leaf,'all'); % Number of cells in leaf

sumbot = nansum((devphase).^2,'all'); % Sum on bottom of calculation = sum of squared deviations from average
LocalGC = ((Ncur-1)*sumtop)./((2)*(sumbot)); % Bring all components together for final calculation of local Geary GearyC
GearyC = nanmean(LocalGC,'all'); % Global GearyC is the average across all cells.
end


function [Index,devphase] = LocalMoran2(phases,radius)
%% Description: Local Moran's MoransI function. Compute a local measure of spatial autocorrelation snapshot of oscillatory phases.
%% Inputs:
% phases: This a matrix containing the phase values from a single time point.
% radius: radius used to generate the spatial weight matrix. e.g. if
%         radius if 1, cells within one cell of eachother are said to
%         interact.
%% Outputs:
% Index: Index is a matrix containing the local Moran values for the
%        inputted phase matrix.
% devphase: This is the distance measure computed for each cell. e.g.
%           the arc length separating the phase of the cell and the
%           population average. This is used later when computing the
%           Moran scatter plots.

%% Code:
x = cos(phases); % Compute x component of phases
y = sin(phases); % Compute y component of phases


K = ones(2*radius + 1,2*radius + 1); % Generate kernel based on convolution: used in spatially weighted matrix calculations
K(radius+1,radius+1) = 0; % Cells are not coupled to themselves to centre of kernel is 0.

phasehat = angle(nanmean(x,'all')+nanmean(y,'all')*i); % Compute average phase as the argument of the complex number given by the mean x and y coords.





leaf = ~isnan(phases); % Identify which cells are part of the leaf.
phases(isnan(phases))=0; % Set NaN values from outside of leaf to zero: Prevents NaN in convolution.



w = conv2(leaf,K,'same').*leaf; % How many cells is each cell coupled to?

devphase = (angle(cos(phases-phasehat)+sin(phases-phasehat)*i)).*leaf; % Compute deviation of phase from average using distance measure. This is the arc length separating the phase and the population average.


W = nansum(w,'all'); % W = sum of the spatially weighted matrix over all indices.

sumbot = (devphase).^2; % Sum on denominator of local Moran computation

% Compute m2 = (sum(devphase^2)))
m2 = nansum(sumbot,'all');
Ncur = nansum(leaf,'all');


% Compute numerator in local Moran computation
w(~leaf) = 0;
sumtop = (conv2((devphase),K,'same').*leaf./w).*leaf.*devphase;
sumtop(~leaf) = NaN;

% Final calculation.
Index = (Ncur)*sumtop./(m2);
end




