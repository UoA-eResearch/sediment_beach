
%This script reconstructs shoreline data using all SVD modes except the first mode
site=2;
addpath('C:\Users\mpul348\Documents\code_paper1')
start_date_request = datetime(1999,08,01,0,0,0); % CoastSat data starts
end_date_request   = datetime(2024,12,30,0,0,0); % Coastsat data ends
n_d = days(end_date_request - start_date_request);
n_days = days(end_date_request - start_date_request);

if site==1
    n=11;
elseif site==2
    n=16;
end

%preallocating variables
yint     = cell(n,1);
yy       = cell(n,1);
ysave    = cell(n,1);
ysave2   = cell(n,1);
day      = cell(n,1);
time     = cell(n,1);
real_time= cell(n,1);
t        = cell(n,1);
s        = cell(1,n); 
n_missing = cell(1,n); 
yfilled = cell(1,n);
U = cell(1,n); 
S = cell(1,n); 
V = cell(1,n); 
singvals = cell(1,n); 
power_each = cell(1,n); 
power_total = cell(1,n); 
explained_each = cell(1,n); 
PC_lst = cell(1,n); 
explained_without_mode1 = cell(1,n); 
explained_rem = cell(1,n); 

%this loop perform the PCA
for k = 1:n
    [s{k}, t1{k}, lat_s{k}, lon_s{k}, tit{k}, sea_long{k}, sea_lat{k}, center_x{k}, center_y{k}]=load_data(k,site);%load and stores the kth results in the kth cell of s, keep all 11 in a cell array so s{3} = s for third beach
    [s{k}] = tr_selection(s{k},site,k);
    [yint{k}, yy{k}, ysave{k}, ysave2{k}, day{k}, time{k}, real_time{k}, t{k}] = smoothgap2(s{k}, t1{k}); %builds a daily time series for every transect
    [Mmean{k}, file_mean{k}, ymean{k}, annual_mean{k}, file_series{k}, coeff{k}, yhat{k}, mint{k}] = shoreline_metrics(yint{k}, time{k}, day{k});
    n_missing{k} = sum(isnan(ymean{k}(:)));
    yfilled{k} = fillmissing(ymean{k}, 'constant', 0);
    [U{k}, S{k}, V{k}] = svd(yfilled{k}, 'econ'); %U is the temporal engeinfunction, S is the singular values and V the spatial eigenfunction
    singvals{k} = diag(S{k});    %select just the diagonal of the S matrix
    power_each{k} = singvals{k}.^2; %Variance explained by each mode
    power_total{k} = sum(power_each{k});  %total variance  = sum(singvals.^2)
    explained_each{k} = power_each{k} / power_total{k};        %fraction explainde by each mode
    PC_lst{k} = U{k}(:,2:end) * S{k}(2:end,2:end) * V{k}(:,2:end)';    
end

%this loop calculate mean per transect and repeat the mean along time
mean_transect = cell(1,n);
mean_PC_lst = cell(1,n);
for k = 1:n
    mean_transect{k} = mean(PC_lst{k}, 1, 'omitnan');
     mean_PC_lst{k} = repmat(mean_transect{k}, size(PC_lst{k},1), 1);
end

% calculating the Variance explained by the PCA, after excluding mode 1
for k = 1:n
    explained_without_mode1{k} = ...
        power_each{k}(2:end) / sum(power_each{k}(2:end));
    explained_rem{k} = ...
        sum(power_each{k}(2:end)) / power_total{k};
end
save('PC_lst_site2.mat', 'PC_lst');
save('mean_PC_lst_site2.mat', 'mean_PC_lst')
save('explained_rem_site2.mat', 'explained_rem')
% %mode1
% figure('Position',[50 50 1600 900])
% tiledlayout(4,4,'TileSpacing','compact','Padding','compact')
% for k = 1:length(V)
%     nexttile
%     mode1 = V{k}(:,1);
%     x = 1:length(mode1);
%     hold on
%     for i = 1:length(mode1)-1
%         if mode1(i) >= 0
%             c = 'b';
%         else
%             c = 'r';
%         end
%         plot(x(i:i+1), mode1(i:i+1), ...
%             'Color', c, 'LineWidth', 1.5)
%     end
%     yline(0,'k-')
%     xlabel('Transect')
%     ylabel('Mode 1')
%     title(sprintf('Beach %d (%.2f%%)', ...
%         k, 100*explained_each{k}(1)))
%     grid on
%     hold off
% end
% 
% %%%
% %mode 2
% figure('Position',[50 50 1600 900])
% tiledlayout(4,4,'TileSpacing','compact','Padding','compact')
% for k = 1:length(V)
%     nexttile
%     mode2 = V{k}(:,2);
%     x = 1:length(mode2);
%     hold on
%     for i = 1:length(mode2)-1
%         if mode2(i) >= 0
%             c = 'b';
%         else
%             c = 'r';
%         end
%         plot(x(i:i+1), mode2(i:i+1), ...
%             'Color', c, 'LineWidth', 1.5)
%     end
%     yline(0,'k-')
%     xlabel('Transect')
%     ylabel('Mode 2')
% 
%     title(sprintf('Beach %d (%.2f%%)', ...
%         k, 100*explained_each{k}(2)))
%     grid on
%     hold off
% end
% %%%%
