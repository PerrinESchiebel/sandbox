npts = 15;

x = data.data(:,3:3:end); % horizontal plane
y = data.data(:,4:3:end); % vertical direction
z = data.data(:,5:3:end); % horizontal plane
% xo = x;yo = y;zo = z;
%Take out points that are there for less than five frames
indfins = sum(isfinite(x))>5;
xm = x(:,indfins);
zm = z(:,indfins);
ym = y(:,indfins);
clear indfins

%%
%%%Have user find wall points
indrail = nan(1,2);
plot(xm(1,:),zm(1,:),'o');hold on;
[xrail,zrail] = getpts;
plot(xrail,zrail,'o','MarkerFaceColor',[78 52 46]./255,'MarkerSize',12);hold off;pause(0.5);
[~,ind] = min(sqrt((xm(1,:)-xrail(1)).^2+(zm(1,:)-zrail(1)).^2));
indrail(1) = ind;
[~,ind] = min(sqrt((xm(1,:)-xrail(2)).^2+(zm(1,:)-zrail(2)).^2));
indrail(2) = ind;

xrail = xm(:,indrail);
zrail = zm(:,indrail);
yrail = ym(:,indrail);

zm(:,indrail) = [];
xm(:,indrail) = [];
ym(:,indrail) = [];
plot(xm,zm,'.k');hold on;plot(xrail,zrail,'s','MarkerFaceColor',[1 0 0.8]);
clear indrail ind;
%%
%%%remove post markers from data and store in x/y/zpegs. Do after getting
%%%the rail since sometimes those points are within this area
pegrange = [-0.2 0.4 0.2 0.5]; %%[xmin zmin xmax zmax]
rangeind = xm(1,:)>pegrange(1) & xm(1,:)<pegrange(3) & zm(1,:)>pegrange(2) & zm(1,:) < pegrange(4);
xpegs = xm(:,rangeind);
zpegs = zm(:,rangeind);
ypegs = ym(:,rangeind);
[~,ix] = sort(xpegs(1,:),'ascend');
xpegs = xpegs(:,ix);
zpegs = zpegs(:,ix);
ypegs = ypegs(:,ix);

zm(:,rangeind) = [];
xm(:,rangeind) = [];
ym(:,rangeind) = [];
clear rangeind pegrange ix
%%
%%%have user pick out which snake midlines they want to look at. This will
%%%take only midlines that are completely present within the rectangle
plot(xm,zm,'.k');hold on;plot(xpegs,zpegs,'o','MarkerFaceColor',[0 1 1]);plot(xrail,zrail,'h','MarkerFaceColor',[1 0 1]);
ax1 = gca;axis([min(min(xm,[],'omitnan'))-0.1,max(max(xm,[],'omitnan'))+0.1,min(min(zm,[],'omitnan'))-0.1,max(max(zm,[],'omitnan'))+0.1])
rect = getrect(ax1);
hold on;
rectangle('Position',rect)
xrange = [rect(1),rect(1)+rect(3)];
zrange = [rect(2),rect(2)+rect(4)];

firstfull = xm>xrange(1) & xm<xrange(2) & zm>zrange(1) & zm<zrange(2) & ym>yrange(1) & ym<yrange(2); %find everywhere all points are inside the rectangle
firstFrame = find(sum(firstfull,2)>npts,1,'first');  %grab the first time all points get in the rectangle
lastFrame = find(sum(firstfull,2)>npts,1,'last');    %grab the last time all points get in the rectangle

plot(xm(lastFrame,:),zm(lastFrame,:),'.c')
plot(xm(firstFrame,:),zm(firstFrame,:),'.c')

xm = xm(firstFrame:lastFrame,:);zm = zm(firstFrame:lastFrame,:);ym = ym(firstFrame:lastFrame,:);

outofbounds = xm<xrange(1) | xm>xrange(2) | zm<zrange(1) | zm>zrange(2) | ym<yrange(1) | ym>yrange(2); %get rid of everything outside the rectangle
xm(outofbounds == 1) = NaN;
zm(outofbounds == 1) = NaN;
ym(outofbounds == 1) = NaN;
clear ax1 rect xrange zrange outofbounds firstfull

%%
% headMarker = find(zm(firstFrame,:)==max(zm(firstFrame,:)));
% tailMarker = find(zm(firstFrame,:)==min(zm(firstFrame,:)));
% finitepts = isfinite(xm(1,:));
% tempx = xm(1,finitepts);tempz = zm(1,finitepts);
% alldist = pdist([tempx',tempz']);
% alldist = squareform(alldist);
% nfinitepts = sum(finitepts);
% [~,sortindex] = sort(tempz,'descend');
% % nmissingtotal = sum(isnan(xtMinus1));
% sortindex = nan(nfinitepts,1);
% sortindex(1) = headMarker;
% % display(ii)
% % if ii==3
% %     keyboard
% % end
% for mm=1:npts
%         alldist(mm,mm) = Inf;
% end
% for rr=1:nfinitepts-1
%     nextpoint = find(alldist(sortindex(rr),:)==min(alldist(sortindex(rr),:)));
%     alldist(:,sortindex(rr)) = Inf;
%     sortindex(rr+1) = nextpoint;
% %     imagesc(alldist);drawnow;pause(0.5);
% end
% xm = xm(:,[sortindex;(nfinitepts+1:(size(ym,2)))']);
% zm = zm(:,[sortindex;(nfinitepts+1:(size(ym,2)))']);
% ym = ym(:,[sortindex;(nfinitepts+1:(size(ym,2)))']);
%      ptx = fillmissing(xm,'spline',2);
%     ptz = fillmissing(zm,'spline',2);
%     pty = fillmissing(ym,'spline',2);
%%
numtimes = size(zm,1);
alwayspoint = find(sum(isfinite(xm))==numtimes,1,'first');
sortx = cell(numtimes,1);
sortz = cell(numtimes,1);
sorty = cell(numtimes,1);
alwayspointtimeFirst = find(sum(isfinite(xm),2)==max(sum(isfinite(xm),2)),1,'first');

%%
finitepts = isfinite(xm(alwayspointtimeFirst,:));
tempx = xm(alwayspointtimeFirst,finitepts);tempz = zm(alwayspointtimeFirst,finitepts);
% tailMarker = find(tempz==min(tempz));
alldist = pdist([tempx',tempz']);
alldist = squareform(alldist);
nfinitepts = sum(finitepts);
sortindex = nan(nfinitepts,1);
 sortindex(1) =  find(zm(alwayspointtimeFirst,:)==max(zm(alwayspointtimeFirst,:)));
for mm=1:nfinitepts
        alldist(mm,mm) = Inf;
end
for rr=1:nfinitepts-1
    nextpoint = find(alldist(sortindex(rr),:)==min(alldist(sortindex(rr),:)));
    alldist(:,sortindex(rr)) = Inf;
    sortindex(rr+1) = nextpoint;
%     imagesc(alldist);drawnow;pause(0.5);
end
alwayspointRealPos = find(sortindex==alwayspoint);
%%

% for tt=1:numtimes
for tt = 124
% headMarker = find(zm(tt,:)==max(zm(tt,:)));
nanbefore = sum(isnan(xm(tt,1:alwayspoint-1)));
finitepts = isfinite(xm(tt,:));
tempx = xm(tt,finitepts);tempz = zm(tt,finitepts);tempy = ym(tt,finitepts);
finitepts = isfinite(xm(tt,:));
% tempx = xm(allwayspointtimeFirst,finitepts);tempz = zm(allwayspointtimeFirst,finitepts);
% % tailMarker = find(tempz==min(tempz));

% maxdistanceallowed = mean(min(alldist)) + std(min(alldist));
nfinitepts = sum(finitepts);
% % [~,sortindex] = sort(tempz,'descend');
% % nmissingtotal = sum(isnan(xtMinus1));
sortindex = nan(nfinitepts,1);

alwayspointthistime = alwayspoint-nanbefore;
alwayspointRealPosthistime = alwayspointRealPos - nanbefore;
sortindexAnt = nan(1,alwayspointRealPosthistime);
sortindexAnt(1) = alwayspointthistime;
sortindexPost = nan(1,nfinitepts - alwayspointRealPosthistime+1);
sortindexPost(1) = alwayspointthistime;
% sortindex(alwayspointthistime) = alwayspoint;
% count = 1;
alldist = pdist([tempx',tempz']);
alldist = squareform(alldist);
for mm=1:nfinitepts
        alldist(mm,mm) = Inf;
end
% AntDist = alldist(1:alwayspointthistime,1:alwayspointthistime);
% PostDist = alldist(alwayspointthistime:end,alwayspointthistime:end);
AntDist = alldist(1:alwayspointRealPosthistime,1:alwayspointRealPosthistime);
PostDist = alldist(alwayspointRealPosthistime:end,alwayspointRealPosthistime:end);
if alwayspointRealPosthistime == 1
else
for rr = 1:alwayspointRealPosthistime-1
%     currentpoint = sortindex(alwayspointthistime - count+1);
%     nextpoint = find(alldist(currentpoint,:)==min(alldist(currentpoint,:)));   

    nextpoint = find(AntDist(sortindexAnt(rr),:)==min(AntDist(sortindexAnt(rr),:)));
    AntDist(:,sortindexAnt(rr)) = Inf;
%     sortindex(rr+1) = nextpoint;
%     alldist(:,currentpoint) = Inf;
    sortindexAnt(rr+1) = nextpoint;
    
% count = count+1;
end
end
for rr = 1:nfinitepts-alwayspointRealPosthistime
%     currentpoint = sortindex(alwayspointthistime - count+1);
%     nextpoint = find(alldist(currentpoint,:)==min(alldist(currentpoint,:)));   

    nextpoint = find(PostDist(sortindexPost(rr),:)==min(PostDist(sortindexPost(rr),:)));
    PostDist(:,sortindexPost(rr)) = Inf;
%     sortindex(rr+1) = nextpoint;
%     alldist(:,currentpoint) = Inf;
    sortindexPost(rr+1) = nextpoint;
%     imagesc(alldist);drawnow;pause(0.5);
% count = count+1;
end
% while count < nfinitepts
%     currentpoint = sortindex(alwayspointthistime - count+1);
%     nextpoint = find(alldist(currentpoint,:)==min(alldist(currentpoint,:)));   
%     alldist(:,currentpoint) = Inf;
%     sortindex(alwayspointthistime - count) = nextpoint;
% %     imagesc(alldist);drawnow;pause(0.5);
% count = count+1;
% end
sortindex(1:alwayspointRealPosthistime) = flipud(sortindexAnt);
sortindex(alwayspointRealPosthistime:end) = sortindexPost;
sortx{tt} = tempx(sortindex);
sortz{tt} = tempz(sortindex);
sorty{tt} = tempy(sortindex);
plot(x(tt+firstFrame-1,:),z(tt+firstFrame-1,:),'o','MarkerFaceColor','k');hold on;
% plot(xm(tt,:),zm(tt,:),'o','Color',[0.8 0 0.8],'LineWidth',2);
scatter(sortx{tt},sortz{tt},50,1:length(sortx{tt}),'o','filled');
colormap(hsv);hold off;axis equal tight;drawnow;
end
%%
% numtimes = size(zm,1);
% sortx = cell(numtimes,1);
% sortz = cell(numtimes,1);
% sorty = cell(numtimes,1);
% for tt=1:numtimes
% % headMarker = find(zm(tt,:)==max(zm(tt,:)));
% 
% finitepts = isfinite(xm(tt,:));
% tempx = xm(tt,finitepts);tempz = zm(tt,finitepts);tempy = ym(tt,finitepts);
% % tailMarker = find(tempz==min(tempz));
% alldist = pdist([tempx',tempz']);
% alldist = squareform(alldist);
% nfinitepts = sum(finitepts);
% % [~,sortindex] = sort(tempz,'descend');
% % nmissingtotal = sum(isnan(xtMinus1));
% sortindex = nan(nfinitepts,1);
% sortindex(1) = tailMarker;
% % display(ii)
% % if ii==3
% %     keyboard
% % end
% for mm=1:nfinitepts
%         alldist(mm,mm) = Inf;
% end
% for rr=1:nfinitepts-1
%     nextpoint = find(alldist(sortindex(rr),:)==min(alldist(sortindex(rr),:)));
%     alldist(:,sortindex(rr)) = Inf;
%     sortindex(rr+1) = nextpoint;
% %     imagesc(alldist);drawnow;pause(0.5);
% end
% % nfps(tt) = nfinitepts;
% sortx{tt} = tempx(sortindex);
% sortz{tt} = tempz(sortindex);
% sorty{tt} = tempy(sortindex);
% plot(x(tt+firstFrame,:),z(tt+firstFrame,:),'o','MarkerFaceColor','k');hold on;
% plot(xm(tt,:),zm(tt,:),'o','Color',[0.8 0 0.8],'LineWidth',2);
% scatter(sortx{tt},sortz{tt},50,1:length(sortx{tt}),'o','filled');hold off;axis equal tight;drawnow;
% end
% % xm = xm(:,[sortindex;(nfinitepts+1:(size(ym,2)))']);
% % zm = zm(:,[sortindex;(nfinitepts+1:(size(ym,2)))']);
% % ym = ym(:,[sortindex;(nfinitepts+1:(size(ym,2)))']);
% %      ptx = fillmissing(xm,'spline',2);
% %     ptz = fillmissing(zm,'spline',2);
% %     pty = fillmissing(ym,'spline',2);
%%
% lastpoint = 19;
numsplinepts = 50;
ptx = nan(numsplinepts,numtimes);
ptz = nan(numsplinepts,numtimes);
for tt = 1:numtimes
pt = interparc(numsplinepts,sortx{tt},sortz{tt},'spline');
ptx(:,tt) = pt(:,1);
ptz(:,tt) = pt(:,2);
end
% Curvature = spatialCurvature(ptz,ptx,5,0.01,1);
% Curvature = spatialCurvature(ptz,ptx,5);