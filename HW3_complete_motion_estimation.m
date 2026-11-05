%% ECE5578 Homework 3 - Motion Estimation (all questions)
% Based on the supplied motion_estimation.m and Li et al., TCSVT 2018.
% Input: frames/Beauty/Beauty-01.png ... Beauty-05.png, and likewise
% Bosphorus and Jockey. Alternatively place Beauty.mp4 etc. in pwd.
% Outputs: figures (PNG) and measurement tables (CSV) in hw3_results.
% Works with base MATLAB + Image Processing Toolbox for imresize/rgb2gray.
% Full-search at 720p is computationally demanding. Allow ample runtime.
clear; clc; close all;
cfg.names = {'Beauty','Bosphorus','Jockey'};
cfg.blockSizes = [8 4];
cfg.search = 16;              % displacements -16,...,+16 (33x33)
cfg.nFrames = 5;              % 4 consecutive P-frame pairs
cfg.outDir = fullfile(pwd,'hw3_results');
cfg.saveFigures = true;
cfg.runDiamond = true;
cfg.runAffine = true;
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
summary = table(); complexity = table(); affineSummary = table();
for seq = 1:numel(cfg.names)
    name = cfg.names{seq};
    frames = loadSequence(name,cfg.nFrames);
    fprintf('\n%s: %d frames at %dx%d\n',name,numel(frames),size(frames{1},2),size(frames{1},1));
    for t = 2:cfg.nFrames
        ref = frames{t-1}; cur = frames{t};
        for bs = cfg.blockSizes
            fprintf('%s %d->%d, block %dx%d: exhaustive...',name,t-1,t,bs,bs);
            [predI, mvI, countI] = exhaustiveME(cur,ref,bs,cfg.search);
            mseI = mean((cur(:)-predI(:)).^2);
            fprintf(' half-pel...');
            [predH, mvH] = halfRefineME(cur,ref,bs,mvI,cfg.search);
            mseH = mean((cur(:)-predH(:)).^2);
            row = table(string(name),t-1,t,bs,mseI,mseH,mseI-mseH,...
                100*(mseI-mseH)/max(mseI,eps),...
                'VariableNames',{'Sequence','ReferenceFrame','CurrentFrame','BlockSize','IntegerMSE','HalfPelMSE','DeltaMSE','ImprovementPercent'});
            summary = [summary;row]; %#ok<AGROW>
            fprintf(' Integer MSE %.5f, half MSE %.5f\n',mseI,mseH);
            if cfg.saveFigures
                showPanels(cur,predI,mvI,bs,name,t,'integer',cfg.outDir);
                showPanels(cur,predH,mvH,bs,name,t,'half',cfg.outDir);
            end
            if cfg.runDiamond && t<=4
                fprintf('   diamond search...');
                [predD,mvD,countD] = diamondME(cur,ref,bs,cfg.search);
                mseD = mean((cur(:)-predD(:)).^2);
                row = table(string(name),t-1,t,bs,mean(countI(:)),mean(countD(:)),...
                    100*(1-mean(countD(:))/mean(countI(:))),mseI,mseD,...
                    'VariableNames',{'Sequence','ReferenceFrame','CurrentFrame','BlockSize','ExhaustiveAvgSSD','DiamondAvgSSD','ReductionPercent','ExhaustiveMSE','DiamondMSE'});
                complexity = [complexity;row]; %#ok<AGROW>
                if cfg.saveFigures,showPanels(cur,predD,mvD,bs,name,t,'diamond',cfg.outDir);end
                fprintf(' %.2f vs %d SSD/block\n',mean(countD(:)),round(mean(countI(:))));
            end
            if cfg.runAffine && t<=4 && bs==8
                fprintf('   affine four-parameter estimation...');
                [predA,mvA] = affineME(cur,ref,bs,mvI);
                mseA = mean((cur(:)-predA(:)).^2);
                affineSummary = [affineSummary;table(string(name),t-1,t,bs,mseI,mseA,mseI-mseA,...
                    'VariableNames',{'Sequence','ReferenceFrame','CurrentFrame','BlockSize','TranslationMSE','AffineMSE','DeltaMSE'})]; %#ok<AGROW>
                if cfg.saveFigures,showPanels(cur,predA,mvA,bs,name,t,'affine',cfg.outDir);end
                fprintf(' %.5f (translation %.5f)\n',mseA,mseI);
            end
        end
    end
end
writetable(summary,fullfile(cfg.outDir,'Q1_Q2_MSE.csv'));
writetable(complexity,fullfile(cfg.outDir,'Q3_search_counts.csv'));
writetable(affineSummary,fullfile(cfg.outDir,'Q4_affine_MSE.csv'));
disp('Q1 and Q2 results'); disp(summary);
disp('Q3 results'); disp(complexity);
disp('Q4 results'); disp(affineSummary);
fprintf('\nSaved all results under %s\n',cfg.outDir);

function frames=loadSequence(name,N)
    frames=cell(N,1);
    for k=1:N
        candidates={fullfile('frames',name,sprintf('%s-%02d.png',name,k)),...
            fullfile('frames',name,sprintf('%s-%d.png',name,k)),...
            sprintf('%s-%02d.png',name,k),sprintf('%s-%d.png',name,k)};
        found='';
        for p=1:numel(candidates)
            if isfile(candidates{p}),found=candidates{p};break;end
        end
        if ~isempty(found)
            a=imread(found);
        else
            videoNames={sprintf('%s.mp4',name),fullfile('frames',sprintf('%s.mp4',name))};
            video='';
            for v=1:numel(videoNames)
                if isfile(videoNames{v}), video=videoNames{v};break;end
            end
            if isempty(video),error('Missing frames for %s. See input comments at top of script.',name);end
            vr=VideoReader(video); a=read(vr,k);
        end
        if size(a,3)==3
            a=0.299*double(a(:,:,1))+0.587*double(a(:,:,2))+0.114*double(a(:,:,3));
        else
            a=double(a);
        end
        if size(a,1)~=720 || size(a,2)~=1280
            a=imresize(a,[720 1280]);
        end
        frames{k}=a;
    end
end

function [pred,mv,nTests]=exhaustiveME(cur,ref,bs,R)
    % All image blocks are compared in parallel for each candidate shift.
    [H,W]=size(cur); assert(mod(H,bs)==0 && mod(W,bs)==0);
    nr=H/bs; nc=W/bs;
    best=inf(nr,nc); dxBest=zeros(nr,nc);dyBest=zeros(nr,nc);
    nTests=zeros(nr,nc);
    yy=1:H; xx=1:W;
    for dy=-R:R
        ridx=min(max(yy+dy,1),H);
        for dx=-R:R
            cidx=min(max(xx+dx,1),W);
            d=(cur-ref(ridx,cidx)).^2;
            d=reshape(d,bs,nr,bs,nc);
            score=reshape(sum(sum(d,1),3),nr,nc);
            isBetter=score<best;
            best(isBetter)=score(isBetter);
            dxBest(isBetter)=dx;dyBest(isBetter)=dy;
            nTests=nTests+1;
        end
    end
    mv.dx=dxBest;mv.dy=dyBest;
    pred=compensateInteger(ref,bs,mv);
end

function pred=compensateInteger(ref,bs,mv)
    [H,W]=size(ref); nr=H/bs;nc=W/bs;pred=zeros(H,W);
    for r=1:nr
        rows=(r-1)*bs+(1:bs);
        for c=1:nc
            cols=(c-1)*bs+(1:bs);
            ry=min(max(rows+mv.dy(r,c),1),H);
            rx=min(max(cols+mv.dx(r,c),1),W);
            pred(rows,cols)=ref(ry,rx);
        end
    end
end

function [pred,mv]=halfRefineME(cur,ref,bs,mvInit,R)
    % The best integer solution plus its eight neighboring half-pel offsets.
    [H,W]=size(ref);nr=H/bs;nc=W/bs;
    pred=zeros(H,W);mv=mvInit;
    [xx,yy]=meshgrid(0:bs-1,0:bs-1);
    for r=1:nr
        row=(r-1)*bs+1; rows=row+(0:bs-1);
        for c=1:nc
            col=(c-1)*bs+1;cols=col+(0:bs-1);
            target=cur(rows,cols);bx=col+xx;by=row+yy;
            best=inf;dx0=mvInit.dx(r,c);dy0=mvInit.dy(r,c);
            for dv=[0 -.5 .5]
                for du=[0 -.5 .5]
                    dx=dx0+du;dy=dy0+dv;
                    if abs(dx)>R || abs(dy)>R,continue;end
                    patch=sampleBlock(ref,bx+dx,by+dy);
                    err=sum((target(:)-patch(:)).^2);
                    if err<best,best=err;pred(rows,cols)=patch;mv.dx(r,c)=dx;mv.dy(r,c)=dy;end
                end
            end
        end
    end
end

function [pred,mv,tests]=diamondME(cur,ref,bs,R)
    % Zhu & Ma's large-diamond search followed by a small-diamond refinement.
    [H,W]=size(cur);nr=H/bs;nc=W/bs;pred=zeros(H,W);
    mv.dx=zeros(nr,nc);mv.dy=zeros(nr,nc);tests=zeros(nr,nc);
    large=[0 0; 0 -2; 1 -1;2 0;1 1;0 2;-1 1;-2 0;-1 -1];
    small=[0 0;0 -1;1 0;0 1;-1 0];
    for r=1:nr
        rows=(r-1)*bs+(1:bs);
        for c=1:nc
            cols=(c-1)*bs+(1:bs); target=cur(rows,cols);
            center=[0 0]; memo=nan(2*R+1,2*R+1);
            % Iterate until the center has minimum SSD (or a safety cap).
            for iter=1:(2*R+5)
                [where,memo]=bestPattern(center,large,target,ref,rows,cols,R,memo);
                if isequal(where,center),break;end
                center=where;
            end
            [where,memo]=bestPattern(center,small,target,ref,rows,cols,R,memo);
            dx=where(1);dy=where(2);
            mv.dx(r,c)=dx;mv.dy(r,c)=dy;
            tests(r,c)=nnz(~isnan(memo));
            pred(rows,cols)=ref(min(max(rows+dy,1),H),min(max(cols+dx,1),W));
        end
    end
end

function [bestXY,memo]=bestPattern(center,pattern,target,ref,rows,cols,R,memo)
    [H,W]=size(ref);bestVal=inf;bestXY=center;
    for i=1:size(pattern,1)
        xy=center+pattern(i,:);dx=xy(1);dy=xy(2);
        if abs(dx)>R || abs(dy)>R,continue;end
        ix=dx+R+1;iy=dy+R+1;
        if isnan(memo(iy,ix))
            p=ref(min(max(rows+dy,1),H),min(max(cols+dx,1),W));
            memo(iy,ix)=sum((target(:)-p(:)).^2);
        end
        if memo(iy,ix)<bestVal
            bestVal=memo(iy,ix);bestXY=xy;
        end
    end
end

function [pred,mv]=affineME(cur,ref,bs,mvInit)
    % Four-parameter field: u=a*x+b*y+c; v=-b*x+a*y+f.
    % Gradient linearization / least-squares as discussed by Li et al.
    % Starts from full-search translation and keeps the best solution found.
    [H,W]=size(ref);nr=H/bs;nc=W/bs;
    [gy,gx]=gradient(ref);
    [lx,ly]=meshgrid(0:bs-1,0:bs-1);
    [X,Y]=meshgrid(-(bs-1)/2:(bs-1)/2,-(bs-1)/2:(bs-1)/2);
    pred=zeros(H,W);mv=mvInit;
    for r=1:nr
        rows=(r-1)*bs+(1:bs);
        for c=1:nc
            cols=(c-1)*bs+(1:bs); x=cols(1)+lx;y=rows(1)+ly;
            target=cur(rows,cols);p=[0 0 mvInit.dx(r,c) mvInit.dy(r,c)];
            bestPatch=sampleAffine(ref,x,y,X,Y,p);
            bestSSD=sum((target(:)-bestPatch(:)).^2);
            for it=1:8
                u=p(1)*X+p(2)*Y+p(3);v=-p(2)*X+p(1)*Y+p(4);
                xq=min(max(x+u,1),W);yq=min(max(y+v,1),H);
                patch=sampleBlock(ref,xq,yq);
                Gx=interp2(gx,xq,yq,'linear',0);
                Gy=interp2(gy,xq,yq,'linear',0);
                A=[Gx(:).*X(:)+Gy(:).*Y(:),Gx(:).*Y(:)-Gy(:).*X(:),Gx(:),Gy(:)];
                residual=target(:)-patch(:);
                dp=(A'*A+1e-3*eye(4))\(A'*residual);
                if any(~isfinite(dp)),break;end
                accepted=false;
                for alpha=[1 .5 .25 .125]
                    trial=p+alpha*dp.';
                    % Restrict excessive local distortion/translation.
                    if abs(trial(1))>.5 || abs(trial(2))>.5 || abs(trial(3))>16 || abs(trial(4))>16
                        continue;
                    end
                    q=sampleAffine(ref,x,y,X,Y,trial);
                    e=sum((target(:)-q(:)).^2);
                    if e<bestSSD-1e-8
                        p=trial;bestSSD=e;bestPatch=q;accepted=true;break;
                    end
                end
                if ~accepted,break;end
            end
            pred(rows,cols)=bestPatch;
            mv.dx(r,c)=p(3);mv.dy(r,c)=p(4); % center displacement for quiver
        end
    end
end

function p=sampleAffine(ref,x,y,X,Y,a)
    u=a(1)*X+a(2)*Y+a(3);
    v=-a(2)*X+a(1)*Y+a(4);
    p=sampleBlock(ref,x+u,y+v);
end

function p=sampleBlock(ref,xq,yq)
    [H,W]=size(ref);
    p=interp2(ref,min(max(xq,1),W),min(max(yq,1),H),'linear');
end

function showPanels(cur,pred,mv,bs,name,t,mode,outDir)
    residual=cur-pred;
    [H,W]=size(cur);nr=H/bs;nc=W/bs;
    [xc,yc]=meshgrid(((1:nc)-.5)*bs,((1:nr)-.5)*bs);
    fig=figure('Visible','off','Color','w','Position',[20 30 1250 850]);
    subplot(2,2,1);imagesc(cur,[0 255]);colormap(gray);axis image off;
    title(sprintf('%s frame %d (current)',name,t));
    subplot(2,2,2);imagesc(pred,[0 255]);axis image off;title(sprintf('Prediction (%s, %dx%d)',mode,bs,bs));
    subplot(2,2,3);imagesc(residual);axis image off;colorbar;title(sprintf('Residual, MSE %.5f',mean(residual(:).^2)));
    subplot(2,2,4);imagesc(cur,[0 255]);axis image off;hold on;
    % Decimate ONLY the vector DISPLAY when vectors would be too crowded.
    stride=max(1,ceil(12/bs));ri=1:stride:nr;ci=1:stride:nc;
    quiver(xc(ri,ci),yc(ri,ci),mv.dx(ri,ci),mv.dy(ri,ci),0,'r');
    title('Motion vectors (quiver)');
    dest=fullfile(outDir,sprintf('%s_F%02d_B%d_%s.png',name,t,bs,mode));
    exportgraphics(fig,dest,'Resolution',120);close(fig);
end
