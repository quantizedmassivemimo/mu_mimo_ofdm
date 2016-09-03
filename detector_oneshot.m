% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : detector_oneshot.m
% -------------------------------------------------------------------------
% Description :
%
%   Simple one-shot detection method.
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================

function [LLR_E1,Stats] = detector_oneshot(TxRx,yQ,tN,sN,N0,HN)

% -- initialization
bin_array = sign(de2bi([0:2^(TxRx.Modulation_order)-1],TxRx.Modulation_order)-0.5);

%% == convert centroids and variances to frequency domain

yF = time2freq(TxRx,yQ.m);
vF = var2freq(TxRx,yQ.v);

%% == channel estimation (if enabled)
switch (TxRx.CSIR)
    case 'no', %  no CSIR -> channel estimation
        tic
        HEst1 = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
        for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
            HEst1(:,:,nto) = yF(:,1:TxRx.Ntrain,nto)*tN(:,:,nto)';
        end
        
        opts = [];
        opts.stopRule='ratioResidual';
        opts.tol = 1e-6;  % Use strict tolerance
        opts.recordObjective = false; %  Record the objective function so we can plot it
        
        % -- perform interpolation per channel entry
        Stats.CHESTtime=0;
        HEst = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
        solvetime = zeros(TxRx.Nrx,TxRx.Ntx);
        for ll = 1:TxRx.Ntx
            parfor kk = 1:TxRx.Nrx
                hest = squeeze(HEst1(kk,ll,:));
                v = squeeze(vF(kk,ll,:));
                
                A = @(x) x;
                At = @(y) y;
                gradf = @(x) vec_chest_gradf(TxRx,x,hest,v);
                proxg = @(x,t) vec_channel_proxy(TxRx,x);
                f = @(x) vec_chest_f(TxRx,x,hest,v);
                g = @(x) 0;
                
                x0 = zeros(TxRx.N,1);
                [solution, outs] = fasta(A,At,f,gradf,g,proxg,x0,opts);
                solvetime(kk,ll) = outs.solvetime;
                HEst(kk,ll,:) = solution;
            end
        end
        
        % -- Extract some statistics
        Stats.CHESTtime = solvetime(:);
        MSE = norm(HN(:)-HEst(:),'fro')^2/(TxRx.Nrx*TxRx.Ntx*TxRx.N);
        MSEdB = 10*log10(MSE);
        Stats.MSEdB = MSEdB;
        
        % corrected gain loss
        MSEdBC = 10*log10(norm(HN(:)-norm(HN(:),'fro')*HEst(:)/norm(HEst(:),'fro'),'fro')^2/(TxRx.Nrx*TxRx.Ntx*TxRx.N));
        Stats.MSEdBC = MSEdBC;
        
        % -- channel rescaling
        switch (TxRx.ChannelRescaling)
            case 'None'
            case 'Variance'
                HEst = TxRx.Nrx*TxRx.Ntx*HEst/norm(HEst(:),'fro');
            case 'Oracle'
                HEst = norm(HN(:),'fro')*HEst/norm(HEst(:),'fro');
        end
        
    case 'legacy', %  no CSIR -> use legacy joint channel estimation
        tic
        HEst = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
        for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
            HEst(:,:,nto) = yF(:,1:TxRx.Ntrain,nto)*tN(:,:,nto)';
        end
        
        opts = [];
        opts.stopRule='ratioResidual';
        opts.tol = 1e-6;  % Use strict tolerance
        opts.recordObjective = false; %  Record the objective function so we can plot it
        
        A = @(X) X;
        At = @(Y) Y;
        gradf = @(X) chest_gradf(TxRx,X,HEst,vF);
        proxg = @(X,t) channel_proxy(TxRx,X);
        f = @(X) chest_f(TxRx,X,HEst,vF);
        g = @(x) 0;
        
        % -- call FASTA solver
        x0 = zeros(TxRx.Nrx*TxRx.Ntx*TxRx.N,1);
        [solution, outs] = fasta(A,At,f,gradf,g,proxg,x0,opts);
        Stats.CHESTtime = outs.solvetime;
        
        % -- Extract some statistics
        HEst = reshape(solution,TxRx.Nrx,TxRx.Ntx,TxRx.N);
        MSE = norm(HN(:)-HEst(:),'fro')^2/(TxRx.Nrx*TxRx.Ntx*TxRx.N);
        MSEdB = 10*log10(MSE);
        Stats.MSEdB = MSEdB;
        
        % corrected gain loss
        MSEdBC = 10*log10(norm(HN(:)-norm(HN(:),'fro')*HEst(:)/norm(HEst(:),'fro'),'fro')^2/(TxRx.Nrx*TxRx.Ntx*TxRx.N));
        Stats.MSEdBC = MSEdBC;
        
        % -- channel rescaling
        switch (TxRx.ChannelRescaling)
            case 'None'
            case 'Variance'
                HEst = TxRx.Nrx*TxRx.Ntx*HEst/norm(HEst(:),'fro');
            case 'Oracle'
                HEst = norm(HN(:),'fro')*HEst/norm(HEst(:),'fro');
        end
        
    case 'simple', %  no CSIR -> channel estimation (simple scheme)
        tic
        HEst = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
        for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
            HEst(:,:,nto) = yF(:,1:TxRx.Ntrain,nto)*tN(:,:,nto)';
        end
        Stats.CHESTtime = toc;
        
        
        
end

%%% -- improve channel estimation (de-noising)

% -- Extract some statistics
MSE = norm(HN(:)-HEst(:),'fro')^2/(TxRx.Nrx*TxRx.Ntx*TxRx.N);
MSEdB = 10*log10(MSE);
Stats.MSEdB = MSEdB;

% -- corrected gain loss
MSEdBC = 10*log10(norm(HN(:)-norm(HN(:),'fro')*HEst(:)/norm(HEst(:),'fro'),'fro')^2/(TxRx.Nrx*TxRx.Ntx*TxRx.N));
Stats.MSEdBC = MSEdBC;

% -- channel rescaling (if desired)
switch (TxRx.ChannelRescaling)
    case 'None'
    case 'Variance'
        HEst = TxRx.Nrx*TxRx.Ntx*HEst/norm(HEst(:),'fro');
    case 'Oracle'
        HEst = norm(HN(:),'fro')*HEst/norm(HEst(:),'fro');
end


%% == Unbiased soft-output MMSE data detection (takes into account quantization)
% this method has less-fine variance information (only obtained in the
% frequency domain. hence it performs slightly worse than the full-
% fledged MMSE detector.
for nda=1:TxRx.Ndata
    for nto=1:length(TxRx.ToneMap); % only estimate for data vectors
        H = HEst(:,:,TxRx.ToneMap(nto));
        D = diag(N0+vF(:,TxRx.Ntrain+nda,TxRx.ToneMap(nto))); % noise variance (including quantization stuff)
        G = mrdivide(H',H*H'+TxRx.Ntx*D*eye(TxRx.Nrx));
        xhat = G*yF(:,TxRx.Ntrain+nda,TxRx.ToneMap(nto));
        M = G*H;
        S = diag(diag(M));
        switch (TxRx.DetectorSoft)
            case 'AP' % ignore SINR scaling
                NPI = ones(TxRx.Ntx,2^TxRx.Modulation_order);
            otherwise
                NPI = real(diag(1/TxRx.Ntx*(M-S)*(M-S)'+G*D*G')) * ones(1,2^TxRx.Modulation_order); %%%%
        end
        dist = abs(ones(TxRx.Ntx,1)*TxRx.Constellations_norm-xhat*ones(1,2^TxRx.Modulation_order)).^2./NPI;
        switch (TxRx.DetectorSoft)
            case {'Yes','AP'}
                for bb=1:TxRx.Modulation_order
                    Q(:,bb) = min(dist(:,bin_array(:,bb)==-1),[],2)-min(dist(:,bin_array(:,bb)==1),[],2);
                end
            case 'No'
                [tmp,idx] = min(dist,[],2);
                Q = bin_array(idx,:);
        end
        LLR_E1(:,(TxRx.Modulation_order*(nto-1)+1):(TxRx.Modulation_order*nto),nda) = Q;
    end
end

end


%%%%
%%%%
%%%%
%%%%

% -- compute proxy for finite tap channel
function y = vec_channel_proxy(TxRx,x)
timedomain =  ifft(ifftshift(x)); % go to time domain
timedomain(TxRx.Channel.MaxNtaps+1:end) = 0; % kill inexistent taps
y = fftshift(fft(timedomain)); % go back to frequency domain
end


% -- compute gradient in frequency comain
function y = vec_chest_gradf(TxRx,x,hest,v)
y = zeros(TxRx.N,1);
idx = union(TxRx.ToneMap,TxRx.PilotMap);
y(idx) = (x(idx)-hest(idx))./v(idx);
end

% -- compute gradient in frequency comain
function f = vec_chest_f(TxRx,x,hest,v)
idx = union(TxRx.ToneMap,TxRx.PilotMap);
f = sum(abs(x(idx)-hest(idx)).^2./v(idx));
end


%%%%
%%%%
%%%%
%%%%


% -- compute proxy for finite tap channel
function Y = channel_proxy(TxRx,X)

% -- convert into matrix form
HN = reshape(X,TxRx.Nrx,TxRx.Ntx,TxRx.N);

% -- apply projection for each spatial stream
for kk = 1:TxRx.Nrx
    for ll = 1:TxRx.Ntx
        timedomain =  ifft(ifftshift(squeeze(HN(kk,ll,:)))); % go to time domain
        timedomain(TxRx.Channel.MaxNtaps+1:end) = 0; % kill inexistent taps
        HN(kk,ll,:) = fftshift(fft(timedomain)); % go back to frequency domain
    end
end

% -- vectorize
Y = HN(:);

end


% -- compute gradient in frequency comain
function Y = chest_gradf(TxRx,X,HEst,vF)
% -- convert into matrix form
XN = reshape(X,TxRx.Nrx,TxRx.Ntx,TxRx.N);
% -- gradient
GN = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
    GN(:,:,nto) = (XN(:,:,nto)-HEst(:,:,nto))./vF(:,1:TxRx.Ntrain,nto);
end

% -- vectorize
Y = GN(:);
end

% -- compute gradient in frequency comain
function v = chest_f(TxRx,X,HEst,vF)
% -- convert into matrix form
XN = reshape(X,TxRx.Nrx,TxRx.Ntx,TxRx.N);
% -- compute cost function f
v = 0;
for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
    v = v + sum(sum(abs(XN(:,:,nto)-HEst(:,:,nto)).^2./vF(:,1:TxRx.Ntrain,nto)));
end
end


%%%%
%%%%
%%%%
%%%%


% conversion from frequency to time domain
function yT = freq2time(TxRx,yF)
% -- from (Nrx x pilot + data symbols x Tones)
% -- to (time x Nrx x pilot + data symbols)
yT = zeros(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N);
for tpd=1:TxRx.Ntrain+TxRx.Ndata
    yT(:,tpd,:) = sqrt(TxRx.N)*ifft(ifftshift(squeeze(yF(:,tpd,:)),2),[],2);
end
end


% conversion from time domain to frequency domain
function yF = time2freq(TxRx,yT)
% -- from (time x Nrx x pilot + data symbols)
% -- to (Nrx x pilot + data symbols x Tones)
yF = zeros(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N);
for tpd=1:TxRx.Ntrain+TxRx.Ndata
    yF(:,tpd,:) = 1/sqrt(TxRx.N)*fftshift(fft(squeeze(yT(:,tpd,:)),[],2),2);
end
end

% conversion from error variances in time domain to variances in frequency
% domain
function vF = var2freq(TxRx,vT)
% -- from (time x Nrx x pilot + data symbols)
% -- to (Nrx x pilot + data symbols x Tones)
vF = zeros(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N);
for tpd=1:TxRx.Ntrain+TxRx.Ndata
    vari = squeeze(vT(:,tpd,:));
    variF = sum(real(vari)+imag(vari),2)/TxRx.N;
    vF(:,tpd,:) = variF*ones(1,TxRx.N);
end
end

%%%%
%%%%
%%%%
%%%%


