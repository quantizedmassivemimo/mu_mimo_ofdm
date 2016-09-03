% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : detector_ref.m
% -------------------------------------------------------------------------
% Description :
%
%   Compute reference curves for unquantized detection.
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================

function [LLR_E1,Stats] = detector_ref(TxRx,yTime,tN,sN,N0,HN)

% -- initialization
bin_array = sign(de2bi([0:2^(TxRx.Modulation_order)-1],TxRx.Modulation_order)-0.5);

%% == channel estimation (if enabled)
switch (TxRx.CSIR)
    case 'no', %  no CSIR -> channel estimation
        
        % -- use FASTA to perform channel estimation
        %opts.stopRule='normalizedResidual';
        opts = [];
        opts.stopRule='ratioResidual';
        opts.tol = 1e-6;  % Use strict tolerance
        opts.recordObjective = false; %  Record the objective function so we can plot it
        
        % -- set up functions
        A = @(x) train_A(TxRx,tN,x);
        At = @(y) train_At(TxRx,tN,y);
        gradf = @(y) train_grad(TxRx,y,yTime,N0); %%%%%%
        proxg = @(x,t) train_proxy(TxRx,t,x);
        f = @(y) train_f(TxRx,y,yTime,N0); %%%%%%\
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
        
        % -- corrected gain loss
        MSEdBC = 10*log10(norm(HN(:)-norm(HN(:),'fro')*HEst(:)/norm(HEst(:),'fro'),'fro')^2/(TxRx.Nrx*TxRx.Ntx*TxRx.N));
        Stats.MSEdBC = MSEdBC;
        
    case 'yes', %  perfect CSIR -> channel estimation
        HEst = HN;
        Stats = NaN;
end

%% == data detection

% -- use FASTA to perform channel estimation
opts = [];
opts.stopRule='ratioResidual';
opts.maxIters = 500;
%opts.verbose = true;
opts.tol = 1e-3;  % Use strict tolerance
%opts.recordObjective = true; %  Record the objective function so we can plot it

% -- set up functions
A = @(x) detect_A(TxRx,HEst,sN,x);
At = @(y) detect_At(TxRx,HEst,sN,y);
Es = 1/TxRx.Ntx;

gradf = @(y) detect_grad(TxRx,y,yTime,N0); %%%%%
f = @(y) detect_f(TxRx,y,yTime,N0); %%%%%

switch (TxRx.DetectorName)
    case 'REFMMSE'
        proxg = @(x,t) x/(1+t/Es); %%% Gaussian prior (MMSE equivalent)
        g = @(x)  norm(x,'fro')^2/Es; % Gaussian prior (MMSE equivalent)
    case 'REFINFTY' % use box constraint
        proxg = @(x,t) detect_proxg(TxRx,t,x) ; % apply box constraint
        g = @(x) 0; %%% detect_g(TxRx,x);
end

% == call FASTA solver
x0 = zeros(TxRx.Ntx*TxRx.Ndata*TxRx.N,1);
[solution, outs] = fasta(A,At,f,gradf,g,proxg,x0,opts);
Stats.DETECTtime = outs.solvetime;

xEst = reshape(solution,TxRx.Ntx,TxRx.Ndata,TxRx.N);
for nda=1:TxRx.Ndata
    for nto=1:length(TxRx.ToneMap); % only estimate for data vectors
        dist = abs(ones(TxRx.Ntx,1)*TxRx.Constellations_norm-xEst(:,nda,TxRx.ToneMap(nto))*ones(1,2^TxRx.Modulation_order)).^2;
        switch (TxRx.DetectorSoft)
            case 'Yes'
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

%% == FUNCTIONS FOR CHANNEL ESTIMATION

% -- compute cost function
function cost = train_f(TxRx,z,yTime,N0)
z = z(1:TxRx.Nrx*TxRx.Ntrain*TxRx.N); % only extract data from used tones
yUsed = reshape(yTime(:,1:TxRx.Ntrain,:,:),TxRx.Nrx*TxRx.Ntrain*TxRx.N,1);
cost = 0.5*norm(z-yUsed,2)^2/N0;
end

% -- compute proxy for finite tap channel
function xprox = train_proxy(TxRx,t,xin)
HEst = reshape(xin,TxRx.Nrx,TxRx.Ntx,TxRx.N);
% -- apply projection for each spatial stream
for kk = 1:TxRx.Nrx
    for ll = 1:TxRx.Ntx
        timedomain =  ifft(ifft1shift(squeeze(HEst(kk,ll,:)))); % go to time domain
        timedomain(TxRx.Channel.MaxNtaps+1:end) = 0; % kill inexistent taps
        HEst(kk,ll,:) = fft1shift(fft(timedomain)); % go back to frequency domain
    end
end
xprox = reshape(HEst,TxRx.Nrx*TxRx.Ntx*TxRx.N,1);
end

% -- compute gradient for channel training
function ygrad=train_grad(TxRx,y,yTime,N0)

% two gradients p(Y|HX) and p(X)

yUsed = reshape(yTime(:,1:TxRx.Ntrain,:,:),TxRx.Nrx*TxRx.Ntrain*TxRx.N,1);
% only take gradient on time-domain component
z = y(1:TxRx.Nrx*TxRx.Ntrain*TxRx.N);
ygrad = [ (z-yUsed)/N0 ; ... %%%%%
    y(TxRx.Nrx*TxRx.Ntrain*TxRx.N+1:end) ];
end


% -- forward transform A(H) for channel training
function y = train_A(TxRx,tN,x)
% apply training to channel matrices from right
HNest = reshape(x,TxRx.Nrx,TxRx.Ntx,TxRx.N);
rHN = zeros(TxRx.Nrx,TxRx.Ntrain,TxRx.N);
for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
    rHN(:,:,nto) = HNest(:,:,nto)*tN(:,:,nto);
end
% conversion from frequency to time domain
rT = zeros(TxRx.Nrx,TxRx.Ntrain,TxRx.N);
for tpd=1:TxRx.Ntrain
    rT(:,tpd,:) = sqrt(TxRx.N)*ifft(ifftshift(squeeze(rHN(:,tpd,:)),2),[],2);
end
% keep the H matrices on the unused tones and don't do anything with them
unusedSet = setdiff(1:TxRx.N,union(TxRx.ToneMap,TxRx.PilotMap));
unusedH = HNest(:,:,unusedSet);
% output time-domain signals and unused parts too
y = [ rT(:); unusedH(:) ];
end

% -- adjoint transform At(Y) for channel training
function x = train_At(TxRx,tN,y)
rT = reshape(y(1:TxRx.Nrx*TxRx.Ntrain*TxRx.N),TxRx.Nrx,TxRx.Ntrain,TxRx.N);
% conversion from time to frequency domain
rF = zeros(TxRx.Nrx,TxRx.Ntrain,TxRx.N);
for tpd=1:TxRx.Ntrain
    rF(:,tpd,:) = 1/sqrt(TxRx.N)*fftshift(fft(squeeze(rT(:,tpd,:)),[],2),2);
end
% apply adjoint to channel from right
xHNest = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
    xHNest(:,:,nto) = rF(:,:,nto)*tN(:,:,nto)';
end
% keep the H matrices on the unused tones and don't do anything with them
unusedSet = setdiff(1:TxRx.N,union(TxRx.ToneMap,TxRx.PilotMap));
xHNest(:,:,unusedSet) = reshape(y(TxRx.Nrx*TxRx.Ntrain*TxRx.N+1:end),TxRx.Nrx,TxRx.Ntx,length(unusedSet));
% output vector
x = xHNest(:);
end

%%%%
%%%%
%%%%
%%%%

% -- forward transform H*X for detection (sN needed for pilot tones)
function y = detect_A(TxRx,HNest,sN,x)
% transmit through channels
xN = reshape(x,TxRx.Ntx,TxRx.Ndata,TxRx.N);
rN = zeros(TxRx.Nrx,TxRx.Ndata,TxRx.N);
% pilots are known, so just replace them with the correct data
for npo=TxRx.PilotMap
    rN(:,:,npo) = HNest(:,:,npo)*sN(:,:,npo); % known pilots
end
for nto=TxRx.ToneMap;
    rN(:,:,nto) = HNest(:,:,nto)*xN(:,:,nto); % unknown data
end
% conversion from frequency to time domain
rT = zeros(TxRx.Nrx,TxRx.Ndata,TxRx.N);
for tpd=1:TxRx.Ndata
    rT(:,tpd,:) = sqrt(TxRx.N)*ifft(ifftshift(squeeze(rN(:,tpd,:)),2),[],2);
end
% output time-domain signals and unused parts too
y = rT(:);
end

% -- adjoint transform H'*Y for detection (sN needed for pilot tones)
function x = detect_At(TxRx,HNest,sN,y)
rT = reshape(y,TxRx.Nrx,TxRx.Ndata,TxRx.N);
% conversion from time to frequency domain
rN = zeros(TxRx.Nrx,TxRx.Ndata,TxRx.N);
for tpd=1:TxRx.Ndata
    rN(:,tpd,:) = 1/sqrt(TxRx.N)*fftshift(fft(squeeze(rT(:,tpd,:)),[],2),2);
end
% apply adjoint to received data from left
xNest = zeros(TxRx.Ntx,TxRx.Ndata,TxRx.N);
for nto=TxRx.ToneMap; % only estimate for data vectors
    xNest(:,:,nto) = HNest(:,:,nto)'*rN(:,:,nto);
end
% we know what has been sent on the pilots!
for npo=TxRx.PilotMap % replace known pilots
    xNest(:,:,npo) = sN(:,:,npo);
end
x = xNest(:);
end

% -- compute gradient for data detection
function ygrad=detect_grad(TxRx,y,yTime,N0)
% vectorize quantization arrays
yUsed = reshape(yTime(:,TxRx.Ntrain+1:end,:,:),TxRx.Nrx*TxRx.Ndata*TxRx.N,1);
% only take gradient on time-domain component
ygrad = (y-yUsed)/N0; %%%%%%
end

% -- compute cost function
function cost = detect_f(TxRx,z,yTime,N0)
yUsed = reshape(yTime(:,TxRx.Ntrain+1:end,:,:),TxRx.Nrx*TxRx.Ndata*TxRx.N,1);
cost = 0.5*norm(z-yUsed,2)^2/N0;
end

% -- compute proxy for constellation linfty-tilde norm
function xproxg = detect_proxg(TxRx,t,x)

constmax = max(real(TxRx.Constellations_norm));

xr = real(x);
ridx = abs(xr)>constmax;
xr(ridx) = constmax*sign(xr(ridx));

xi = imag(x);
iidx = abs(xi)>constmax;
xi(iidx) = constmax*sign(xi(iidx));

xproxg = xr + 1i*xi;

end

%%%%
%%%%
%%%%
%%%%

function out = fft1shift(in)
m = length(in);
p = ceil(m/2);
out = [ in(p+1:m) ; in(1:p) ];
end

function out = ifft1shift(in)
m = length(in);
p = floor(m/2);
out = [ in(p+1:m) ; in(1:p) ];
end


%%%%
%%%%
%%%%
%%%%



