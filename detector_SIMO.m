% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : detector_SIMO.m
% -------------------------------------------------------------------------
% Description :
%
%   Compute SIMO lower bound curves for (un)quantized detection.
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================

function [LLR_E1,Stats] = detector_SIMO(TxRx,yFreq,tN,sN,N0,HN);

% -- initialization
bin_array = sign(de2bi([0:2^(TxRx.Modulation_order)-1],TxRx.Modulation_order)-0.5);

%% == channel estimation (if enabled)
switch (TxRx.CSIR)
    case 'no', %  no CSIR -> channel estimation
        tic
        HEst = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
        for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
            HEst(:,:,nto) = yFreq(:,1:TxRx.Ntrain,nto)*tN(:,:,nto)';
        end
        
        opts = [];
        opts.stopRule='ratioResidual';
        opts.tol = 1e-6;  % Use strict tolerance
        opts.recordObjective = false; %  Record the objective function so we can plot it
        
        A = @(X) X;
        At = @(Y) Y;
        gradf = @(X) chest_gradf(TxRx,X,HEst);
        proxg = @(X,t) channel_proxy(TxRx,X);
        f = @(X) chest_f(TxRx,X,HEst);
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
        
    case 'yes', % perfect CSIR -> no channel estimation
        HEst = HN;
        Stats = NaN;
end




%% == SIMO bound detector

for nda=1:TxRx.Ndata
    for nto=1:length(TxRx.ToneMap); % only estimate for data vectors
        H = HEst(:,:,TxRx.ToneMap(nto));
        sTrue = sN(:,nda,TxRx.ToneMap(nto));
        y = yFreq(:,TxRx.Ntrain+nda,TxRx.ToneMap(nto));
        y_tilde = y-H*sTrue;
        
        % -- SIMO detection main loop
        for nn=[TxRx.Ntx:-1:1]
            
            % -- interference cancellation (with known data)
            y_SIMO = y_tilde + H(:,nn)*sTrue(nn);
            
            % -- compute noise-plus-interference terms
            w = H(:,nn);
            g = w'*w;
            NPI = norm(w,2).^2*N0;
            dist = abs(w'*y_SIMO-g*TxRx.Constellations_norm).^2/NPI;
            
            switch (TxRx.DetectorSoft)
                case 'Yes'
                    for bb=1:TxRx.Modulation_order
                        Q(nn,bb) = min(dist(:,bin_array(:,bb)==-1),[],2)-min(dist(:,bin_array(:,bb)==1),[],2);
                    end
                case 'No'
                    [tmp,idx] = min(dist,[],2);
                    Q(nn,:) = bin_array(idx,:);
            end
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
function Y = chest_gradf(TxRx,X,HEst)
% -- convert into matrix form
XN = reshape(X,TxRx.Nrx,TxRx.Ntx,TxRx.N);
% -- gradient
GN = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
    GN(:,:,nto) = XN(:,:,nto)-HEst(:,:,nto);
end
% -- vectorize
Y = GN(:);
end

% -- compute gradient in frequency comain
function v = chest_f(TxRx,X,HEst)
% -- convert into matrix form
XN = reshape(X,TxRx.Nrx,TxRx.Ntx,TxRx.N);
% -- compute cost function f
v = 0;
for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only estimate for trained tones
    v = v + norm(XN(:,:,nto)-HEst(:,:,nto),'fro')^2;
end
end

