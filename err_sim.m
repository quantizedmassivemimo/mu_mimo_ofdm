% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : err_sim.m
% -------------------------------------------------------------------------
% Description :
%
%   Main simulation file for Quantized Massive MU-MIMO-OFDM Uplink as
%   described in the paper:
%
%   [A] C. Studer and G. Durisi, "Quantized Massive MU-MIMO-OFDM Uplink,"
%   IEEE Transactions on Communications, 2016
%
%   * If you use this simulator in any way, then you must cite [A] *
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.3      studer  cleanup and release on github
%   01-dec-08  1.2      studer  major cleanup and speedup
%   14-sep-06  1.1      studer  enhanced with iterative MIMO decoding
%   10-aug-06  1.0      studer  adapted from ofdm_sim environment
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================

function err_sim(TxRx,RunID)

% -- initialization
[TxRx,Result] = init(TxRx,RunID);
K_list = 1:length(TxRx.Sim.SNR_list);

% -- main simulation loop (packets)
tic;
for pack=1:TxRx.Sim.nr_of_packets
    
    % -- generate channel matrices according to specified channel model
    HN = channel(TxRx);
    
    % -- transmitter: generates pilots and bits, maps to symbol vectors,
    % -- and transmits symbol vectors over (noiseless) channel
    [data_bin,data_bin_coded,rxN,tN,sN] = tx(TxRx,HN);
    
    % -- generate noise realization
    noise = sqrt(0.5)*(randn(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N)+1i*randn(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N));
    
    for k=K_list
        
        % -- set noise variance
        N0 = 1/TxRx.Sim.SNR_list(k);
        
        % -- add Gaussian noise (FD channel model)
        yFreq = rxN + noise*sqrt(N0);
        
        % -- TD channel with quantization (time x Nrx x pilot + data symbols)
        yTime = zeros(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N);
        for tpd=1:TxRx.Ntrain+TxRx.Ndata
            yTime(:,tpd,:) = sqrt(TxRx.N)*ifft(ifftshift(squeeze(yFreq(:,tpd,:)),2),[],2);
        end
        
        % -- quantizer (quantize directly to bin position information)
        
        % quantize with equiprobable bins
        w = sort([ real(yTime(:)); imag(yTime(:)) ],'ascend');
        q = w(round(linspace(1,length(w),1+TxRx.Q.Bins)));
        q(1) = -inf; q(end) = +inf;
        TxRx.Q.levels = [q(1:end-1) q(2:end)];
        % find centroid within each bin (used in simpler, less accurate recovery method)
        lsp = round(linspace(1,length(w),1+TxRx.Q.Bins));
        for bb=1:TxRx.Q.Bins
            widx = w(lsp(bb):lsp(bb+1));
            TxRx.Q.centers(bb) = mean(widx); % conditional mean (centroid)
            TxRx.Q.variances(bb) = mean(abs(widx-TxRx.Q.centers(bb)).^2); % variance for real or imaginary part
        end
        
        % quantize!
        yQ.r = zeros(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N,2); % real part
        yQ.i = zeros(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N,2); % imaginary part
        for ii=1:TxRx.Nrx
            for jj=1:TxRx.Ntrain+TxRx.Ndata
                for kk=1:TxRx.N
                    % bin quantizer
                    ridx = sum((real(yTime(ii,jj,kk))>TxRx.Q.levels),2)==1;
                    iidx = sum((imag(yTime(ii,jj,kk))>TxRx.Q.levels),2)==1;
                    yQ.r(ii,jj,kk,:) = TxRx.Q.levels(ridx,:);
                    yQ.i(ii,jj,kk,:) = TxRx.Q.levels(iidx,:);
                    % centroid quantizer
                    yQ.m(ii,jj,kk) = TxRx.Q.centers(ridx) + 1i*TxRx.Q.centers(iidx);
                    yQ.v(ii,jj,kk) = TxRx.Q.variances(ridx) + 1i*TxRx.Q.variances(iidx); % contains variance of real and imaginary part
                end
            end
        end
        
        % -- wrapper for detectors (detector needs to know training and pilots!)
        switch (TxRx.DetectorName)
            case {'Infty','MMSE'} % quantized recovery via convex optimzation
                [LLR_E1,Stats] = detector_quant(TxRx,yQ,tN,sN,N0,HN);
            case 'OneShot' % low-complexity, one-shot detection method
                [LLR_E1,Stats] = detector_oneshot(TxRx,yQ,tN,sN,N0,HN);
            case 'SIMO' % Simulate SIMO bound with perfect CSIR
                [LLR_E1,Stats] = detector_SIMO(TxRx,yFreq,tN,sN,N0,HN);
            case {'REFINFTY','REFMMSE','REFOS'} % Reference curves assume unquantized system
                [LLR_E1,Stats] = detector_ref(TxRx,yTime,tN,sN,N0,HN);
        end
        
        % -- wrapper for channel decoders
        [LLR_P1,data_bin_hat] = decoder(TxRx,LLR_E1);
        
        errs_uncoded = sum(sum(abs((LLR_E1>0)-data_bin_coded),3),2)';
        errs_coded = sum(sum(abs(data_bin_hat-data_bin),3),2)'; % sum(abs(data_bin_hat-data_bin)');
        
        % -- compute coded BER/FER
        %tmp = sum(sum(abs(data_bin_hat-data_bin)))
        %[ complexity , tmp ]
        Result.BERunc(k,:) = Result.BERunc(k,:) + errs_uncoded;
        Result.PERunc(k,:) = Result.PERunc(k,:) + (errs_uncoded>0);
        Result.BER(k,:) = Result.BER(k,:) + errs_coded;
        Result.PER(k,:) = Result.PER(k,:) + (errs_coded>0);
        Result.Stats(k,pack) = Stats;
        
    end % end SNR loop
    
    % -- save intermediate results
    Result.total_time = Result.total_time + toc;
    save_results(TxRx,RunID,Result,pack);
    
    % -- show # of elapsed packets
    disp(sprintf('\nPacket %i of %i transmitted',pack,TxRx.Sim.nr_of_packets))
    
end % end packet loop

% -- simulation completed!
disp(' =====================================================');
disp(sprintf(' ==================== Total Simulation Time: %f sec. ',Result.total_time));
disp(' =====================================================');

% -- store FINAL results to disk
save_results(TxRx,RunID,Result,TxRx.Sim.nr_of_packets);

return

% -------------------------------------------------------------------------
% -- save (intermediate/final) results to disk
% -------------------------------------------------------------------------
function save_results(TxRx,RunID,Res,NrOfSimPackets)

% -- compose simulation results structure
Results.TxRx = TxRx;
Results.RunName = TxRx.Sim.name;
Results.RunID = RunID;
Results.simulationTime = Res.total_time;
Results.NrOfSimPackets = NrOfSimPackets; % store actual number of simulated channels
if (NrOfSimPackets==TxRx.Sim.nr_of_packets)
    Results.Final = 'yes';
else
    Results.Final = 'no';
end
% -- for each iteration
Results.BERunc = Res.BERunc./(Results.NrOfSimPackets*TxRx.Ndata*TxRx.Nused*TxRx.Modulation_order);
Results.PERunc = Res.PERunc./Results.NrOfSimPackets;
Results.BER = Res.BER./(Results.NrOfSimPackets*TxRx.Ndata*TxRx.Nused*TxRx.Modulation_order*TxRx.Code.Rate);
Results.PER = Res.PER./Results.NrOfSimPackets;
Results.Stats = Res.Stats;

% -- store Results structure to a separate file
Results.FileName = sprintf('results/%s_%d.mat',TxRx.Sim.name,RunID);
save(Results.FileName,'Results');

return