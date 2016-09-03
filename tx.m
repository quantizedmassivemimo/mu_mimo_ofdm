% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : tx.m
% -------------------------------------------------------------------------
% Description :
%
%   This file generates random bits, applies forward error-correction with
%   a convolutional encoder (note: convenc reqiures communication toolbox).
%   Furthermore, interleaving and modulation is performed and the symbol
%   vectors are transmitted over the channel matrices (without noise).
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% ========================================================================= 

function [data_bin,data_bin_coded_stream,rxN,tN,sN] = tx(TxRx,HN)

  % -- generate random bits and apply forward error correction   
  data_bin = randi([0 1],TxRx.Ntx,ceil(TxRx.Nused*TxRx.Modulation_order*TxRx.Code.Rate),TxRx.Ndata);
  data_bin(1:TxRx.Ntx,[size(data_bin,2)+1-(TxRx.Code.K-1)-TxRx.Code.PADbits:size(data_bin,2)],:) = 0;
  
  % -- separate enconding for each terminal
  for nda=1:TxRx.Ndata
    for ntx=1:TxRx.Ntx
      data_bin_coded = convenc(data_bin(ntx,:,nda),TxRx.Code.trellis); % requires communication toolbox  
      % -- puncturing
      data_bin_punctured = data_bin_coded(TxRx.Code.Puncturing.Index);  
      % -- interleaving
      data_bin_coded_stream(ntx,:,nda) = data_bin_punctured(TxRx.Code.InterleaverPerm(ntx,:));
    end
  end
  
  % -- modulation: training pilots (Hadamard)
  H = sqrt(1/TxRx.Ntx)*hadamard(2^ceil(log2(max(TxRx.Ntrain,TxRx.Ntx))));
  tN = zeros(TxRx.Ntx,TxRx.Ntrain,TxRx.N);
  for nto=union(TxRx.ToneMap,TxRx.PilotMap); % only train used tones    
    tN(:,:,nto) = diag(sign(randn(TxRx.Ntx,1)))*H*diag(sign(randn(TxRx.Ntx,1))) ;
  end
    
  % -- modulation: data symbol vectors
  sN = zeros(TxRx.Ntx,TxRx.Ndata,TxRx.N);
  for nda=1:TxRx.Ndata
    for ntx=1:TxRx.Ntx
      tmp = reshape(data_bin_coded_stream(ntx,:,nda),TxRx.Modulation_order,TxRx.Nused);
      tmp_dec = bi2de(tmp.').';
      sN(ntx,nda,TxRx.ToneMap) = TxRx.Constellations_norm(tmp_dec+1);
    end
  end
  % and dont forget the pilot tones (fill with random stuff)
  for npi=TxRx.PilotMap
    sN(:,:,npi) = sqrt(1/TxRx.Ntx)*(2*randi([0 1],TxRx.Ntx,TxRx.Ndata)-1);
  end
  
  % -- transmit over the noiseless channel (in frequency domain)
  rxN = zeros(TxRx.Nrx,TxRx.Ntrain+TxRx.Ndata,TxRx.N);
  for n=1:TxRx.N
    rxN(:,:,n) = HN(:,:,n)*[ tN(:,:,n) , sN(:,:,n)];
  end 
      
end

  