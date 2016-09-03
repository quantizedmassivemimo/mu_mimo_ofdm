% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : decoder.m
% -------------------------------------------------------------------------
% Description :
%
%   Performs channel decoding using a max-log BCJR decoder.
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================

function [LLR_P1,binary_data_hat] = decoder(TxRx,LLR_E1)

  for nda=1:TxRx.Ndata
    for ntx=1:TxRx.Ntx

      % -- de-interleaving
      LLR_A2(1,TxRx.Code.InterleaverPerm(ntx,:)) = LLR_E1(ntx,:,nda);   
          
      % -- de-puncturing
      LLR_A2_punctured = zeros(1,2*length(LLR_A2(1,:))*TxRx.Code.Rate);    
      LLR_A2_punctured(TxRx.Code.Puncturing.Index) = LLR_A2;

      [LLR_P1(ntx,:,nda),binary_data_hat(ntx,:,nda)] = BCJR(TxRx.Code.trellis,LLR_A2_punctured); 
  
    end
  end
  
return
