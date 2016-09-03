% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : channel.m
% -------------------------------------------------------------------------
% Description :
%
%   Wrapper for channel models. This file only includes a simple n-Tap
%   channel model. The TGn channel models have been omitted due to 
%   copyright issues. 
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================

function [HN]=channel(TxRx)

  HN = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N);
    
  switch(TxRx.Channel.Model)
    case 'TGn', % -- TGN channel models
      error('TGn models not supported due to Copyright issues.')
    case 'Tap', % -- tap channel model      
      % -- generate channels in time domain
      Ht = zeros(TxRx.Nrx,TxRx.Ntx,TxRx.N); % Mt x Mr x time
      Ht(:,:,1:TxRx.Channel.Ntaps) = sqrt(0.5/TxRx.Channel.Ntaps)*(randn(TxRx.Nrx,TxRx.Ntx,TxRx.Channel.Ntaps) ...
                                     + 1i*randn(TxRx.Nrx,TxRx.Ntx,TxRx.Channel.Ntaps) );

      % -- convert to frequency domain                                        
      for kk = 1:TxRx.Nrx
        for ll = 1:TxRx.Ntx
          HN(kk,ll,:) = fftshift(fft(squeeze(Ht(kk,ll,:))));
        end
      end       
      
    otherwise,
      error('No valid TxRx.Channel.Model in channel.m')
  end 
  
end
