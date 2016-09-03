% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : init.m
% -------------------------------------------------------------------------
% Description :
%
%   Standard detector framework. Initializes all parameters.
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================             

function [TxRx,Result] = init(TxRx,RunID)

  % -- initialization of random seed
  rng(RunID);
   
  % -- construct IEEE 802.11a-compliant mapping/constellation points
  switch (TxRx.Modulation_order)
    case 1,
      TxRx.Constellations = [ -1 1 ];
    case 2,
      TxRx.Constellations = [ -1-1j,-1+1j, ...
                              +1-1j,+1+1j ];
    case 4,
      TxRx.Constellations = [ -3-3j,-3-1j,-3+3j,-3+1j, ...
                              -1-3j,-1-1j,-1+3j,-1+1j, ...
                              +3-3j,+3-1j,+3+3j,+3+1j, ...
                              +1-3j,+1-1j,+1+3j,+1+1j ];
    case 6,
      TxRx.Constellations = [ -7-7j,-7-5j,-7-1j,-7-3j,-7+7j,-7+5j,-7+1j,-7+3j, ...
                              -5-7j,-5-5j,-5-1j,-5-3j,-5+7j,-5+5j,-5+1j,-5+3j, ...
                              -1-7j,-1-5j,-1-1j,-1-3j,-1+7j,-1+5j,-1+1j,-1+3j, ...
                              -3-7j,-3-5j,-3-1j,-3-3j,-3+7j,-3+5j,-3+1j,-3+3j, ...
                              +7-7j,+7-5j,+7-1j,+7-3j,+7+7j,+7+5j,+7+1j,+7+3j, ...
                              +5-7j,+5-5j,+5-1j,+5-3j,+5+7j,+5+5j,+5+1j,+5+3j, ...
                              +1-7j,+1-5j,+1-1j,+1-3j,+1+7j,+1+5j,+1+1j,+1+3j, ...
                              +3-7j,+3-5j,+3-1j,+3-3j,+3+7j,+3+5j,+3+1j,+3+3j ];
    otherwise,  
      error('Modulation order for 11a mapping not supported.')
  end
  % -- normalize power such that Ex[|s|^2]=1/TxRX.Ntx
  const_power = sum(abs(TxRx.Constellations).^2)/length(TxRx.Constellations);
  TxRx.Constellations_norm = TxRx.Constellations/sqrt(const_power*TxRx.Ntx);
    
  % -- count number of used OFDM tones
  TxRx.Nused = length(TxRx.ToneMap);
  
  % -- generate code and interleaver data
  TxRx.Code.trellis = poly2trellis(TxRx.Code.K, TxRx.Code.generators); % requires communication toolbox
  for nrx=1:TxRx.Ntx
    TxRx.Code.InterleaverPerm(nrx,:) = randperm(TxRx.Nused*TxRx.Modulation_order); % random interleaver
  end

  % -- puncturing according to IEEE 802.11n
  switch (TxRx.Code.Rate) % total code rate (after puncturing)
    case 1/2, % -- no puncturing
      TxRx.Code.Puncturing.Period = 1;
      TxRx.Code.Puncturing.Pattern = [ 1 ];
    case 3/4,
      TxRx.Code.Puncturing.Period = 18;
      TxRx.Code.Puncturing.Pattern = [ 1 2 3 6 7 8 9 12 13 14 15 18 ];
    case 2/3, 
      TxRx.Code.Puncturing.Period = 12;
      TxRx.Code.Puncturing.Pattern = [ 1 2 3 5 6 7 9 10 11];
    case 5/6,
      TxRx.Code.Puncturing.Period = 10;
      TxRx.Code.Puncturing.Pattern = [ 1 2 3 6 7 10 ];
    otherwise,  
      error('Rate not supported')    
  end
  codedBits = length(TxRx.ToneMap)*TxRx.Modulation_order; % bits per terminal
  patLen = length(TxRx.Code.Puncturing.Pattern); 
  for i=1:ceil(codedBits/patLen)       
    TxRx.Code.Puncturing.Index(patLen*(i-1)+1:patLen*i) = ...
          (i-1)*TxRx.Code.Puncturing.Period+TxRx.Code.Puncturing.Pattern;   
  end
  TxRx.Code.Puncturing.Index(:,codedBits+1:end) = [];
  TxRx.Code.PADbits = ceil(TxRx.Ntx*TxRx.Nused*TxRx.Modulation_order*TxRx.Code.Rate)-floor(TxRx.Ntx*TxRx.Nused*TxRx.Modulation_order*TxRx.Code.Rate);
  
  % -- Compute SNR in dB
  TxRx.Sim.SNR_list = 10.^(TxRx.Sim.SNR_dB_list/10);
  
  % -- unused tones
  TxRx.ToneMapC = setdiff(1:TxRx.N,union(TxRx.ToneMap,TxRx.PilotMap));
  
  % -- prepare Result structure
  Result.total_time = 0;
  Result.BERunc = zeros(length(TxRx.Sim.SNR_list),TxRx.Ntx); % bit error rate
  Result.PERunc = zeros(length(TxRx.Sim.SNR_list),TxRx.Ntx); % packet error rate  
  Result.BER = zeros(length(TxRx.Sim.SNR_list),TxRx.Ntx); % bit error rate
  Result.PER = zeros(length(TxRx.Sim.SNR_list),TxRx.Ntx); % packet error rate
    
return
  