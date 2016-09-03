% =========================================================================
% Title       : Simulator for Quanized Massive MU-MIMO-OFDM Uplink
% File        : consolidate.m
% -------------------------------------------------------------------------
% Description :
%
%   Collects at most "number" simulation results and computes avereage
%   error rates and complexities.  
%
% -------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   03-sep-16  1.0      studer  cleanup and release on github
% -------------------------------------------------------------------------
%   (C) 2016 Christoph Studer (email: studer@cornell.edu)
% =========================================================================

function Out = consolidate(filename,number);
  disp(sprintf('consolidate %s ...',filename));
  k = 0;
  for i=1:number
    fullname = ['results/',filename,'_',num2str(i-1),'.mat'];
    % -- check whether file exists
    if exist(fullname)>0
      tmp = load(fullname);
      if k==0
        ptype = tmp;
        BERunc = ptype.Results.BERunc;
        PERunc = ptype.Results.PERunc;
        BER = ptype.Results.BER;
        PER = ptype.Results.PER;
        NrOfSimPackets = ptype.Results.NrOfSimPackets;
      else
        BERunc = BERunc + tmp.Results.BERunc;
        PERunc = PERunc + tmp.Results.PERunc;
        BER = BER + tmp.Results.BER;
        PER = PER + tmp.Results.PER;
        NrOfSimPackets = NrOfSimPackets + tmp.Results.NrOfSimPackets;
      end
      k=k+1;
    end
  end
  % -- average and compute output
  ptype.Results.BERunc = BERunc/k;
  ptype.Results.PERunc = PERunc/k; 
  ptype.Results.NrOfSimPackets = NrOfSimPackets;
  disp(sprintf('%i packets simulated.\n',NrOfSimPackets));
  Out = ptype.Results;
return
