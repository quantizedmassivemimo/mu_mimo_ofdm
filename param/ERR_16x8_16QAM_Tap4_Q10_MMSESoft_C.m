function ERR_16x8_16QAM_Tap4_Q10_MMSESoft_C(RunID)

  % -- simulation setup
  TxRx.Sim.name = mfilename;     % Get script name
  TxRx.Sim.MaxPER = 1000;        % Maximum required PER 
  TxRx.Sim.nr_of_packets = 1000; % Number of packets
  TxRx.Sim.SNR_dB_list = [14:2:16];    
  
  % -- system configuration
  TxRx.Modulation_order = 4; % Modulation scheme: BPSK (1), QPSK (2), 16-QAM (4), 64-QAM (6)
  TxRx.Nrx = 16;             % Number of receivers (base station antennas)
  TxRx.Ntx = 8;              % Number of user antennas (=terminals)
  TxRx.N = 128;              % Nr. of subcarriers
  TxRx.Ntrain = TxRx.Ntx;    % train everything once
  TxRx.Ndata = 2;            % data symbols (the more the better)
  TxRx.ToneMap = [-58:-54,-52:-26,-24:-12,-10:-2,2:10,12:24,26:52,54:58]+64; % 40MHz IEEE 802.1n   
  TxRx.PilotMap = [-53,-25,-11,11,25,53]+64; % 40MHz IEEE 802.1n
    
  % -- code and decoder properties
  TxRx.Code.K = 7;                   % Constraint Length
  TxRx.Code.generators = [133  171]; % Generator Polynomial      
  TxRx.Code.Rate = 5/6;              % code rates '1/2','3/4','2/3','5/6'
  
  % -- channel model
  TxRx.Channel.Model = 'Tap'; % 'Tap' (uniform profile)
  TxRx.Channel.Ntaps = 4;
  TxRx.Channel.MaxNtaps = 16;
  
  % -- detector type
  TxRx.DetectorName = 'MMSE'; % 'Infty', 'MMSE'
  TxRx.DetectorSoft = 'Yes'; % 'Yes', 'No
  TxRx.ChannelRescaling = 'None'; % 'None', 'Variance', 'Oracle'
  TxRx.CSIR = 'no'; % 'yes','no','simple'
  
  % -- quantizer
  TxRx.Q.Type = 'centroid'; % 'bins' or 'centroid'
  TxRx.Q.Bins = 2^10; % number of bins

  % -- start simulation
  err_sim(TxRx,RunID);

return
  
