% =========================================================================
% Title       : Plots for Quanized Massive MU-MIMO-OFDM Uplink
% File        : err_plot.m
% -------------------------------------------------------------------------
% Description :
%
%   Plots the key figures for the Quantized Massive MU-MIMO-OFDM Uplink 
%   as described in the paper:
%
%   [A] C. Studer and G. Durisi, "Quanized Massive MU-MIMO-OFDM Uplink,"
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

function err_plot

clc

% =========================================================================

TargetPER = 1e-2; % set operating point to 1% PER

% -- collect max. 10 files
number = 10;

% -- load simulation data and consolidate (*** does not consolidate stats ***)
SIMOCSIR = consolidate('ERR_16x8_16QAM_Tap4_SIMO_CSIR',number);%%%%%%%
SIMO = consolidate('ERR_16x8_16QAM_Tap4_SIMO',number);%%%%%%%
REFMMSE = consolidate('ERR_16x8_16QAM_Tap4_REFMMSE',number);%%%%%%%

M3 = consolidate('ERR_16x8_16QAM_Tap4_Q3_MMSESoft',number);
M4 = consolidate('ERR_16x8_16QAM_Tap4_Q4_MMSESoft',number);
M5 = consolidate('ERR_16x8_16QAM_Tap4_Q5_MMSESoft',number);
M6 = consolidate('ERR_16x8_16QAM_Tap4_Q6_MMSESoft',number);
M7 = consolidate('ERR_16x8_16QAM_Tap4_Q7_MMSESoft',number);
M8 = consolidate('ERR_16x8_16QAM_Tap4_Q8_MMSESoft',number);
M9 = consolidate('ERR_16x8_16QAM_Tap4_Q9_MMSESoft',number);
M10 = consolidate('ERR_16x8_16QAM_Tap4_Q10_MMSESoft',number);

M3C = consolidate('ERR_16x8_16QAM_Tap4_Q3_MMSESoft_C',number);
M4C = consolidate('ERR_16x8_16QAM_Tap4_Q4_MMSESoft_C',number);
M5C = consolidate('ERR_16x8_16QAM_Tap4_Q5_MMSESoft_C',number);
M6C = consolidate('ERR_16x8_16QAM_Tap4_Q6_MMSESoft_C',number);
M7C = consolidate('ERR_16x8_16QAM_Tap4_Q7_MMSESoft_C',number);
M8C = consolidate('ERR_16x8_16QAM_Tap4_Q8_MMSESoft_C',number);
M9C = consolidate('ERR_16x8_16QAM_Tap4_Q9_MMSESoft_C',number);
M10C = consolidate('ERR_16x8_16QAM_Tap4_Q10_MMSESoft_C',number);

O3 = consolidate('ERR_16x8_16QAM_Tap4_Q3_OneShot',number);
O4 = consolidate('ERR_16x8_16QAM_Tap4_Q4_OneShot',number);
O5 = consolidate('ERR_16x8_16QAM_Tap4_Q5_OneShot',number);
O6 = consolidate('ERR_16x8_16QAM_Tap4_Q6_OneShot',number);
O7 = consolidate('ERR_16x8_16QAM_Tap4_Q7_OneShot',number);
O8 = consolidate('ERR_16x8_16QAM_Tap4_Q8_OneShot',number);
O9 = consolidate('ERR_16x8_16QAM_Tap4_Q9_OneShot',number);
O10 = consolidate('ERR_16x8_16QAM_Tap4_Q10_OneShot',number);

% =========================================================================

OP_16x8_MMSE_Q = [ getOperatingPoint(M3,TargetPER) ; ...
    getOperatingPoint(M4,TargetPER) ; ...
    getOperatingPoint(M5,TargetPER) ; ...
    getOperatingPoint(M6,TargetPER) ; ...
    getOperatingPoint(M7,TargetPER) ; ...
    getOperatingPoint(M8,TargetPER) ; ...
    getOperatingPoint(M9,TargetPER) ; ...
    getOperatingPoint(M10,TargetPER) ];

OP_16x8_MMSE_C = [ getOperatingPoint(M3C,TargetPER) ; ...
    getOperatingPoint(M4C,TargetPER) ; ...
    getOperatingPoint(M5C,TargetPER) ; ...
    getOperatingPoint(M6C,TargetPER) ; ...
    getOperatingPoint(M7C,TargetPER) ; ...
    getOperatingPoint(M8C,TargetPER) ; ...
    getOperatingPoint(M9C,TargetPER) ; ...
    getOperatingPoint(M10C,TargetPER) ];


OP_16x8_OS = [ getOperatingPoint(O3,TargetPER) ; ...
    getOperatingPoint(O4,TargetPER) ; ...
    getOperatingPoint(O5,TargetPER) ; ...
    getOperatingPoint(O6,TargetPER) ; ...
    getOperatingPoint(O7,TargetPER) ; ...
    getOperatingPoint(O8,TargetPER) ; ...
    getOperatingPoint(O9,TargetPER) ; ...
    getOperatingPoint(O10,TargetPER) ];

% -- tradeoff plot
figure(1)
plot(getOperatingPoint(SIMOCSIR,TargetPER)*[1 1],[0 100],'g-.','Linewidth',2)
hold on
plot(getOperatingPoint(SIMO,TargetPER)*[1 1],[0 100],'m:','Linewidth',2)
%
plot(OP_16x8_MMSE_Q,[3 4 5 6 7 8 9 10],'bo-','Linewidth',2)
%
plot(OP_16x8_MMSE_C,[3 4 5 6 7 8 9 10 ],'cv-.','Linewidth',2)
%
plot(OP_16x8_OS,[ 3 4 5 6 7 8 9 10 ],'rs--','Linewidth',2)
%
plot(getOperatingPoint(REFMMSE,TargetPER)*[1 1],[0 100],'b:','Linewidth',2)
plot(getOperatingPoint(M10C,TargetPER)*[1 1],[0 100],'c:','Linewidth',2)
plot(getOperatingPoint(O10,TargetPER)*[1 1],[0 100],'r:','Linewidth',2)
hold off
axis([10 22 0 12])
grid on
set(gca,'FontSize',12)
ylabel('Quantization bits Q_b','FontSize',12)
xlabel('Minimum SNR for 1% PER','FontSize',12)
legend('SIMO bound, CSIR','SIMO bound, CHEST','Quantizer, MMSE','Mismatch 1, MMSE','Mismatch 2, MMSE',1)
title('16 BS antennas, 8 users')

%%%
%%%
%%%
B32U8SISOCSIR = consolidate('ERR_32x8_16QAM_Tap4_SIMO_CSIR',number);
B32U8SISO  = consolidate('ERR_32x8_16QAM_Tap4_SIMO',number);

B32U8Q1 = consolidate('ERR_32x8_16QAM_Tap4_Q1_MMSESoft',number);
B32U8Q2 = consolidate('ERR_32x8_16QAM_Tap4_Q2_MMSESoft',number);
B32U8Q3 = consolidate('ERR_32x8_16QAM_Tap4_Q3_MMSESoft',number);
B32U8Q4 = consolidate('ERR_32x8_16QAM_Tap4_Q4_MMSESoft',number);
B32U8Q5 = consolidate('ERR_32x8_16QAM_Tap4_Q5_MMSESoft',number);
B32U8Q6 = consolidate('ERR_32x8_16QAM_Tap4_Q6_MMSESoft',number);
B32U8Q7 = consolidate('ERR_32x8_16QAM_Tap4_Q7_MMSESoft',number);
B32U8Q8 = consolidate('ERR_32x8_16QAM_Tap4_Q8_MMSESoft',number);
B32U8Q9 = consolidate('ERR_32x8_16QAM_Tap4_Q9_MMSESoft',number);
B32U8Q10 = consolidate('ERR_32x8_16QAM_Tap4_Q10_MMSESoft',number);

B32U8OS1 = consolidate('ERR_32x8_16QAM_Tap4_Q1_OneShot',number);
B32U8OS2 = consolidate('ERR_32x8_16QAM_Tap4_Q2_OneShot',number);
B32U8OS3 = consolidate('ERR_32x8_16QAM_Tap4_Q3_OneShot',number);
B32U8OS4 = consolidate('ERR_32x8_16QAM_Tap4_Q4_OneShot',number);
B32U8OS5 = consolidate('ERR_32x8_16QAM_Tap4_Q5_OneShot',number);
B32U8OS6 = consolidate('ERR_32x8_16QAM_Tap4_Q6_OneShot',number);
B32U8OS7 = consolidate('ERR_32x8_16QAM_Tap4_Q7_OneShot',number);
B32U8OS8 = consolidate('ERR_32x8_16QAM_Tap4_Q8_OneShot',number);
B32U8OS9 = consolidate('ERR_32x8_16QAM_Tap4_Q9_OneShot',number);
B32U8OS10 = consolidate('ERR_32x8_16QAM_Tap4_Q10_OneShot',number);

OP_32x8 = [ getOperatingPoint(B32U8Q1,TargetPER) ; ...
    getOperatingPoint(B32U8Q2,TargetPER) ; ...
    getOperatingPoint(B32U8Q3,TargetPER) ; ...
    getOperatingPoint(B32U8Q4,TargetPER) ; ...
    getOperatingPoint(B32U8Q5,TargetPER) ; ...
    getOperatingPoint(B32U8Q6,TargetPER) ; ...
    getOperatingPoint(B32U8Q7,TargetPER) ; ...
    getOperatingPoint(B32U8Q8,TargetPER) ; ...
    getOperatingPoint(B32U8Q9,TargetPER) ; ...
    getOperatingPoint(B32U8Q10,TargetPER) ];

OP_32x8_OS = [ getOperatingPoint(B32U8OS1,TargetPER) ; ...
    getOperatingPoint(B32U8OS2,TargetPER) ; ...
    getOperatingPoint(B32U8OS3,TargetPER) ; ...
    getOperatingPoint(B32U8OS4,TargetPER) ; ...
    getOperatingPoint(B32U8OS5,TargetPER) ; ...
    getOperatingPoint(B32U8OS6,TargetPER) ; ...
    getOperatingPoint(B32U8OS7,TargetPER) ; ...
    getOperatingPoint(B32U8OS8,TargetPER) ; ...
    getOperatingPoint(B32U8OS9,TargetPER) ; ...
    getOperatingPoint(B32U8OS10,TargetPER) ];

% -- tradeoff plot
figure(2)
plot(getOperatingPoint(B32U8SISOCSIR,TargetPER)*[1 1],[0 100],'g-.','Linewidth',2)
hold on
plot(getOperatingPoint(B32U8SISO,TargetPER)*[1 1],[0 100],'m:','Linewidth',2)
%
plot(OP_32x8,[1 2 3 4 5 6 7 8 9 10 ],'bo-','Linewidth',2)
%
plot(OP_32x8_OS,[1 2 3 4 5 6 7 8 9 10 ],'rs--','Linewidth',2)
%
plot(getOperatingPoint(B32U8Q10,TargetPER)*[1 1],[0 100],'b:','Linewidth',2)
plot(getOperatingPoint(B32U8OS10,TargetPER)*[1 1],[0 100],'r:','Linewidth',2)
%
hold off
axis([7 16 0 12])
grid on
set(gca,'FontSize',12)
ylabel('Quantization bits Q_b','FontSize',12)
xlabel('Minimum SNR for 1% PER','FontSize',12)
legend('SIMO bound, CSIR','SIMO bound, CHEST','Quantizer, MMSE','Mismatch 2, MMSE',1)
title('32 BS antennas, 8 users')

%%%
%%%
%%%

B64U8SISOCSIR = consolidate('ERR_64x8_16QAM_Tap4_Q1_SIMO_CSIR',number);
B64U8SISO  = consolidate('ERR_64x8_16QAM_Tap4_Q1_SIMO',number);

B64U8Q1 = consolidate('ERR_64x8_16QAM_Tap4_Q1_MMSESoft',number);
B64U8Q2 = consolidate('ERR_64x8_16QAM_Tap4_Q2_MMSESoft',number);
B64U8Q3 = consolidate('ERR_64x8_16QAM_Tap4_Q3_MMSESoft',number);
B64U8Q4 = consolidate('ERR_64x8_16QAM_Tap4_Q4_MMSESoft',number);
B64U8Q5 = consolidate('ERR_64x8_16QAM_Tap4_Q5_MMSESoft',number);
B64U8Q6 = consolidate('ERR_64x8_16QAM_Tap4_Q6_MMSESoft',number);
B64U8Q7 = consolidate('ERR_64x8_16QAM_Tap4_Q7_MMSESoft',number);
B64U8Q8 = consolidate('ERR_64x8_16QAM_Tap4_Q8_MMSESoft',number);
B64U8Q9 = consolidate('ERR_64x8_16QAM_Tap4_Q9_MMSESoft',number);
B64U8Q10 = consolidate('ERR_64x8_16QAM_Tap4_Q10_MMSESoft',number);

B64U8Q2OS = consolidate('ERR_64x8_16QAM_Tap4_Q2_OneShot',number);
B64U8Q3OS = consolidate('ERR_64x8_16QAM_Tap4_Q3_OneShot',number);
B64U8Q4OS = consolidate('ERR_64x8_16QAM_Tap4_Q4_OneShot',number);
B64U8Q5OS = consolidate('ERR_64x8_16QAM_Tap4_Q5_OneShot',number);
B64U8Q6OS = consolidate('ERR_64x8_16QAM_Tap4_Q6_OneShot_r',number);
B64U8Q7OS = consolidate('ERR_64x8_16QAM_Tap4_Q7_OneShot_r',number);
B64U8Q8OS = consolidate('ERR_64x8_16QAM_Tap4_Q8_OneShot_r',number);
B64U8Q9OS = consolidate('ERR_64x8_16QAM_Tap4_Q9_OneShot_r',number);
B64U8Q10OS = consolidate('ERR_64x8_16QAM_Tap4_Q10_OneShot_r',number);


OP_64x8 = [ getOperatingPoint(B64U8Q1,TargetPER) ; ...
    getOperatingPoint(B64U8Q2,TargetPER) ; ...
    getOperatingPoint(B64U8Q3,TargetPER) ; ...
    getOperatingPoint(B64U8Q4,TargetPER) ; ...
    getOperatingPoint(B64U8Q5,TargetPER) ; ...
    getOperatingPoint(B64U8Q6,TargetPER) ; ...
    getOperatingPoint(B64U8Q7,TargetPER) ; ...
    getOperatingPoint(B64U8Q8,TargetPER) ; ...
    getOperatingPoint(B64U8Q9,TargetPER) ; ...
    getOperatingPoint(B64U8Q10,TargetPER) ];

OP_64x8_OS = [ getOperatingPoint(B64U8Q2OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q3OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q4OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q5OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q6OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q7OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q8OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q9OS,TargetPER) ; ...
    getOperatingPoint(B64U8Q10OS,TargetPER) ];


figure(3)

% -- tradeoff plot
plot(getOperatingPoint(B64U8SISOCSIR,TargetPER)*[1 1],[0 100],'g-.','Linewidth',2)
hold on
plot(getOperatingPoint(B64U8SISO,TargetPER)*[1 1],[0 100],'m:','Linewidth',2)
%
plot(OP_64x8,[1 2 3 4 5 6 7 8 9 10 ],'bo-','Linewidth',2)
%
plot(OP_64x8_OS,[ 2 3 4 5 6 7 8 9 10 ],'rs--','Linewidth',2)
%
plot(getOperatingPoint(B64U8Q10,TargetPER)*[1 1],[0 100],'b:','Linewidth',2)
plot(getOperatingPoint(B64U8Q10OS,TargetPER)*[1 1],[0 100],'r:','Linewidth',2)
%
hold off
axis([4 12 0 12])
grid on
set(gca,'FontSize',12)
ylabel('Quantization bits Q_b','FontSize',12)
xlabel('Minimum SNR for 1% PER','FontSize',12)
legend('SIMO bound, CSIR','SIMO bound, CHEST','Quantizer, MMSE','Mismatch 2, MMSE',1)
title('64 BS antennas, 8 users')


%%%


B128U8SISOCSIR = consolidate('ERR_128x8_16QAM_Tap4_Q1_SIMO_CSIR',number);
B128U8SISO  = consolidate('ERR_128x8_16QAM_Tap4_Q1_SIMO',number);

B128U8Q1 = consolidate('ERR_128x8_16QAM_Tap4_Q1_MMSESoft',number);
B128U8Q2 = consolidate('ERR_128x8_16QAM_Tap4_Q2_MMSESoft_XXX',number);
B128U8Q3 = consolidate('ERR_128x8_16QAM_Tap4_Q3_MMSESoft_XXX',number);
B128U8Q4 = consolidate('ERR_128x8_16QAM_Tap4_Q4_MMSESoft_XXX',number);
B128U8Q5 = consolidate('ERR_128x8_16QAM_Tap4_Q5_MMSESoft_XXX',number);
B128U8Q6 = consolidate('ERR_128x8_16QAM_Tap4_Q6_MMSESoft_XXX',number);
B128U8Q7 = consolidate('ERR_128x8_16QAM_Tap4_Q7_MMSESoft_XXX',number);
B128U8Q8 = consolidate('ERR_128x8_16QAM_Tap4_Q8_MMSESoft_XXX',number);
B128U8Q9 = consolidate('ERR_128x8_16QAM_Tap4_Q9_MMSESoft_XXX',number);
B128U8Q10 = consolidate('ERR_128x8_16QAM_Tap4_Q10_MMSESoft_XXX',number);

B128U8Q2_OS = consolidate('ERR_128x8_16QAM_Tap4_Q2_OneShot_XXX',number);
B128U8Q3_OS = consolidate('ERR_128x8_16QAM_Tap4_Q3_OneShot_XXX',number);
B128U8Q4_OS = consolidate('ERR_128x8_16QAM_Tap4_Q4_OneShot_XXX',number);
B128U8Q5_OS = consolidate('ERR_128x8_16QAM_Tap4_Q5_OneShot_XXX',number);
B128U8Q6_OS = consolidate('ERR_128x8_16QAM_Tap4_Q6_OneShot_XXX',number);
B128U8Q7_OS = consolidate('ERR_128x8_16QAM_Tap4_Q7_OneShot_XXX',number);
B128U8Q8_OS = consolidate('ERR_128x8_16QAM_Tap4_Q8_OneShot_XXX',number);
B128U8Q9_OS = consolidate('ERR_128x8_16QAM_Tap4_Q9_OneShot_XXX',number);
B128U8Q10_OS = consolidate('ERR_128x8_16QAM_Tap4_Q10_OneShot_XXX',number);

OP_128x8 = [ getOperatingPoint(B128U8Q1,TargetPER) ; ...
    getOperatingPoint(B128U8Q2,TargetPER) ; ...
    getOperatingPoint(B128U8Q3,TargetPER) ; ...
    getOperatingPoint(B128U8Q4,TargetPER) ; ...
    getOperatingPoint(B128U8Q5,TargetPER) ; ...
    getOperatingPoint(B128U8Q6,TargetPER) ; ...
    getOperatingPoint(B128U8Q7,TargetPER) ; ...
    getOperatingPoint(B128U8Q8,TargetPER) ; ...
    getOperatingPoint(B128U8Q9,TargetPER) ; ...
    getOperatingPoint(B128U8Q10,TargetPER)  ];

OP_128x8_OS = [ getOperatingPoint(B128U8Q2_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q3_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q4_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q5_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q6_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q7_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q8_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q9_OS,TargetPER) ; ...
    getOperatingPoint(B128U8Q10_OS,TargetPER)];

figure(4)
% -- tradeoff plot

plot(getOperatingPoint(B128U8SISOCSIR,TargetPER)*[1 1],[0 100],'g-.','Linewidth',2)
hold on
plot(getOperatingPoint(B128U8SISO,TargetPER)*[1 1],[0 100],'m:','Linewidth',2)
%
plot(OP_128x8,[1 2 3 4 5 6 7 8 9 10 ],'bo-','Linewidth',2)
%
plot(OP_128x8_OS,[ 2 3 4 5 6 7 8 9 10 ],'rs--','Linewidth',2)
%
plot(getOperatingPoint(B128U8Q10,TargetPER)*[1 1],[0 100],'b:','Linewidth',2)
plot(getOperatingPoint(B128U8Q10_OS,TargetPER)*[1 1],[0 100],'r:','Linewidth',2)
%
hold off
axis([1 7 0 12])
grid on
set(gca,'FontSize',12)
ylabel('Quantization bits Q_b','FontSize',12)
xlabel('Minimum SNR for 1% PER','FontSize',12)
legend('SIMO bound, CSIR','SIMO bound, CHEST','Quantizer, MMSE','Mismatch 2, MMSE',1)
title('128 BS antennas, 8 users')

end



function [SNR] = getOperatingPoint(In,TargetPER)

% -- extract average PER
PER = mean(In.PER,2);

% -- Find the target BER
idx2 = find(PER<TargetPER);
if ~isempty(idx2)
    idx2=idx2(1);
    d_log_PER=log10(PER(idx2-1))-log10(PER(idx2));
    d_SNR = In.TxRx.Sim.SNR_dB_list(idx2-1)-In.TxRx.Sim.SNR_dB_list(idx2);
    SNR = (log10(TargetPER)-log10(PER(idx2-1)))*(d_SNR/d_log_PER)+In.TxRx.Sim.SNR_dB_list(idx2-1);
else
    warning('operating point inexsistent')
    SNR = inf;
end

end

