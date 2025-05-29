
function [LTE_sb_2,LTE_sb_5,small_band_offset] = main_gen_LTE_sig(num_of_small_bands,small_band_num)

%numOFDMsymb = 28;
numOFDMsymb = 14;

plotEnabled = 0;

plotFFTup = 4;

subc_spacing = 60; % kHz % 60 % 15

IRFFilterSubBand = loadResampleFilters768();
res_mat = [ 1 1 ; 4 3; 2 1; 2 1 ];
n_append = 0;
% 
n_order = 1; %4; %10; % order factor n of the resampling filter 2 × n × max(p,q)
              %used, if IRFFilterSubBand = []

EvmCalcOffset = 4; %70 %140; %4; % 208;


g = 0.0; %0.4;

%N = 1;                      %number of subbands (for FFT and cp)
N = 5;                      %number of subbands (for FFT and cp)

numFFTLocal = 384;          % Number of FFT points local

numFFTGlobal = 2048;        % Number of FFT points in the whole band

numRBsLocal = 27;           % Number of resource blocks (15 in LTE)

%cpLenLocal = 39;            % Cyclic prefix length in samples
%cpLenLocal = 33;            % Cyclic prefix length in samples
cpLenLocal = 27;            % Cyclic prefix length in samples
%cpLenLocal = 0;            % Cyclic prefix length in samples

BandShift = [ 648 324 0 -324 -648 ];
%BandShift = [ 2*333 333 0 -333 -2*333 ];

%BandShift = [ -648 ];
%BandShift = [ -648 324 ];

numFFTLocal = 60 / subc_spacing * numFFTLocal;

numFFTGlobal = 60 / subc_spacing * numFFTGlobal;

numRBsLocal = 60 / subc_spacing * numRBsLocal;

cpLenLocal = 60 / subc_spacing * cpLenLocal;

BandShift = 60 / subc_spacing * BandShift;

bitsPerSubCarrier = 2;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM

UpRate = 16;

DownRate = 3;
        
rbSize = 12;                % Number of subcarriers per resource block

plt_off = 0;

numDataCarriersLocal = numRBsLocal*rbSize;    % number of data subcarriers in subband

cpLenGlobal = cpLenLocal * UpRate / DownRate;

numDataCarriersSmallBand = numDataCarriersLocal/num_of_small_bands;

if mod(numDataCarriersSmallBand,1)
    
    error('Wrong number of small bands');
    
end


small_band_offset = (small_band_num - 0.5 * (num_of_small_bands + 1)) * ...
            numDataCarriersSmallBand * subc_spacing * 1e3;
    
% QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator( ...
    'ModulationOrder', 2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% EVM object
evm_obj = comm.EVM;


% Generate data symbols

rng(1234);


symbolsInAcc = [];

rxSigAcc = [];

txSigOFDMSumAcc = [];

txSigOFDMAcc = [];

BufferLength = 28; %28; %28; % Buffer length

txBuffer = zeros( BufferLength * ( numFFTLocal + cpLenLocal ), N);

txBufferPnt = 1;

tmp0 = zeros(numDataCarriersLocal,1);
%tmp1 = ones(numDataCarriersLocal,1);
tmp1 = [ones(numDataCarriersSmallBand,1);zeros(numDataCarriersLocal-numDataCarriersSmallBand,1)];
tmp1 = circshift(tmp1,(small_band_num-1)*numDataCarriersSmallBand);
tmp01 =[tmp0;tmp1;tmp0;tmp0;tmp1];

for k = 1 : numOFDMsymb

    % Symbols gen

    bitsIn = randi( [0 1], N*bitsPerSubCarrier*numDataCarriersLocal, 1 );

    symbolsIn = step( qamMapper, bitsIn );
    
    symbolsIn = symbolsIn .* tmp01;
                   
    symbolsInAcc = [ symbolsInAcc; symbolsIn ];
    
    symbolsIn = reshape( symbolsIn, [ numDataCarriersLocal N ]);
        
    txSigOFDUpOld = zeros( numFFTGlobal + cpLenGlobal, N );
    
    txSigOFDMup = zeros( numFFTGlobal + cpLenGlobal, N );
        
    intSig = zeros( numFFTGlobal + cpLenGlobal, N );

    intSigTmp = zeros( numDataCarriersLocal, N, N );
    
    txSigOFDMSum = 0;
    
    txSigOFDUpR = [];
        
    for sB = 1 : N
        
        % Phase compensation
        
        txSigOFDMtmp = symbolsIn( :, sB ) * exp( 2 * pi * 1i * 2 * cpLenLocal * BandShift(sB) / numFFTLocal );
        
        txSigOFDM = ifft_cp( txSigOFDMtmp, numFFTLocal, 'middle', cpLenLocal );
        
        txBuffer( ( txBufferPnt - 1 ) * ( numFFTLocal + cpLenLocal ) + 1 : txBufferPnt * ( numFFTLocal + cpLenLocal ), sB ) = txSigOFDM;
                               
         txSigOFDUp = slice_cont_flt( txBuffer( :, sB ), txBufferPnt, numFFTLocal, cpLenGlobal, BandShift(sB), plotEnabled, plotFFTup, UpRate, DownRate, numFFTGlobal, IRFFilterSubBand, res_mat, n_append, n_order );
         
         txSigOFDMSum = txSigOFDMSum + txSigOFDUp;
         
         txSigOFDUpR(:,sB) = txSigOFDUp( 2 * cpLenGlobal + numFFTGlobal - EvmCalcOffset + 1 : end - EvmCalcOffset);
                        
    end
    
    txSigOFDMAcc = [ txSigOFDMAcc; txSigOFDUpR];

    txBufferPnt = txBufferPnt + 1;
    
    % Clean buffer when full
    if ~mod( k, BufferLength )
        
        txBuffer = zeros( BufferLength * ( numFFTLocal + cpLenLocal ), N);
        
        txBufferPnt = 1;
    
    end
    
    txSigOFDMSumAcc = [ txSigOFDMSumAcc; txSigOFDMSum( 2 * cpLenGlobal + numFFTGlobal - EvmCalcOffset + 1 : end - EvmCalcOffset)];

    txSigOFDMTranc = txSigOFDMSum( 2 * cpLenGlobal + numFFTGlobal - EvmCalcOffset + 1 : end - EvmCalcOffset);

    rxSigTmp = fft_cp_old( txSigOFDMTranc, numFFTGlobal, N*numDataCarriersLocal, EvmCalcOffset, BandShift );
    
    % Phase correction at receiver
    

    rxSigAcc = [ rxSigAcc; rxSigTmp  ];

end


symbols_resh = reshape( symbolsInAcc, [ rbSize N*numRBsLocal numOFDMsymb ]  );

symbols_resh = permute( symbols_resh, [ 1 3 2 ] );

symbols_resh = reshape( symbols_resh, [ rbSize*numOFDMsymb N*numRBsLocal ]);

rxSig_resh = reshape( rxSigAcc, [ rbSize N*numRBsLocal numOFDMsymb ]  );

rxSig_resh = permute( rxSig_resh, [ 1 3 2 ] );

rxSig_resh = reshape( rxSig_resh, [ rbSize*numOFDMsymb N*numRBsLocal ]);

evmes = [];

for k = 1 : size(symbols_resh,2)
    
    evmes_tmp = step( evm_obj, symbols_resh(:,k), numFFTLocal / numFFTGlobal * rxSig_resh(:,k) );
    
    evmes = [ evmes evmes_tmp ];
    
end

% figure; plot( 1 : 5 / N : 5 / N * length( evmes ) + plt_off, evmes );

max_evmes = max(evmes);

rxSig_resh = numFFTLocal / numFFTGlobal * rxSig_resh;


release(evm_obj);

step( evm_obj, symbols_resh(:), rxSig_resh(:) );

ang = angle(symbols_resh(:) - rxSig_resh(:));

symbolsIn = reshape( symbolsIn, [rbSize N * numRBsLocal] );    

figure; PSDoFreqplot(1e5*txSigOFDMSumAcc,'r',122.88e6,100);

LTE_sb_2 = squeeze(txSigOFDMAcc(:,2));
% save(out_file_name_1,'LTE_sb_2');
LTE_sb_5 = squeeze(txSigOFDMAcc(:,5));
% save(out_file_name_2,'LTE_sb_5');

end