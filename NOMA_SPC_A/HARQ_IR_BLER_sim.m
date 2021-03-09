clc; clear variables; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power from 0 to 40 dbm 
% Send 1000 blocks at every power level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic

%%%%%% Transmit power setting %%%%%%
% Transmit power in dBm
Pt = 50:5:100;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);   

%%%%% Number of block sent %%%%%
nBlock = 10000;

%%%%%% User distance setting %%%%%%
d1 = 400;
eta = 4;

%%%%% DLSCH paramemter setting %%%%%
dlsch.blkLen = 5120;
dlsch.targetCodeRate = 567/1024;
dlsch.rv = [0 2 3 1];
dlsch.mod = '16QAM';
dlsch.nLayers = 1;
dlsch.outlen = 10240;

% DLSCH encoder setting 
dlschEncoder = nrDLSCH;
dlschEncoder.MultipleHARQProcesses = false;
dlschEncoder.TargetCodeRate = dlsch.targetCodeRate;

% DLSCH decoder setting 
dlschDecoder = nrDLSCHDecoder;
dlschDecoder.TargetCodeRate = dlsch.targetCodeRate;
dlschDecoder.TransportBlockLength = dlsch.blkLen;
dlschDecoder.MaximumLDPCIterationCount = 6;

%%%%% PDSCH paramemter setting %%%%%
pdsch.ncellid = 42;
pdsch.rnti = 6143;
pdsch.mod = dlsch.mod;
pdsch.nLayers = dlsch.nLayers;




%%%%% Noise Generation %%%%%
BW = 10^7;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)


totalBlockerr = zeros(length(Pt),1);

for u = 1:length(Pt)
    u
    for mm = 1:nBlock
        % Channel Generation
        N =2560;
        h1 = sqrt(d1^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
        % Generate noise samples for both users
        w1 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
        
        % Generate transport block
        trBlk = randi([0 1],dlsch.blkLen,1,'int8');

        % Load transport block into encoder
        setTransportBlock(dlschEncoder, trBlk);
        rv_idx = 1;
        
        while 1
            % DLSCH encoding
            codedTrBlock = dlschEncoder(dlsch.mod, dlsch.nLayers, dlsch.outlen, dlsch.rv(rv_idx));

            % PDSCH generation
            xmod = nrPDSCH(codedTrBlock, pdsch.mod, pdsch.nLayers, pdsch.ncellid, pdsch.rnti);

            % Do superposition coding
            x = sqrt(pt(u))*xmod;

            % Received signals
            y1 = (h1.').*x + w1.';

            % Equalize 
            y_eq = y1./h1.';
%             y_eq = awgn(xmod,Pt(u),'measured');
%             sigp = 10*log10(norm(y_eq,2)^2/numel(y_eq));
%             snr = sigp - 1;
%             noisep = 10^(snr/10);

%             [ydemod, symbol] = nrPDSCHDecode(y_eq, pdsch.mod, pdsch.ncellid, pdsch.rnti, noisep);
            [ydemod, symbol] = nrPDSCHDecode(y_eq, pdsch.mod, pdsch.ncellid, pdsch.rnti, var(w1));

            [decbits,blkerr] = dlschDecoder(ydemod, dlsch.mod, dlsch.nLayers, dlsch.rv(rv_idx));
            bler(u,mm,rv_idx) = blkerr;

%             
%             if  dlsch.rv(rv_idx) == 1
%                 resetSoftBuffer (dlschDecoder,1);
%                 break;

            if  rv_idx == 4 % maximum transmission time
                resetSoftBuffer (dlschDecoder,0);
                break;
            elseif blkerr == 0
                resetSoftBuffer (dlschDecoder,0);
                break;
            else
                rv_idx = rv_idx + 1;
            end

        end
        totalBlockerr(u) = totalBlockerr(u) + blkerr;
    end

end

toc

blerr = totalBlockerr./nBlock

figure (1)
semilogy(Pt, blerr,'*r');

title('BLER vs Transmit Power');
xlabel('Transmit power (P in dBm)');
ylabel('BLER');

figure (2)
plot(Pt, blerr,'*r');

title('BLER vs Transmit Power');
xlabel('Transmit power (P in dBm)');
ylabel('BLER');
