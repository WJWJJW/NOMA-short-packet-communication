clc; clear variables; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power from 0 to 40 dbm 
% Send 1000 blocks at every power level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic

%%%%%% Transmit power setting %%%%%%
% Transmit power in dBm
Pt = 30:2:50;                
% Transmit power in linear scale
pt = (10^-3)*10.^(Pt/10);   

%%%%% Number of block sent %%%%%
nBlock = 100;

%%%%%% User distance setting %%%%%%
d1 = 400;
d2 = 700;
a1 = 0.001;
a2 = 1 - a1;
eta = 4;

%%%%% DLSCH paramemter setting %%%%%
dlsch.blkLen = 5120;
dlsch.targetCodeRate = 567/1024;
dlsch.rv = [0 2 3 1];
dlsch.mod = '16QAM';
dlsch.nLayers = 1;
dlsch.outlen = 10240;

% DLSCH encoder setting 
% User 1
dlschEncoder1 = nrDLSCH;
dlschEncoder1.MultipleHARQProcesses = false;
dlschEncoder1.TargetCodeRate = dlsch.targetCodeRate;

% User 2
dlschEncoder2 = nrDLSCH;
dlschEncoder2.MultipleHARQProcesses = false;
dlschEncoder2.TargetCodeRate = dlsch.targetCodeRate;

% DLSCH decoder setting 
% User 1
dlschDecoder1 = nrDLSCHDecoder;
dlschDecoder1.TargetCodeRate = dlsch.targetCodeRate;
dlschDecoder1.TransportBlockLength = dlsch.blkLen;
dlschDecoder1.MaximumLDPCIterationCount = 6;

dlschDecoder12 = nrDLSCHDecoder;
dlschDecoder12.TargetCodeRate = dlsch.targetCodeRate;
dlschDecoder12.TransportBlockLength = dlsch.blkLen;
dlschDecoder12.MaximumLDPCIterationCount = 6;

% User 2
dlschDecoder2 = nrDLSCHDecoder;
dlschDecoder2.TargetCodeRate = dlsch.targetCodeRate;
dlschDecoder2.TransportBlockLength = dlsch.blkLen;
dlschDecoder2.MaximumLDPCIterationCount = 6;

%%%%% PDSCH paramemter setting %%%%%
pdsch.ncellid = 42;
pdsch.rnti = 6143;
pdsch.mod = dlsch.mod;
pdsch.nLayers = dlsch.nLayers;




%%%%% Noise Generation %%%%%
BW = 10^6;                  %System bandwidth
No = -174 + 10*log10(BW);   %Noise power (dBm)
no = (10^-3)*10.^(No/10);   %Noise power (linear scale)


totalBlockerr12 = zeros(length(Pt),1);
totalBlockerr11 = zeros(length(Pt),1);
totalBlockerr2 = zeros(length(Pt),1);

bler12 = zeros(length(Pt), nBlock, 4);
bler11 = zeros(length(Pt), nBlock, 4);
bler2 = zeros(length(Pt), nBlock, 4);

for u = 1:length(Pt)
    u
    for mm = 1:nBlock
        % Channel Generation
        N =2560;
        h1 = sqrt(d1^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
        h2 = sqrt(d2^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
        % Generate noise samples for both users
        w1 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
        w2 = sqrt(no)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
        
        % Generate transport block
        trBlk1 = randi([0 1],dlsch.blkLen,1,'int8');
        trBlk2 = randi([0 1],dlsch.blkLen,1,'int8');

        % Load transport block into encoder
        setTransportBlock(dlschEncoder1, trBlk1);
        setTransportBlock(dlschEncoder2, trBlk2);
        rv_idx = 1;
        blkerr11 = 1;
        blkerr2 = 1;
        
        while 1
            % DLSCH encoding
            codedTrBlock1 = dlschEncoder1(dlsch.mod, dlsch.nLayers, dlsch.outlen, dlsch.rv(rv_idx));
            codedTrBlock2 = dlschEncoder2(dlsch.mod, dlsch.nLayers, dlsch.outlen, dlsch.rv(rv_idx));

            % PDSCH generation
            x1 = nrPDSCH(codedTrBlock1, pdsch.mod, pdsch.nLayers, pdsch.ncellid, pdsch.rnti);
            x2 = nrPDSCH(codedTrBlock2, pdsch.mod, pdsch.nLayers, pdsch.ncellid, pdsch.rnti);

            % Do superposition coding
            x = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2);

            % Received signals
            y1 = (h1.').*x + w1.';
            y2 = (h2.').*x + w2.';

            % Equalize 
            y_eq1 = y1./h1.';
            y_eq2 = y2./h2.';
            
            
            % User 1
            % Decoding user 2 data first
            % PDSCH decoding
            if blkerr11 == 1
                [dlschLLR12, symbol12] = nrPDSCHDecode(y_eq1, pdsch.mod, pdsch.ncellid, pdsch.rnti, var(w1));

                % DLSCH decoding
                [decbits12,blkerr12] = dlschDecoder12(dlschLLR12, dlsch.mod, dlsch.nLayers, dlsch.rv(rv_idx));
                bler12(u,mm,rv_idx) = blkerr12;

                % Subtract user 2's data form received signal
                
                % DLSCH encoding
                setTransportBlock(dlschEncoder1, decbits12);
                codedTrBlock12 = dlschEncoder1(dlsch.mod, dlsch.nLayers, dlsch.outlen, dlsch.rv(rv_idx));

                % PDSCH encoding
                x12 = nrPDSCH(codedTrBlock12, pdsch.mod, pdsch.nLayers, pdsch.ncellid, pdsch.rnti);
               
                % Substraction
                y_s = y_eq1 - sqrt(a2*pt(u)).*x12;

                % Decoding user 1 data directly
                % PDSCH decoding
                [dlschLLR11, symbol11] = nrPDSCHDecode(y_s, pdsch.mod, pdsch.ncellid, pdsch.rnti, var(w1));

                % DLSCH decoding
                [decbits11,blkerr11] = dlschDecoder1(dlschLLR11, dlsch.mod, dlsch.nLayers, dlsch.rv(rv_idx));

                bler11(u,mm,rv_idx) = blkerr11;
            end
            
            if blkerr2 == 1
                % User 2
                % PDSCH decoding
                [dlschLLR2, symbol2] = nrPDSCHDecode(y_eq2, pdsch.mod, pdsch.ncellid, pdsch.rnti, var(w2));

                % DLSCH decoding
                [decbits2,blkerr2] = dlschDecoder2(dlschLLR2, dlsch.mod, dlsch.nLayers, dlsch.rv(rv_idx));
                bler2(u,mm,rv_idx) = blkerr2;
            end
            
            if  rv_idx == 3 % maximum transmission time
                resetSoftBuffer (dlschDecoder1,0);
                resetSoftBuffer (dlschDecoder2,0);
                resetSoftBuffer (dlschDecoder12,0);
                break;
            elseif blkerr11 == 0 && blkerr2 == 0
                resetSoftBuffer (dlschDecoder1,0);
                resetSoftBuffer (dlschDecoder2,0);
                resetSoftBuffer (dlschDecoder12,0);
                break;
            else
                rv_idx = rv_idx + 1;
            end

        end
        totalBlockerr11(u) = totalBlockerr11(u) + blkerr11;
        totalBlockerr2(u) = totalBlockerr2(u) + blkerr2;
    end

end

toc

blErr1 = totalBlockerr11./nBlock;
blErr2 = totalBlockerr2./nBlock;

figure (1)
semilogy(Pt, blErr1,'or', 'linewidth', 1.5);
hold on;grid on;
semilogy(Pt, blErr2,'ob', 'linewidth', 1.5);

title('BLER vs Transmit Power');
xlabel('Transmit power (P in dBm)');
ylabel('BLER');
legend('Sim. User 1/Near user','Sim. User 2/Far user');

figure (2)
plot(Pt, blErr1,'*r', 'linewidth', 1.5);
hold on;grid on;
plot(Pt, blErr2,'*b', 'linewidth', 1.5);

title('BLER vs Transmit Power');
xlabel('Transmit power (P in dBm)');
ylabel('BLER');
legend('Sim. User 1/Near user','Sim. User 2/Far user');





