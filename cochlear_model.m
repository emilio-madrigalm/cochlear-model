% COCHLEAR_MODEL  Cascaded filter-bank model of the human cochlea.
%   512 sections, each a pole-pair resonator (displacement) followed by a
%   zero-pair notch (pressure). Spatial differentiation sharpens tuning.
%   Run this file; figures are saved to ./figures.
clear; clc; close all;
outDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~isfolder(outDir), mkdir(outDir); end

%% Calculate the Filter Frequencies

N  = 512;
dx = 3.5/N;                        % = 0.00684 cm
x  = (0:N-1)*dx;                 % position of each piece along the membrane
fp = 7500*10.^(-0.667*x);

fz = zeros(1,N);
fz(1) = 7580;
for i = 2:N
    fz(i) = 0.985*fp(i-1); 
end                              % the loop sets the notch frequency for all 512 filter sections.


%% Quality Factors and Bandwidths 

Qp = 8;                           % poles Q factor
Qz = 2*Qp;                        % zeros Q factor = 16
BWp = fp/Qp;                         % poles bandwidth 
BWz = fz/Qz;                       % zeros bandwidth

%% Filter Coefficients and Gain Normalisation

fs = 16000; % sampling freq

% resonator (poles)
thp = 2*pi*fp/fs;              % converts each Hz freq into an angle
rp  = 1 - pi*BWp/fs;              % converts each bw into a pole radius
b1  = 2*rp.*cos(thp);            % coefficient 1
b2  = rp.^2;                     % coefficient 2
Gp  = 1 - b1 + b2;             % gain (= 1 at DC) 

% notch (zeros)
thz = 2*pi*fz/fs;             
rz  = 1 - pi*BWz/fs;         
a1  = 2*rz.*cos(thz);           % coefficient 1
a2  = rz.^2;                  % coefficient 2
Gz  = 1./(1 - a1 + a2);        % gain (1 at DC)

%% Time-Domain Simulation

% Cascade with impulse input
L   = 16000;                    % 1 s, long enough for apex sections to decay
imp = [1; zeros(L-1,1)];           % unit impulse
[D, P] = runCochlea(imp, Gp, b1, b2, Gz, a1, a2);

% Spatial differentiation
d = D(:,1:end-1) - D(:,2:end);      % 1st stage: d_m[n] = s_m[n] - s_{m+1}[n]
e = d(:,1:end-1) - d(:,2:end);      % 2nd stage: e_m[n] = d_m[n] - d_{m+1}[n]

% Impulse responses: before vs after differentiation
secs = [88 226 422];             % sections to show
lens = [300 1500 5000];            % samples to show for each
figure;
for k = 1:3
    s = secs(k);
    subplot(3,2,2*k-1);
    plot(D(1:lens(k), s));
    title(sprintf('Section %d (%.1f Hz): no differentiation', s, fp(s)));
    ylabel('Displacement');
    subplot(3,2,2*k);
    plot(e(1:lens(k), s));
    title(sprintf('Section %d: two differentiations', s));
end
subplot(3,2,5); xlabel('Sample number');
subplot(3,2,6); xlabel('Sample number');
exportgraphics(gcf, fullfile(outDir, 'impulse_responses.png'), 'Resolution', 150);

% Magnitude responses via FFT
fAx  = (0:L-1)*fs/L;               % frequency of each FFT bin (Hz)
keep = 2:L/2;                    % positive freqs, skip 0 Hz (log axis can't show it)
s    = 192;                      % section tuned near 1 kHz
sigs  = {D(:,s), d(:,s), e(:,s)};
names = {'Before', 'After first', 'After second'};
figure;
for k = 1:3
    H   = abs(fft(sigs{k}));          % magnitude response
    HdB = 20*log10(H/max(H));       % normalise peak to 0 dB
    subplot(3,1,k);
    semilogx(fAx(keep), HdB(keep));
    grid on; ylim([-100 5]); xlim([5 2000]);
    title(sprintf('%s spatial differentiation: section %d (%.2f kHz)', ...
          names{k}, s, fp(s)/1000));
    ylabel('Magnitude (dB)');
end
xlabel('Frequency (Hz)');
exportgraphics(gcf, fullfile(outDir, 'magnitude_response_1kHz.png'), 'Resolution', 150);

% Sinusoid input
f0    = 1000;                       % input frequency (Hz) - change to test others
Ls    = 2000;
n     = (0:Ls-1)';
sinIn = sin(2*pi*f0*n/fs);
[Ds, Ps] = runCochlea(sinIn, Gp, b1, b2, Gz, a1, a2);
ds = Ds(:,1:end-1) - Ds(:,2:end);
es = ds(:,1:end-1) - ds(:,2:end);

% One time instant, after transient has decayed
row = 1000;
figure;
subplot(2,1,1); plot(Ds(row,:));
title(sprintf('Displacement, no differentiation (%d Hz, sample %d)', f0, row));
ylabel('Displacement'); xlim([1 N]);
subplot(2,1,2); plot(es(row,:));
title('Displacement, two spatial differentiations');
xlabel('Section number'); ylabel('Displacement'); xlim([1 N]);
exportgraphics(gcf, fullfile(outDir, 'tone_snapshot_1kHz.png'), 'Resolution', 150);

% Several time instants superimposed -> envelope
rows = 1000:1015;                   % 16 samples = one full period of 1 kHz
figure;
subplot(1,2,1); plot(Ps(rows,:)');
title('Pressure (superimposed)'); xlabel('Section number'); ylabel('Pressure'); xlim([1 N]);
subplot(1,2,2); plot(es(rows,:)');
title('Displacement, two differentiations (superimposed)'); xlabel('Section number'); xlim([1 N]);
exportgraphics(gcf, fullfile(outDir, 'travelling_wave_envelope.png'), 'Resolution', 150);

%% Final Validation
% Three-tone input
f3  = [300 1000 4000];              % low, middle, high frequency (Hz)
L3  = 3000;                      % number of samples
n3  = (0:L3-1)';                   % sample numbers, as a column
in3 = sin(2*pi*f3(1)*n3/fs) + sin(2*pi*f3(2)*n3/fs) + sin(2*pi*f3(3)*n3/fs);
[D3, P3] = runCochlea(in3, Gp, b1, b2, Gz, a1, a2);

% Displacement vs section number at one time instant
d3 = D3(:,1:end-1) - D3(:,2:end);   % one stage of spatial differentiation
e3 = d3(:,1:end-1) - d3(:,2:end);   % two stages

row3 = 2500;                   % time instant (transient has decayed)
figure;
subplot(3,1,1); plot(D3(row3,:));
title(sprintf('Before spatial differentiation (sample %d)', row3));
ylabel('Displacement'); xlim([1 N]);
subplot(3,1,2); plot(d3(row3,:));
title('After one stage of spatial differentiation');
ylabel('Displacement'); xlim([1 N]);
subplot(3,1,3); plot(e3(row3,:));
title('After two stages of spatial differentiation');
xlabel('Section number'); ylabel('Displacement'); xlim([1 N]);
exportgraphics(gcf, fullfile(outDir, 'three_tone_separation.png'), 'Resolution', 150);

% Magnitude responses of ten sections
secs10 = round(linspace(30, 480, 10));    % ten sections spread base -> apex
figure;
for k = 1:10
    s = secs10(k);

    H = abs(fft(D(:,s)));              % before differentiation
    subplot(2,1,1);
    semilogx(fAx(keep), 20*log10(H(keep)/max(H))); hold on;

    H = abs(fft(e(:,s)));               % after two stages
    subplot(2,1,2);
    semilogx(fAx(keep), 20*log10(H(keep)/max(H))); hold on;
end
subplot(2,1,1);
title('(a) Before spatial differentiation');
ylabel('Normalised magnitude (dB)');
grid on; xlim([20 8000]); ylim([-80 5]);
subplot(2,1,2);
title('(b) After two stages of spatial differentiation');
xlabel('Frequency (Hz)'); ylabel('Normalised magnitude (dB)');
grid on; xlim([20 8000]); ylim([-80 5]);
exportgraphics(gcf, fullfile(outDir, 'filterbank_responses.png'), 'Resolution', 150);

% Switch the 1000 Hz tone off at sample 1500
nOff = 1500;
gate = double(n3 < nOff);              % 1 before sample 1500, 0 after
in6  = sin(2*pi*300*n3/fs) + gate.*sin(2*pi*1000*n3/fs) + sin(2*pi*4000*n3/fs);
[D6, P6] = runCochlea(in6, Gp, b1, b2, Gz, a1, a2);
d6 = D6(:,1:end-1) - D6(:,2:end);
e6 = d6(:,1:end-1) - d6(:,2:end);

secs6  = [54 186 301];             % the places of 4000, 1000 and 300 Hz
names6 = {'4000 Hz place', '1000 Hz place (switched off)', '300 Hz place'};
figure;
for k = 1:3
    subplot(3,1,k); plot(e6(:, secs6(k))); xline(nOff, 'r');
    title(sprintf('Section %d: %s', secs6(k), names6{k}));
    ylabel('Displacement');
end
xlabel('Sample number');
exportgraphics(gcf, fullfile(outDir, 'tone_offset_response.png'), 'Resolution', 150);


function [D, P] = runCochlea(in, Gp, b1, b2, Gz, a1, a2)
    L = length(in);  N = length(Gp);
    D = zeros(L, N);  P = zeros(L, N);
    for i = 1:N
        D(:,i) = filter(Gp(i), [1, -b1(i), b2(i)], in);          % resonator -> displacement
        P(:,i) = filter(Gz(i)*[1, -a1(i), a2(i)], 1, D(:,i));      % notch -> pressure
        in = P(:,i);                                          % feeds next section
    end
end