% -------------------------------------------------------------------------
% Script for simulating virtual analog model of Buchla 259 wavefolder [1].
%
% This script considers the use of the bandlimited ramp function (BLAMP)
% to reduce aliasing present in the wavefolder output due to the 
% discontinuities in the signal derivative. BLAMP is implemented as a
% 2-point polyBLAMP (polynomial approximation of BLAMP) residual 
% compensation [2].
%
% Antialiased model still requires oversampling (by around 8x) to obtain 
% desired signal quality (SNR > 60 dB for all frequencies up to 5 kHz and
% input amplitudes up to 5V as in [1]) but there is a significant reduction 
% compared to the regular algorithm: 64x oversampling loses to 
% 8x + antialiasing by around 15-20 dB in SNR! 
%
% Here SNR is defined as a power ratio between desired harmonics and
% aliased components.
%
% References:
% [1] F. Esqueda, H. Pöntynen, J. D. Parker and V. Välimäki, "Virtual 
%     Analog Buchla 259 Wavefolder", in Proc. 20th Int. Conf. Digital Audio 
%     Effects (DAFx-17), Edinburgh, UK, Sept. 5-9, pp. 192-199.
%     https://dafx17.eca.ed.ac.uk/papers/DAFx17_paper_82.pdf
% [2] Esqueda, F., Bilbao, S., & Valimaki, V. (2016). Aliasing reduction in 
%     clipped signals. IEEE Transactions on Signal Processing, 64(20), 
%     5255-5267. 
%     https://doi.org/10.1109/TSP.2016.2585091
%
% Author: Victor Zheleznov
% Date: 29/03/2023
% -------------------------------------------------------------------------

clc; close all; clear all;

% simulation parameters
audio_rate = 44.1e3; % default audio rate [Hz]
factor = 8;          % oversampling factor
dur = 1;             % duration of simulation [sec]

% input - create sine wave input
f0 = 1e3; % frequency [Hz]
amp = 5;  % amplitude

% physical parameters
Vs = 6;                                          % op-amp supply voltage [V]
r = [10e3,   100e3, 100e3;                       % resistances for folding cells [Ohms]
     49.9e3, 100e3, 43.2e3;
     91e3,   100e3, 56e3;
     30e3,   100e3, 68e3;
     68e3,   100e3, 33e3];
w = -[12; 27.777; 21.428; -17.647; -36.363; -5]; % output weights
rf2 = 1.2e6;                                     % low-pass filter resistance [Ohms] 
c = 100e-12;                                     % low-pass filter capacitance [F]

% calculate derived parameters
SR = factor*audio_rate;    % sample rate [Hz]
Nf = round(SR*dur);        % number of samples in simulation
k = 1/SR;                  % time step [sec]
tvec = dur*(0:Nf-1).'/Nf;  % time vector [sec]
u = amp*sin(2*pi*f0*tvec); % input vector

% regular model
V = zeros(Nf,6);
for i = 1:size(V,2)-1
    V(:,i) = folding_cell(u, amp, f0, SR, r(i,:), Vs, false);
end
V(:,6) = u;

% antialiased model
Vaa = zeros(Nf,6);
for i = 1:size(Vaa,2)-1
    Vaa(:,i) = folding_cell(u, amp, f0, SR, r(i,:), Vs, true);
end
Vaa(:,6) = u;

% mixing stage
Vmix_aa = Vaa*w;
Vmix = V*w;

% filtering stage
wc = 1/(rf2*c);
b = [wc*k/(2 + wc*k), wc*k/(2 + wc*k)];
a = [1, (wc*k-2)/(wc*k+2)];
Vout_aa = filter(b,a,Vmix_aa);
Vout = filter(b,a,Vmix);

% plot time domain
fig_t = figure;
hold on;
plot(tvec, u, 'k--');
plot(tvec, Vout, 'b');
plot(tvec, Vout_aa, 'r');
xlim([0 2/f0])
xlabel('Time [sec]', 'Interpreter', 'latex');
ylabel('Amplitude [V]', 'Interpreter', 'latex');
title(sprintf('Buchla Wavefolder waveform with a sine wave input signal with amplitude $A = %.1f$ and frequency $f = %.1f$ Hz', amp, f0), 'Interpreter', 'latex');
legend({'Input', 'Regular waveform', 'Antialiasing with polyBLAMP'}, 'Interpreter', 'latex');

% calculate FFT parameters
NFFT = 2^(ceil(log(Nf)/log(2)));         % fft size
NFFT_2 = NFFT / 2 + 1;
L = Nf;                                  % window length
win = 0.5*(1 - cos(2*pi.*(0:L-1)./L)).'; % hanning window

% calculate FFT
y_frame = Vout.*win;
Y = fft(y_frame, NFFT);
Y = Y(1:NFFT_2);
yaa_frame = Vout_aa.*win;
Yaa = fft(yaa_frame, NFFT);
Yaa = Yaa(1:NFFT_2);

% cut to audio rate
fmax = audio_rate/2;
fvec = (0:SR/NFFT:fmax).';
Nmax = length(fvec);
Y_dB = 20*log10(abs(Y(1:Nmax))./max(abs(Y(1:Nmax))));
Yaa_dB = 20*log10(abs(Yaa(1:Nmax))./max(abs(Yaa(1:Nmax))));

% plot frequency domain
fig_fft = figure;
hold on
plot(fvec, Y_dB, 'b', fvec, Yaa_dB, 'r');
xlim([0 fmax]);
ylim([-100 5]);
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Magnitude [dB]', 'Interpreter', 'latex');
title(sprintf('Buchla Wavefolder spectrum with a sine wave input signal with amplitude $A = %.1f$ and frequency $f = %.1f$ Hz', amp, f0), 'Interpreter', 'latex');

% find odd harmonics and aliased components
odd_harm = (f0:2*f0:fmax);
odd_harm_idx = find_peaks(Yaa_dB, 100, -100);       % find high spectrum peaks
diff = abs(fvec(odd_harm_idx)-odd_harm) < SR/NFFT;  % check odd harmonics values
odd_harm_idx = odd_harm_idx(sum(diff,2) == 1);      % filter odd harmonics indexes
alias_idx_aa = find_peaks(Yaa_dB, 10, -1000);       % find aliased spectrum components
alias_idx_aa = setdiff(alias_idx_aa, odd_harm_idx); 
alias_idx = find_peaks(Y_dB, 10, -1000);      
alias_idx = setdiff(alias_idx, odd_harm_idx);

% check harmonics and aliased components on spectrum graph
plot(fvec(odd_harm_idx), Yaa_dB(odd_harm_idx), 'ro');
plot(fvec(alias_idx), Y_dB(alias_idx), 'bx');
plot(fvec(alias_idx_aa), Yaa_dB(alias_idx_aa), 'rx');

% add legend
legend({'Regular waveform', 'Antialiasing with polyBLAMP',... 
        'Desired harmonics', 'Aliased components', 'Aliased components'},... 
        'Interpreter', 'latex');

% calculate SNR
SNR_aa = 10*log10(sum(abs(Yaa(odd_harm_idx)).^2) / sum(abs(Yaa(alias_idx_aa)).^2));
SNR = 10*log10(sum(abs(Y(odd_harm_idx)).^2) / sum(abs(Y(alias_idx)).^2));
disp("SNR (regular) = " + SNR + " dB");
disp("SNR (antialiased) = " + SNR_aa + " dB");

% save .wav
audiowrite("buchla_wavefolder_input.wav", u./(max(abs(u))), SR);
audiowrite("buchla_wavefolder_output.wav", Vout_aa./(max(abs(Vout_aa))), SR);

%% FUNCTIONS
% function to calculate output of a single folding cell
% input:
%   u - input signal (sine with constant frequency and amplitude);
%   amp - sine amplitude;
%   f0 - sine frequency;
%   SR - sampling rate;
%   r - resistances [Ohms];
%   Vs - op-amp supply voltage [V];
%   USE_AA - flag for anti-aliasing using two-point polyBLAMP;
%   PLOT_BLAMP (optional) - flag to plot polyBLAMP values.
% output:
%   Vk - output voltage [V].
function Vk = folding_cell(u, amp, f0, SR, r, Vs, USE_AA, varargin)
    if length(varargin) >= 1
        PLOT_BLAMP = varargin{1};
    else
        PLOT_BLAMP = false;
    end
    g = r(2)*r(3)/(r(1)*r(2) + r(2)*r(3) + r(1)*r(2));
    Nf = length(u);
    Np = floor(Nf/SR*f0);

    thres = r(1)*Vs/r(2);
    cond = (abs(u) > thres);
    Vk = zeros(size(u));
    if amp > thres
        if USE_AA == true
            % inverse clipper
            Vk = cond.*u + (1-cond).*sign(u)*thres;

            % find clipping points
            t1 = asin(thres/amp) / (2*pi*f0);
            t2 = 1/(2*f0) - t1;
            t3 = 1/(2*f0) + t1;
            t4 = 1/f0 - t1;
            clp_t = [t1, t2, t3, t4] + (0:Np-1).'*(1/f0);
            clp_t = reshape(clp_t.', [], 1);

            % fit BLAMP samples to clipping points
            clp_n = floor(clp_t*SR);
            clp_t(clp_n == 0 | clp_n >= length(u)-1) = [];
            clp_n(clp_n == 0 | clp_n >= length(u)-1) = [];
            d = (clp_t - clp_n/SR)*SR;
            mu = abs(2*pi*f0*amp*cos(2*pi*f0*clp_t)/SR);
            pol = sign(amp*sin(2*pi*f0*clp_t));
            res = [(-d.^3/6+d.^2/2-d/2+1/6), (d.^3/6)];
            
            % plot
            if PLOT_BLAMP == true
                % calculate variables for plots
                k = 1/SR;
                tvec = 0:k:10/f0;
                N = nnz(clp_t < max(tvec));
                
                % waveform
                figure;
                subplot(2,1,1)
                hold on;
                plot(tvec, u(1:length(tvec)), 'b');
                plot(tvec, Vk(1:length(tvec)), 'ko-');
                xline(clp_t(1:N), '--');
                xline(clp_t(1:N)+(1-d(1:N))/SR, 'r--');
                xline(clp_t(1:N)-d(1:N)/SR, 'r--');
                yline(thres, 'b--')
                yline(-thres, 'b--')
                xlim([0 1/f0]);
                xlabel('Time [sec]', 'Interpreter', 'latex');
                ylabel('Amplitude [V]', 'Interpreter', 'latex');
                title('Inverse clipper', 'Interpreter', 'latex')
                legend({'Input', 'Clipped signal', 'Clipping points'}, 'Interpreter', 'latex')
                
                % BLAMP samples
                subplot(2,1,2)
                hold on;
                blampl = mu.*res(:,1).*pol;
                blampr = mu.*res(:,2).*pol;
                plot(clp_n(1:N)/SR, blampl(1:N), 'kx');
                plot((clp_n(1:N)+1)/SR, blampr(1:N), 'kx');
                xlim([0 1/f0]);
                ylim([-max(abs([blampl; blampr])), max(abs([blampl; blampr]))]);
                xlabel('Time [sec]', 'Interpreter', 'latex');
                ylabel('Amplitude', 'Interpreter', 'latex');
                title('PolyBLAMP samples', 'Interpreter', 'latex')
            end

            % compensate residual in clipped signal
            clp_n = clp_n+1;
            Vk(clp_n)   = Vk(clp_n)   + mu.*res(:,1).*pol;
            Vk(clp_n+1) = Vk(clp_n+1) + mu.*res(:,2).*pol;

            % output voltage
            Vk = g*(Vk - sign(Vk)*thres);
        else
            % output voltage
            Vk = cond.*(g*(u - sign(u)*thres));
        end
    end
end

% find all peak indexes in an array
% input:
%   x - 1d array;
%   l - number of neighboring values to check;
%   threshold - threshold for a minimum acceptable value of a peak (optional).
% output:
%   peaks_idx - array with peaks indexes.
function peaks_idx = find_peaks(x, l, threshold)
    if nargin < 3
        threshold = -inf;
    end
    k = 1;
    peaks_idx = zeros(length(x),1);
    for i = 1:length(x)
        if check_peak(x, i, l, threshold) == true
            peaks_idx(k) = i;
            k = k+1;
        end
    end
    peaks_idx = peaks_idx(1:k-1);
end

% check if a value in the array is a peak
% input:
%   x - 1d array;
%   idx - value index;
%   l - number of neighboring values to check;
%   threshold - threshold for a minimum acceptable value of a peak (optional).
% output:
%   out - true if value is a peak, false otherwise.
function out = check_peak(x, idx, l, threshold)
    out = true;
    % check threshold
    if nargin < 3
        threshold = -inf;
    end
    if x(idx) <= threshold
        out = false;
        return
    end
    % check bounds
    l_2 = floor(l/2);
    if idx <= l_2 || idx > length(x) - l_2
        out = false;
        return
    end
    % check neighbours
    for j = -l_2:l_2
        if j == 0
            continue;
        end
        if x(idx) <= x(idx+j)
            out = false;
            break
        end
    end
end