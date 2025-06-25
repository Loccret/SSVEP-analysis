% version 2: change to flexible morlet wave freq
clear

data = readmatrix('C:\MY\matlab\SSVEP\seesion2\data\50hz\ear_11hz.csv');
% data = readmatrix('C:\MY\matlab\SSVEP\seesion2\data\50hz\ear_.csv');
% data = data(:,1);

dataR = reshape(data, 1, []);
srate = 100;



time = (0 : length(data)/100 * srate - 1)/srate;

% %% sanity check
% data = sin(2*pi*5*time(1:end-1))+ sin(2*pi*11*time(1:end-1));
% dataR = reshape(data, 1, []);
% %%

time = time - mean(time);

min_freq = 2;
max_freq = 26;
num_freq = 60;

% max_freq = 50;
% num_freq = 50;

frex = linspace(min_freq, max_freq, num_freq); 

tf = zeros(num_freq, length(data));

ndata = length(dataR);
nkern = length(time);
nConv = ndata + nkern - 1;
halfK = floor(nkern/2);

dataX = fft(dataR,nConv);


% range_cycles = [ 5 15 ];
% nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);

for fi=1:num_freq
    
%     s = nCycles(fi)/(2*pi*frex(fi));
%     cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
        
    cmw  = exp(1i*2*pi*frex(fi)*time) .*...
        exp( -4*log(2)*time.^2 / 0.6^2 );
        
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % the rest of convolution
    as = ifft( dataX.*cmwX ); %analytic signal
    as = as(halfK+1:end-halfK+1);
    as = reshape(as,size(data));
    
    % extract power
    aspow = abs(as).^2;
    
    % average over trials and put in matrix
    tf(fi,:) = mean(aspow,2);
end


%%% and plot!
% plot(sum(tf,2))
figure(203), clf
contourf(time,frex,tf,40,'linecolor','none')

% figure(14), clf
% contourf(tf)
% figure(15),clf
% contourf(tf,40,'linecolor','none')
% figure(16),clf
% contourf(tf,20,'linecolor','none')

% set(gca,'clim',[0 1]*10000,'xlim',[-.1 1.4])

% set(gca,'clim',[0. 0.3])
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(jet)
colorbar

% %% Morlet wavelet
% time = -1:1/1000:1;
% fwhm = 0.8;
% gaus_win = exp((-4 * log(2) * time.^2)/fwhm^2);
% plot(time, gaus_win);
%% deal with original data


% S1_5Hz1 = readmatrix('C:\MY\matlab\SSVEP\orig\S1\orig\S1_5Hz1_data.csv');
% 
% temp = S1_5Hz1(2:end,50:55);
% band            = [4, 15.]; 
% srate = 2000;
% 
% temp_filtered = bandpass(temp, band, srate);
% temp_filtered = temp_filtered(1:20:end,:);
% 
% figure(1), clf
% tt = 0;
% for i = 1:16
%     tt = 1+ tt;
%     subplot(4, 4, tt);
%     plot(temp_filtered(1+(i-1)*100 : i*100, 1))
%     subtitle(num2str(tt))
% end
% 
% 
% temp_filtered_2 = band_ele_2(1:20:end,1,1,1,1);
% 
% figure(2), clf
% tt = 0;
% for i = 1:9
%     tt = 1+ tt;
%     subplot(3, 3, tt);
%     plot(temp_filtered_2(1+(i-1)*100 : i*100))
%     subtitle(num2str(tt))
% end


%%


























