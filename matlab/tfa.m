data = readmatrix('C:\MY\matlab\SSVEP\seesion2\data\50hz\ear_11hz.csv');
% data = readmatrix('C:\MY\matlab\SSVEP\seesion2\data\50hz\head_.csv');
% data = data(:,1);

dataR = reshape(data, 1, []);
srate = 100;



time = (0 : length(data)/100 * srate)/srate;

% %% sanity check
% data = sin(2*pi*5*time(1:end-1))+ sin(2*pi*11*time(1:end-1));
% dataR = reshape(data, 1, []);
% %%

time = time - mean(time);

min_freq = 4;
max_freq = 15;
num_freq = 36;

frex = linspace(min_freq, max_freq, num_freq); 

tf = zeros(num_freq, length(data));

ndata = length(dataR);
nkern = length(time);
nConv = ndata + nkern - 2;
halfK = floor(nkern/2);

dataX = fft(dataR,nConv);

for fi=1:num_freq
    
    
    cmw  = exp(1i*3.7*pi*frex(fi)*time) .*...
            exp( -2.5*log(2)*time.^2 / .3^2 );
%     % 这个参数可以看到清晰的频率响应    
%     cmw  = exp(1i*4*pi*frex(fi)*time) .*...
%         exp( -2*log(2)*time.^2 / .3^2 );        
        
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
contourf(time(1:end-1),frex,tf,40,'linecolor','none')

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
























