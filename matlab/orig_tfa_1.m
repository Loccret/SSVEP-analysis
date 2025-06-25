% 
% S1_5Hz1 = readmatrix('C:\MY\matlab\SSVEP\orig\S2\orig\S2_9Hz2_data.csv');
% 
% temp = S1_5Hz1(2:end,50:55);
% % temp = S1_5Hz1(2:end,2:9);
% band            = [4, 15.]; 
% srate = 2000;
% 
% temp_filtered = bandpass(temp, band, srate);
% temp_filtered = temp_filtered(1:20:end - 6000,:);

%% read haba
temp_filtered = squeeze(band_ele_2(2001:20:3000,12,4,1,8:end));


%%
data = temp_filtered;

dataR = reshape(data, 1, []);
srate = 100;



time = (0 : length(data)/100 * srate)/srate;

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
nConv = ndata + nkern - 2;
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
contourf(time(1:end-1),frex,tf,40,'linecolor','none')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colormap(jet)
colorbar




