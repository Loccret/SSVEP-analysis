% version 2: change to flexible morlet wave freq
% version 3: read all ele data as (8, 50, 2150, 4)
clear

dataORI = readmatrix('C:\MY\matlab\SSVEP\seesion2\data\pre_pro_head_500ms.csv');
data = reshape(dataORI, 4, 4, 540, 50, 6);    % head; (pre-pro, stimulus)

% % % data = reshape(dataORI, 8, 50, 540, 4, 4);  % (..., stimulus, pre-pro)


srate = 100;



% %% sanity check
% data = sin(2*pi*5*time(1:end-1))+ sin(2*pi*11*time(1:end-1));
% dataR = reshape(data, 1, []);
% %%



min_freq = 2 ;
max_freq = 26;
num_freq = 60;

% max_freq = 50;
% num_freq = 50;



% range_cycles = [ 18 30 ];
% nCycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_freq);



figure(1009),clf
tt = 0;
for i=1:4
    for j = 1:4
        temp_data = squeeze(data( i, j, 101, :, :));
        % 200 9Hz的刺激看起来和5Hz那么像？
        dataR = reshape(temp_data, 1, []);
        
        time = (0 : length(temp_data)/100 * srate)/srate;
        time = time - mean(time);
        
        frex = linspace(min_freq, max_freq, num_freq); 
        tf = zeros(num_freq, length(temp_data));

        ndata = length(dataR);
        nkern = length(time);
        nConv = ndata + nkern - 2;
        halfK = floor(nkern/2);
        
        dataX = fft(dataR,nConv);
        
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
            as = reshape(as,size(temp_data));

            % extract power
            aspow = abs(as).^2;

            % average over trials and put in matrix
            tf(fi,:) = mean(aspow,2);
            
        end
        
        
        tt  = tt + 1;
        subplot(4,4,tt);
        contourf(time(1:end-1),frex,tf,40,'linecolor','none')
        xlabel('Time (s)'), ylabel('Frequency (Hz)')
        subtitle(['pre-processing: ', num2str(5+(i-1)*2),'| stimulus: ', num2str(5+(j-1)*2)]);
        colormap(jet)
        colorbar
    end
end


























