close all
clear
clc

clc
M = 128; % signal dimension
N = M; % nr of atoms in the basis (this will be different when using redundant set of generators)

DCT = nan(M, N); % matrix containing the standard basis (a kronecker delta in each column)
D = nan(M, N); % matrix containing the DCT basis (a DCT function in each column)

%% generate 1D DCT basis (DCT type II)
disp('generating 1-D DCT basis')

for k = 1 : M
    % take the formula from slides and remember to normalize. Each atom
    % goes in a column of DCT matrix
    % DCT(:, k) = 

    % display the atom
    figure(51), plot(DCT(:, k), 'b-')
    title(['element: ', num2str(k) ' from DCT basis']);

end
% check orthogonality
% is_DCT_orth =

% display the basis in the matrix
figure(3), imagesc(DCT), title('DCT basis, atoms in the column'), axis equal, axis tight


%% generate 1D DCT basis using the function idct which is the inverse dct transform
% stack this in the matrix D
for k = 1 : M
    % define the atom
    a = zeros(1, M);
    a(k) = 1;
    D(:, k) = idct(a)';
end
% display the basis in the matrix
figure(3), imagesc(D), title('DCT basis, atoms in the column'), axis equal, axis tight

% check that D and DCT defined above coincide
% is_DCT_equal_D = 

% check orthogonality
% is_D_orth = 

%% Analysis: compute the representation of an input ECG signal

% load few ECG signals to be processed
temp = load('data/croppedECG.mat'); 
nBeats = 10;
S = temp.S(:, 1 : nBeats);
coeff_D = nan(M, nBeats); % nitialize the matrix of representations of S w.r.t. D

% display the signal and the representation coefficients
for ii = 1 : size(S, 2)
    
    % compute the representation of each beat w.r.t. the two basis 
    % coeff_D(:, ii) = 
    
    figure(5)
    subplot(2,1,1)
    plot(S(:, ii), 'r'), title(['original beat ', num2str(ii)]);
    subplot(2,1,2)
    plot(coeff_D(:, ii), 'm'), title('coefficients w.r.t. DCT basis');
    pause(0.2)
    set(gcf, 'Position', [680   2   372   815])
    %     disp('hit a key')
    %     pause()
end

%% Synthesis: reconstruct all the ECG signals from their representations

% reconstruct the two signals (express them w.r.t. the standar dbasis)
% S_hat_D = 

% check if there is perfect reconstruction
% trivial because S_hat_D = D * D' * S and D * D' = eye(M)  since D is orthonormal basis (the same applies to C)
disp('perfect reconstruction?')
% is_reconstruction_perfect = 

%% Add noise to ECG data and inspect the representations

% add AWGN noise to ECG signals
nBeats = 10;    
temp = load('data/croppedECG.mat'); % load few ECG signals to be processed
S0 = temp.S(:, 1 : nBeats);

sigma_noise = 0.1;
S = S0 + sigma_noise * randn(size(S0));

coeff_D = nan(M, nBeats); % nitialize the matrix of representations of S w.r.t. D

% display the signal and the representation coefficients
for ii = 1 : size(S, 2)
    
    % compute the representation of each beat w.r.t. the two basis 
    % coeff_D(:, ii) 
    
    figure(5)
    subplot(2,1,1)
    plot(S(:, ii), 'r'), title(['original beat nr ', num2str(ii)]);
    subplot(2,1,2)
    plot(coeff_D(:, ii), 'g'), title('coefficients w.r.t. DCT basis');
    set(gcf, 'Position',  [680   2   372   815])
    drawnow
    pause(0.05)
    
    % disp('hit a key')
    % pause()
end

%% Hard thresholding
% noise affects all the coefficients of our transformation
% keep only few coefficients having largest magnitude

% add AWGN noise to ECG signals
nBeats = 10;
temp = load('data/croppedECG.mat'); % load few ECG signals to be processed
S0 = temp.S(:, 1 : nBeats);
sigma_noise = 0.1;

% add noise
S = S0 + sigma_noise * randn(size(S0));

L = 21; % sparsity level (try different values) size(S, 1)

for ii = 1 : 2 %nBeats
    
    % origSignal = 
    % noisySignal = 
    
    % transform each signal separately (analysis) 
    % coeff = 
    
    % keep only the L largest coefficients (absolute value)
    % [...]
    % coeff_HT = ;
    
    % invert the transformation
    % s_hat = 
    
    LN_WDT = 3;
    MRK_SZ = 10;
    
    % show the signals
    figure(5)
    subplot(1,2,1)
    plot(noisySignal, 'r.', 'MarkerSize', MRK_SZ), 
    hold on
    plot(origSignal, 'b--', 'LineWidth', LN_WDT),  
    plot(s_hat, 'k-', 'LineWidth', LN_WDT),
    title(['original beat nr ', num2str(ii)]);
    hold off
    legend('noisy', 'original', 'hard-thresholded')
    
    subplot(1,2,2)
    plot(coeff, 'r.', 'MarkerSize', MRK_SZ), 
    hold on
    % coefficients of the noise free signal
    plot(D' * origSignal, 'b--', 'LineWidth', LN_WDT + 1),  
    stem(find((coeff_HT ~= 0)), coeff_HT(coeff_HT ~= 0), 'k-', 'LineWidth', LN_WDT, 'MarkerSize', MRK_SZ),  
    title('DCT Coefficients')
    hold off
    legend('noisy', 'original', 'hard-thresholded')

    set(gcf, 'Position', [423        97        1206         636])

    pause(0.05)
    % disp('hit a key')
    % pause()
end



