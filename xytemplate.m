clc
clear all
fs=360;
Fs=360;
gr=0;
[record,sig]=rdsamp('116',[],2000);
ecg= record(:, 1);
[annotation,type] = rdann('116','atr',[],2000);
type=type(2:end-1); 
annotation=annotation(2:end-1);
win_size = 0.5; 
win_len = round(win_size * fs);
ecg_median = medfilt1(ecg, win_len);
ecg_corrected = ecg - ecg_median;
Fs = 360; 
fc = 15; 
[b,a] = butter(2,fc/(Fs/2),'low');

filtered_signal = filter(b,a,ecg_corrected);

notch_freq = 60;
[b_notch,a_notch] = butter(2,[notch_freq-1, notch_freq+1]/(Fs/2),'stop');

ecg = filter(b_notch,a_notch,filtered_signal);

if ~isvector(ecg)
  error('ecg must be a row or column vector');
end
ecg = ecg(:); % vectorize
flag = 1;

Signal_mean = [];
Signal_mean_index = [];
Thresh_Signal_mean = 0;
Thresh_Noise_mean = 0;

Noise_mean = [];
Noise_mean_index = [];

Signal_filter = [];
Signal_filter_index = [];
Thresh_Signal_filter = 0;
Thresh_Noise_filter = 0;

Signal_temp_mean = 0;
Noise_temp_mean = 0;
Signal_temp_filter = 0;
Noise_temp_filter = 0;

Signal_mean_buf = [];
Noise_mean_buf = [];

Signal_filter_buf = [];
Noise_filter_buf = [];
Thresh_Signal_mean_buf = [];
Thresh_Signal_filter_buf = [];

regular_RR_interval = 0;

mean_interval = 0;

RR_interval =  0;

R_Flag = 0;

if Fs == 360
    % H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
    b = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a = [1 -2 1];
    h_l = filter(b,a, [1,zeros(1,12)]);
    ecg_l = conv(ecg, h_l);
    ecg_l = ecg_l / max(abs(ecg_l));
    % H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
    b = [-1 zeros(1,15) 32 -32 zeros(1,14) 1];
    a = [1 -1];
    h_h = filter(b,a, [1, zeros(1,32)]);
    ecg_filter = conv(ecg_l, h_h);
    ecg_filter = ecg_filter / max(abs(ecg_filter));

else
    wp = [5 15] * 2 / Fs;
    [B,A] = butter(3, wp);
    ecg_filter = filtfilt(B,A, ecg);
    ecg_filter = ecg_filter / max(abs(ecg_filter));
end
%% H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
    h_d = [-1 -2 0 2 1];
    ecg_deri = conv(h_d, ecg_filter);
    ecg_deri = ecg_deri / max(abs(ecg_deri));
    ecg_square = ecg_deri .^ 2; 
    Window_width = 0.15 * Fs;
    ecg_mean = conv(ecg_square, ones(1, Window_width) / Window_width);

[peaks, locs] = findpeaks(ecg_mean, 'MINPEAKDISTANCE', round(0.2 * Fs));

Thresh_Signal_mean = max(ecg_mean(1:2*Fs)) / 3;
Thresh_Noise_mean = mean(ecg_mean(1:2*Fs)) / 2;
Signal_temp_mean = Thresh_Signal_mean;
Noise_temp_mean = Thresh_Noise_mean;

Thresh_Signal_filter = max(ecg_filter(1:2*Fs)) / 3;
Thresh_Noise_filter = mean(ecg_filter(1:2*Fs)) / 2;
Signal_temp_filter = Thresh_Signal_filter;
Noise_temp_filter = Thresh_Noise_filter;

for i = 1 : length(peaks)   
    if locs(i) - round(0.15 * Fs) >= 1 && locs(i) <= length(ecg_filter)
        [peak_filter, index_filter] = max( ecg_filter(locs(i) - round(0.15 * Fs) : locs(i)));
        index_filter = locs(i) - round(0.15 * Fs) + index_filter - 1;
    else
        if i == 1 
            [peak_filter, index_filter] = max( ecg_filter(1 : locs(i)));
        else 
            [peak_filter, index_filter] = max( ecg_filter(locs(i) - round(0.15 * Fs) : end));
            index_filter = locs(i) - round(0.15 * Fs) + index_filter - 1 ;
        end
    end
    gr=0;
    
    if length(Signal_mean) >= 9
        RR_diff = diff(Signal_mean_index(end - 8:end));
        mean_interval = mean(RR_diff);
        temp_interval = Signal_mean_index(end) - Signal_mean_index(end - 1); 
        if(temp_interval >= 0.92 * mean_interval || temp_interval <= 1.16 * mean_interval)
            Thresh_Signal_mean = Thresh_Signal_mean / 2;
            Thresh_Signal_filter = Thresh_Signal_filter / 2;
        else
            regular_RR_interval = mean_interval;
        end
    end
    if regular_RR_interval
        RR_interval = regular_RR_interval;
    elseif mean_interval
        RR_interval = mean_interval;
    else
        RR_interval = 0;
    end
  
    if RR_interval
        
        if locs(i) - Signal_mean_index(end) >= 1.66 * RR_interval
            [pk_temp, pk_temp_ind] = max( ecg_mean(Signal_mean_index(end) + round(0.2 * Fs) : locs(i) - round(0.2 * Fs)));
            pk_temp_ind = Signal_mean_index(end) + round(0.2 * Fs) + pk_temp_ind - 1;
            
            if pk_temp >= Thresh_Noise_mean
                Signal_mean = [Signal_mean pk_temp];
                Signal_mean_index = [Signal_mean_index pk_temp_ind];       
                
                if pk_temp_ind <= length(ecg_filter)
                    [pk_filter_temp, pk_filter_temp_ind] = max(ecg_filter(pk_temp_ind - round(0.15 * Fs) :  pk_temp_ind));
                else
                    [pk_filter_temp, pk_filter_temp_ind] = max(ecg_filter(pk_temp_ind - round(0.15 * Fs) :  length(ecg_filter)));
                end        
                pk_filter_temp_ind = ecg_filter(pk_temp_ind - round(0.15 * Fs) + pk_filter_temp_ind - 1);
                   
                if pk_filter_temp >= Thresh_Noise_filter
                    Signal_filter = [Signal_filter pk_filter_temp];
                    Signal_filter_index = [Signal_filter_index pk_filter_temp_ind];      
                    Signal_temp_filter = 0.25 * pk_filter_temp + 0.75 * Signal_temp_filter;
                end
                Signal_temp_mean = 0.25 * pk_temp + 0.75 * Signal_temp_mean;
            end
        end
    end
    
    if peaks(i) >= Thresh_Signal_mean
        
        if(length(Signal_mean) >= 3)        
            
            if locs(i) - Signal_mean_index(end) < round(0.36 * Fs)
                Slope_1 = mean(diff(ecg_mean(locs(i) - round(0.075 * Fs) : locs(i))));
                Slope_2 = mean(diff(ecg_mean(Signal_mean_index(end) - round(0.075 * Fs) : Signal_mean_index(end))));              
                
                if Slope_1 <= Slope_2 / 2              
                    Noise_mean = [Noise_mean peaks(i)];
                    Noise_mean_index = [Noise_mean_index locs(i)];
                    Noise_temp_mean = 0.25 * peaks(i) + 0.75 * Noise_temp_mean;
                    Noise_temp_filter = 0.25 * peak_filter + 0.75 * Noise_temp_filter;
                    R_Flag = 1;              
                else
                    R_Flag = 0;
                end
            end
        end  
        if R_Flag == 0
            Signal_mean = [Signal_mean peaks(i)];
            Signal_mean_index = [Signal_mean_index locs(i)];
            Signal_temp_mean = 0.125 * peaks(i) + 0.875 * Signal_temp_mean;
            if peak_filter >= Thresh_Signal_filter
                Signal_filter = [Signal_filter peak_filter];
                Signal_filter_index = [Signal_filter_index index_filter];
                Signal_temp_filter = 0.125 * peak_filter + 0.875 * Signal_temp_filter;
            end
        end              
    
    elseif peaks(i) < Thresh_Signal_mean && peaks(i) >= Thresh_Noise_mean
        Noise_temp_mean = 0.125 * peaks(i) + 0.875 * Noise_temp_mean;
        Noise_temp_filter = 0.125 * peak_filter + 0.875 * Noise_temp_filter;
    
    else
        Noise_temp_mean = 0.125 * peaks(i) + 0.875 * Noise_temp_mean;
        Noise_temp_filter = 0.125 * peak_filter + 0.875 * Noise_temp_filter;
        Noise_mean = [Noise_mean peaks(i)];
        Noise_mean_index = [Noise_mean_index locs(i)];
    end
    
    Thresh_Signal_mean = Noise_temp_mean + 0.25 * (Signal_temp_mean - Noise_temp_mean);
    Thresh_Noise_mean = Thresh_Signal_mean / 2;
    Thresh_Signal_filter = Noise_temp_filter + 0.25 * (Signal_temp_filter - Noise_temp_filter);
    Thresh_Noise_filter = Thresh_Signal_filter / 2;
    
    Signal_mean_buf = [Signal_mean_buf Signal_temp_mean];
    Noise_mean_buf = [Noise_mean_buf Noise_temp_mean];
    
    Signal_filter_buf = [Signal_filter_buf Signal_temp_filter];
    Noise_filter_buf = [Noise_filter_buf Noise_temp_filter];
    Thresh_Signal_mean_buf = [Thresh_Signal_mean_buf Thresh_Signal_mean];
    Thresh_Signal_filter_buf = [Thresh_Signal_filter_buf Thresh_Signal_filter];
   
    R_Flag = 0;
    
end
r_locs =Signal_filter_index-28;
r_locs = r_locs(2:end-1);

n_templates = 0;
templates = zeros(3,311); 
rr_intervals = [];

t_amplifications = [];
r_amplification_ratios = [];
qrs_durations=[];

[qrs_start, qrs_end] = find_qrs_boundaries(ecg, r_locs, fs);
qrs_durations=qrs_end-qrs_start;

window_start = zeros(1,length(r_locs));
window_end = zeros(1,length(r_locs));
waveform_matrix = zeros(311,length(r_locs));

i = 1;
while i <length(r_locs)
 
    window_start(i) = r_locs(i) - 160;
    window_end(i) = r_locs(i) + 150;
    waveform = ecg(window_start(i):window_end(i));
    rr_interval = diff([r_locs(i), r_locs(i+1)])/fs;

    t_amplification = max(waveform(201:end));
    r_amplification_ratio = max(waveform(1:200))/max(waveform); 
     
    qrs_duration = qrs_durations(i);
    qrs_duration=qrs_duration/fs;
    
    if rr_interval <0.8 && t_amplification > 0 && r_amplification_ratio > 0.3 && qrs_duration<0.12 
        
        if n_templates < 3
            n_templates = n_templates + 1;      
            templates(n_templates,:) = waveform;
        
        else
            ncc_values = zeros(1, size(templates, 1));
            for j = 1:size(templates, 1)
                for k = 1:size(templates, 1)
                x_mean = mean(templates(j,:));
                y_mean =  mean(templates(k,:));
                numerator = sum((templates(j,:) - x_mean) .* (templates(k,:) - y_mean));
                denominator = sqrt(sum((templates(j,:) - x_mean).^2)) .* sqrt(sum((templates(k,:) - y_mean).^2));
                ncc_values(j) = numerator / denominator;
                end
            end
            
            [~, idx] = min(ncc_values);
            
            templates(idx,:) = waveform;          
        end
    end
    
    waveform_matrix(:,i) = waveform;
    rr_intervals(i) = rr_interval;
    t_amplifications(i) = t_amplification;
    r_amplification_ratios(i) = r_amplification_ratio;
    
    i = i + 1;
end

n_templates = size(templates, 1);
ncc_matrix = zeros(n_templates, n_templates);

for i = 1:n_templates
    for j = 1:n_templates     
        
        x_mean = mean(templates(i,:));
        y_mean = mean(templates(j,:));     
        numerator = sum((templates(i,:) - x_mean) .* (templates(j,:) - y_mean));
        denominator = sqrt(sum((templates(i,:) - x_mean).^2) ) .* sqrt( sum((templates(j,:) - y_mean).^2));%分母
        ncc_matrix(i,j) = numerator / denominator;
    end
end
template_width = 311;
figure;
for i = 1:n_templates
    for j = 1:n_templates
        
        text((i-1)*template_width+150, (j-1)*template_width+150, sprintf('%.2f', ncc_matrix(i,j)), 'HorizontalAlignment', 'center');
    end
end
imagesc(ncc_matrix);
title('NCC Matrix');
xlabel('Template Index','FontWeight', 'bold');
ylabel('Template Index','FontWeight', 'bold');
colorbar;
avg_ncc = mean(ncc_matrix, 2);
max_ncc = max(avg_ncc);
min_ncc = min(avg_ncc);
disp(['三个模板节拍的平均最大相关系数：' num2str(max_ncc )]);
disp(['三个模板节拍的平均最小相关系数：' num2str(min_ncc)]);

x_labels = {'Template 1', 'Template 2', 'Template 3', 'Template 4', 'Template 5'};

figure;
bar(avg_ncc);
title('Average NCC Values');
ylabel('Average NCC','FontWeight', 'bold');
set(gca, 'XTickLabel', x_labels);

hold on;
plot(find(avg_ncc == max_ncc), max_ncc, 'r*', 'MarkerSize', 10);
plot(find(avg_ncc == min_ncc), min_ncc, 'g*', 'MarkerSize', 10);
hold off;
legend('Average NCC', 'Maximum NCC', 'Minimum NCC');

for i = 1:size(templates, 1)
  figure
    
   plot(templates(i,:));
   title(sprintf('Template %d', i));
end

template1= templates(1,:);
template2= templates(2,:);
template3= templates(3,:);

Templates = [template1; template2; template3];
new_templates = size(Templates, 1);
wave =zeros(length(r_locs),311);
wave_matrix = zeros(length(r_locs),311);
ncc_values = zeros(length(r_locs), new_templates);

for i = 1:length(r_locs)

    wave_start(i) = r_locs(i) - 160;
    wave_end(i) = r_locs(i) + 150;
    wave = ecg(wave_start(i):wave_end(i));
    wave_matrix(i,:) = wave;
    for j = 1:new_templates
 
        template_mean = mean(Templates(j,:));
        wave_mean = mean(wave);        
   
        fenzi = sum((Templates(j,:) - template_mean) .* (wave_matrix(i,:) - wave_mean));
        fenmu = sqrt(sum((Templates(j,:) - template_mean).^2)) .* sqrt(sum((wave_matrix(i,:) - wave_mean).^2));
        ncc_values(i,j) = fenzi / fenmu;
    end
end

i=1;
frr_intervals = [];
ft_amplifications = [];
fr_amplification_ratios = [];
fwindow_start = zeros(1,length(r_locs));
fwindow_end = zeros(1,length(r_locs));
fwaveform_matrix = zeros(311,length(r_locs));
results =zeros(size(type));
fqrs_durations=[];

[qrs_start, qrs_end] = find_qrs_boundaries(ecg, r_locs, fs);
fqrs_durations=qrs_end-qrs_start;

for i = 1:length(r_locs)-1
  
    if max(ncc_values(i,:)) < 0.8 && min(ncc_values(i,:))<min_ncc
        
    fwindow_start(i) = r_locs(i) - 160;
    fwindow_end(i) = r_locs(i) + 150;
   
    fwaveform = ecg(fwindow_start(i):fwindow_end(i));
    
    frr_interval = diff([r_locs(i), r_locs(i+1)])/fs;

    fqrs_duration = fqrs_durations(i);
    fqrs_duration=fqrs_duration/fs;     
    ft_amplification = max(fwaveform(201:end)); 
    
    if frr_interval >1.2|| ft_amplification <0 ||  min(ncc_values(i,:))<min_ncc || fqrs_duration>0.12
        results(i) = 1;
       disp(['当前第', num2str(i), '个心拍为PVC']);
  
    fwaveform_matrix(:,i) = fwaveform;
    frr_intervals(i) = frr_interval;
   
    ft_amplifications(i) = ft_amplification;
   
    i = i + 1;
    else
        results(i) = 0;
    end  
    
    else
       results(i) = 0;  
    end
end


TP = 0; 
TN = 0; 
FP = 0; 
FN = 0; 
for i = 1:length(type)
    if type(i) == 'V' 
        if results(i) == 1 
            TP = TP + 1;
        else 
            FN = FN + 1;
        end
    else 
        if results(i) == 0 
            TN = TN + 1;
        else 
            FP = FP + 1;
        end
    end
end
accuracy = (TP + TN) / (TP + TN + FP + FN);
sensitivity = TP / (TP + FN);
specificity = TN / (TN + FP);
accuracy_percent = accuracy * 100;
sensitivity_percent = sensitivity * 100;
specificity_percent = specificity * 100;
fprintf('准确度：%f%%\n', accuracy_percent);
fprintf('灵敏度：%f%%\n', sensitivity_percent);
fprintf('特异性：%f%%\n', specificity_percent);

waveform_color = 'b';
figure;
hold on;
for i = 1:size(templates, 1)
    x_shift = (i-1)*template_width;
    if i < size(templates, 1)
        templates(i,end) = templates(i+1,1);
    end
    plot((1:length(templates(i,:)))+x_shift, templates(i,:), 'color', waveform_color);
end
hold off;
xlabel('Sample','FontWeight', 'bold');
ylabel('Amplitude','FontWeight', 'bold');

function [qrs_start, qrs_end] = find_qrs_boundaries(ecg_signal, r_locs, fs)

search_range = 30; 
num_points = 55; 
threshold_factor = 1/5; 
min_slope_points = 5; 
num_r_locs = length(r_locs);
qrs_start = zeros(1, num_r_locs);
qrs_end = zeros(1, num_r_locs);
for i = 1:num_r_locs
    r_loc = r_locs(i);    

    start_search = max(1, r_loc - search_range);
    end_search = min(length(ecg_signal), r_loc + search_range);
    [~, start_idx] = max(diff(ecg_signal(start_search:r_loc)));
    [~, end_idx] = max(diff(ecg_signal(r_loc:end_search)));    
    start_idx = start_idx + start_search - 1;
    end_idx = end_idx + r_loc - 1;   

    search_start = max(1, start_idx - num_points);

     search_end = max(end_idx, end_idx + num_points);
    threshold1 = threshold_factor * max(abs(diff(ecg_signal(search_start:start_idx))));
    threshold2 = threshold_factor *  max(abs(diff(ecg_signal(end_idx:search_end))));
    for j = start_idx:-1:search_start
        if all(diff(ecg_signal(j:j+min_slope_points-1)) < threshold1)
            qrs_start(i) = j;
            break;
        end
    end    
    
    for k= end_idx:1:search_end
        if all(diff(ecg_signal(k-min_slope_points+1:k)) < threshold2)
            qrs_end(i) = k;
            break;
        end
    end    
end
end

 
   