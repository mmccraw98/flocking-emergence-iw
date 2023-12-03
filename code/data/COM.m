% clc; clear; close all;
% Data_A = readtable("A_nodes.csv");
% time_A = Data_A.time;
% x_A = Data_A.x;
% y_A = Data_A.y;
% z_A = Data_A.z;
% unique_time_A = unique(time_A);
% x_com_A = zeros(length(unique_time_A));
% y_com_A = zeros(length(unique_time_A));
% %z_com = zeros(l);
% index_A = zeros(length(unique_time_A));
% for ii=1:length(unique_time_A);
%     for jj=1:length(time_A)
%         if unique_time_A(ii) == time_A(jj)
%             x_com_A(ii) = x_com_A(ii) + x_A(jj);
%             y_com_A(ii) = y_com_A(ii) + y_A(jj);
%             %z_com(ii) = z_com(ii) + z(jj);
%             index_A(ii) = index_A(ii) + 1;
%         end
%     end
%     x_com_A(ii) = x_com_A(ii)/index_A(ii);
%     y_com_A(ii) = y_com_A(ii)/index_A(ii);
% end
% r_com_A = sqrt(x_com_A.^2 + y_com_A.^2);
% 
% figure(11); hold on;
% r_com_A_=r_com_A-mean(r_com_A);
% four_A = fft(r_com_A_)/numel(time_A);
% four_A(1)=0;
% fourshift_A=fftshift(four_A);
% freq_A = (-length(fourshift_A)/2:length(fourshift_A)/2 -1)*1/length(fourshift_A);
% plot(freq_A,abs(fourshift_A));

% figure(1); hold on;
% plot(x_com_A, y_com_A);
% 
% figure(2); hold on;
% plot(unique_time_A, r_com_A);
% xlabel('time');
% ylabel('COM positions');
% title('For Data set A');
%%
clc; clear; close all;
Data_B = readtable("B_nodes.csv");
time_B = Data_B.time;
x_B = Data_B.x;
y_B = Data_B.y;
unique_time_B = unique(time_B);
x_com_B = zeros(length(unique_time_B));
y_com_B = zeros(length(unique_time_B));
index_B = zeros(length(unique_time_B));
for ii=1:length(unique_time_B)
    for jj=1:length(time_B)
        if unique_time_B(ii) == time_B(jj)
            x_com_B(ii) = x_com_B(ii) + x_B(jj);
            y_com_B(ii) = y_com_B(ii) + y_B(jj);
            index_B(ii) = index_B(ii) + 1;
        end
    end
    x_com_B(ii) = x_com_B(ii)/index_B(ii);
    y_com_B(ii) = y_com_B(ii)/index_B(ii);
end
r_com_B = sqrt(x_com_B.^2 + y_com_B.^2);


figure(10); hold on;
r_com_B_=r_com_B-mean(r_com_B);
four_B = fft(r_com_B_)/numel(time_B);
four_B(1)=0;
fourshift_B=fftshift(four_B);
freq_B = (-length(fourshift_B)/2:length(fourshift_B)/2 -1)*1/length(fourshift_B);
plot(freq_B,abs(fourshift_B));

figure(3); hold on;
plot(x_com_B, y_com_B);

figure(4); hold on;
plot(unique_time_B, r_com_B);
xlabel('time');
ylabel('COM positions');
title('For Data set B');
% %%
% clc; clear; close all;
% Data_C = readtable("C_nodes.csv");
% time_C = Data_C.time;
% x_C = Data_C.x;
% y_C = Data_C.y;
% unique_time_C = unique(time_C);
% x_com_C = zeros(length(unique_time_C));
% y_com_C = zeros(length(unique_time_C));
% index_C = zeros(length(unique_time_C));
% for ii=1:length(unique_time_C)
%     for jj=1:length(time_C)
%         if unique_time_C(ii) == time_C(jj)
%             x_com_C(ii) = x_com_C(ii) + x_C(jj);
%             y_com_C(ii) = y_com_C(ii) + y_C(jj);
%             index_C(ii) = index_C(ii) + 1;
%         end
%     end
%     x_com_C(ii) = x_com_C(ii)/index_C(ii);
%     y_com_C(ii) = y_com_C(ii)/index_C(ii);
% end
% r_com_C = (x_com_C.^2 + y_com_C.^2);
%
% figure(12); hold on;
% r_com_C_=r_com_C-mean(r_com_C);
% four_C = fft(r_com_C_)/numel(time_C);
% four_C(1)=0;
% fourshift_C=fftshift(four_C);
% freq_C = (-length(fourshift_C)/2:length(fourshift_C)/2 -1)*1/length(fourshift_C);
% plot(freq_C,abs(fourshift_C));
%
% figure(5); hold on;
% plot(x_com_C, y_com_C);
% 
% figure(6); hold on;
% plot(unique_time_C, r_com_C);
% xlabel('time');
% ylabel('COM positions');
% title('For Data set C');
