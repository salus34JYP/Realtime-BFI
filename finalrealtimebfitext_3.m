%%
clear 
clc
close all

%% pixel location
%                   b     g     o     r     s     y     p
pixel_location = [ 469   515   473   512   535   449   492;
                   489   535   493   532   555   469   512;
                   566   638   640   564   599   605   602;
                   586   658   660   584   619   625   622;
];


sz = size(pixel_location(:,:));
load('background_array.mat');
% load('number_box.mat');
%%
frame_label = zeros(1,5400,'double');
final_BFI = zeros(7,5400,'double');
final_BFI_std = zeros(7,5400,'double');
background_array = zeros(1024,1280,'double');
baseline_array = zeros(1024,1280,'double');
frame_array = zeros(1024,1280,'double');
BSsumimage = zeros(1024,1280,'double');

%%
% prompt = 'What is the frame number? ';
% totalframe = input(prompt);

% RT_READTIFF.M
%
% Basler realtime tiff reader
% make sure the data directory doesn't have any previous tiff data
% before starting the acquisition
prompt = 'What is the frame number? ';
totalframe = input(prompt);
DataDir = pwd;
DataDir_2 = strcat(DataDir,'\');
% imds = imageDatastore(DataDir);
%% acquire header of tiff filenames

while 1
    fn0 = dir([DataDir_2 '*0000.tiff']);
    if size(fn0,1) && (fn0(1).bytes ~= 0)
        Header = extractBefore(fn0.name,'0000');
        break;
    else
       pause(0.001)
    end
end
disp('tiff file detected. Starting the main loop ...')

%% main loop

N = totalframe; % total number of frames
tic
for i = 0:N-1
    str4 = sprintf('%04.0f',i);
    filename = [Header str4 '.tiff'];

    while 1
        if isfile(filename)
%             disp(filename)
             pause(0.035);
          

            tiff_image = imread(filename);
            tiff_image = cast(tiff_image,'double');
            tiff_image = tiff_image - background_array;
            frame_label(1,i+1) = i+1;

            for label = 1 : sz(1,2)
                BFI_box = zeros(9,1);

                for count = 1 : 9
                    y = pixel_location(1,label) : 7 : pixel_location(1,label) + 20;
                    x = pixel_location(3,label) : 7 : pixel_location(3,label) + 20;
                    [Y,X] = meshgrid(y,x);
                    y_n = numel(Y);
                    X_m = numel(X);
                    Y_pixel = reshape(Y, [y_n ,1]);
                    X_pixel = reshape(X, [X_m ,1]);
                    MEAN = mean(tiff_image(Y_pixel(count):Y_pixel(count)+6,X_pixel(count):X_pixel(count)+6),"all");
                    STD = std(tiff_image(Y_pixel(count):Y_pixel(count)+6,X_pixel(count):X_pixel(count)+6),1,'all');
                    K = STD/MEAN;
                    BFI = 1/(K)^2;
                    BFI_box(count,1) = BFI;
                end

                BFI_box_mean = mean(BFI_box(:,1));
                BFI_box_std = std(BFI_box(:,1));

                final_BFI(label,i+1) = BFI_box_mean;
                final_BFI_std(label,i+1)= BFI_box_std;
           end
   
%          hold on
         plot(frame_label,final_BFI(7,:),'.',...
        'Color','#7E2F8E')
         xlabel("times (s)")
         ylabel("BFI (A. U)")
         xticks([0 300  600 900 1200 1500 1800 2100 2400 2700 3000 3300 3600 3900 4200 4500 4800 5100 5400])
         xticklabels({'0','30','60','90','120','150','180','210','240','270','300','330','360','390','420','450','480','510','540'})
         title("Blood Flow Index vs times")

         text(350, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.1, 'baseline','Color','#000000','FontSize',15)
         text(1350, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.1, 'insert','Color','#000000','FontSize',15)
         text(1250, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.05, 'acupuncture','Color','#000000','FontSize',12)
         text(1850, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.1, 'stimulate','Color','#000000','FontSize',15)
         text(2500, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.1, '1st rest','Color','#000000','FontSize',15)
         text(3200, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.1, 'laser','Color','#000000','FontSize',15)
         text(3030, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.05, 'acupuncture','Color','#000000','FontSize',12)
         text(3700, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.1, '2nd rest','Color','#000000','FontSize',15)
         text(4700, max(final_BFI(7,:))+ max(final_BFI(7,:))*0.1, 'release','Color','#000000','FontSize',15)

         xline([1200 1800 2400 3000 3600 4200],'--','LineWidth',2,'Color','#3832a8')
%          hold off
%          text(3300,100,'BFI =','FontSize',25,'FontWeight','bold');
%          if number_box(i+1) == fix((i+1)/10)
%              round_bfi = round(final_BFI(label,number_box(i+1)+1),1);
%              txt = sprintf('%.0f', round_bfi);
%              final_bfi_text = num2str(txt);
%              text(3550, max(final_BFI(7,:))*1.25, final_bfi_text,'FontSize',25,'FontWeight','bold');
%              
%          end
         
         axis([-300 5700 0 max(final_BFI(7,:))+ max(final_BFI(7,:))*0.2])
         drawnow;
         

            break
        else
            disp(filename)
            pause(0.001);
        end
    end
end
toc
%%
baseline_mean_BFI = mean(final_BFI(:,1:1200),2);
normalized_BFI = final_BFI./baseline_mean_BFI;
fulltrans = transpose(normalized_BFI);
save('0903BFI','fulltrans')