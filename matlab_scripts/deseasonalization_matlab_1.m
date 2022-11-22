clc;
clear;
%==== Import data ==========
filename = 'hourly_series_2014.csv'; % If in NYC
%filename = 'hourly_series_2015_bris.csv'; If in SEQ
%filename = 'hourly_series_2009_den.csv'; If in CPH
opts = detectImportOptions(filename);
opts.VariableTypes{1} = 'char';
% opts = setvaropts(opts,'time', ...
%                        'DatetimeFormat','yyyy-MM-dd HH:mm:ss');
A = readtable(filename, opts);
format longG

%===== Estimation ====
t_raw=datenum(A.time,'yyyy-mm-dd HH:MM:SS' );


[~,id] = unique(t_raw);
sl_raw=A.wl;
cnstit='auto';
lat = 40.42; % If in NYC
% lat =-24.76; If in SEQ
% lat= 56.10; If in CPH
coef = ut_solv ( t_raw(id), sl_raw(id),[],lat, cnstit, 'NoTrend'); % we use 'NoTrend' since data is already detrended.

% 	coef = UT_SOLV ( t_raw, sl_raw, [], lat, cnstit , {options} ); 
%   [ sl_fit, ~ ] = UT_RECONSTR ( t_fit, coef , {options} ); 

%======= Prediction ====
filename = 'hourly_series_bris.csv';
opts = detectImportOptions(filename);
opts.VariableTypes{1} = 'char';
B = readtable(filename, opts);
t_fit=datenum(B.time,'yyyy-mm-dd HH:MM:SS' );
[ sl_fit, ~ ] = ut_reconstr ( t_fit, coef ); % sl_fit is the tide component predicted for data set 'hourly_series.csv'.
B.wlh = sl_fit;
writetable(B,'hourly_series_predict.csv','Delimiter',','); % if in NYC
%writetable(B,'hourly_series_predict_bris.csv','Delimiter',','); %if in SEQ
%writetable(B,'hourly_series_predict_den.csv','Delimiter',','); %if in CPH 


plot(B.wlh)

