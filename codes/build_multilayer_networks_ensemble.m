%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% script to build the ensemble of equivalent multilayer networks %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace

clear
clc
close all

%% set paths and directories

dir_data = 'D:\Mary\work\Lifespan\Data';
savedir_net = 'D:\Mary\work\Lifespan\Data\Network_bootstrap';


%% load data

load(fullfile(dir_data,'MLNetwork'))
% load a matrix of dimension [N*N*T]: N=114(nodes), L=620(subjects)
% [after the preprocessing we concatenated all the anatomical networks (one
% for each subject) in a 3D matrix from the youngest subject to the oldest] 

loadname(fullfile(dir_data,'Age'))
% vector [L*1] containing ages associated to each network


%%

