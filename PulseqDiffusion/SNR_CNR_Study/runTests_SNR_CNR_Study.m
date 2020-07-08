%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior Tecnico - Universidade de Lisboa
%
% name: runTests_SNR_CNR_Study
% type: Script
% description:
% It runs all unit test for every function on this folder.

%% Initialize
addpath(genpath('./SNR_CNR_Study'));

clear all; clc; close all

testResults=runtests('.')