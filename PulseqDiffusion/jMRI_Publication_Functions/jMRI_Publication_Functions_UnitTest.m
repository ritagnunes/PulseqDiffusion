%% Read Me
% 
% by:
% T.T. Fernandes, September 2019
% LarSys - Instituto Superior Tecnico - Universidade de Lisboa
%
% name: jMRI_Publication_Functions_UnitTest
% type: Script
% description:
% It runs all unit test for every function on this folder.

%% Initialize
addpath(genpath('./SjMRI_Publication_Function'));

clear all; clc; close all

%% Test 1 - TE_TimeTest
fprintf('--- 1 - Test ---\n');
TE_timeTest
fprintf('\n\n');

%% Test 2 - aux_deltaTest
fprintf('--- 2 - Test ---\n');
aux_deltaTest
fprintf('\n\n');

%% Test 3 - deltaTest
fprintf('--- 3 - Test ---\n');
deltaTest
fprintf('\n\n');

%% Test 4 - DWSignalTest
fprintf('--- 4 - Test ---\n');
DWSignalTest
fprintf('\n\n');

%% Test 5 - Ernstangle_dTest
fprintf('--- 5 - Test ---\n');
Ernstangle_dTest
fprintf('\n\n');


%% Test 6 - LongSS_SignalTest
fprintf('--- 6 - Test ---\n');
LongSS_SignalTest
fprintf('\n\n');


