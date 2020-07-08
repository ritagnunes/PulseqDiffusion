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
addpath(genpath('./jMRI_Publication_Functions'));

clear all; clc; close all

testResults=runtests('.')