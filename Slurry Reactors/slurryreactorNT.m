function[]=slurryreactorNT

%Fischer Tropch 
% CO+2H2 = -CH2- +H20

clear all
clc

T=515 %Kelvins
P=30; %Bars
krxn=1.5e-4; %reaction rate constant (molCO/kgcat s)
