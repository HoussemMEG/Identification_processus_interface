function [M] = CalcMoment(n,t,y)
%Fonction pour le calcul des moments
% n : Le nombre de moments que vous voulez calculer
% t : Vecteur temps
% y : Vecteur contenant la fonction f(t)
M=[];       
for i=0:n   
   M=[M (-1)^i/factorial(i)*TP2_IDP_Sy(t'.^i.*y,dx)];
end
end

