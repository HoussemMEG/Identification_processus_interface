%Contient toutes les fonctions utilisées
classdef fonctions
    methods(Static)

        function [M] = CalcMoment(n,t,y)
        %Fonction pour le calcul des moments
        % n : Le nombre de moments que vous voulez calculer
        % t : Vecteur temps (uni colonne)
        % y : Vecteur contenant la fonction f(t) (uni colonne)
        % M : Vecteur contenant les moments calculés
        if nargin~=3
            error('Vous n''avez pas introduit les bons paramètres');
        else
        M=[];
        for i=0:n-1   
           M=[M fonctions.IntSimpson(t'.^i.*y,t(2))];
        end
        end
        enddan

        function [Sy] = IntSimpson(Y,DX)
        % Fonction qui retourne l'intégrale de Simpson,
        % Y : Les Abcisses & DX : Le pas d'intégration.
        if not(mod(length(Y),2)) 
            error('Le Nombre de points doit être impaire.');
        end
        Sy = DX*( sum(Y(1:end-1)) + sum(Y(2:end)) + 2*sum(Y(2:2:end-1)) )/3 ;
        end

        function [M] = CalcMomentGen(Me,Ms)
        % Fonction qui calcule les moments impulstionnels
        % A partir des moments de l'entrée et de la sortie
        % Me : Moments de l'entrée
        % Ms : Moments de la sortie
        n=length(Me);
        M=Ms(1)/Me(1);
        for k=1:n-1
            S=0;
            for i=0:k-1
                S=S+nchoosek(k,i)*M(i+1)*Me(k-i+1);
            end
        M=[M (Ms(k+1)-S)/Me(1)];
        end
        end
        
        function [S] = CalcParamPolynome(M,n,a_val)
            syms a
            V=[];
            G=zeros(2*n,n+1);
            for k=0:n 
            V=[V a^k];
            end
            for k=1:n
            V=[V;diff(V(end,:),a)];
            end
            V=subs(V,a,a_val);

            for i=1:2*n
                for k=1:n+1 
                 if i+k-1<2*n+1 %il y a le shift de M(i)=M(i+1)
                 G(i+k-1,:)=G(k+i-1,:)+nchoosek(i+k-2,k-1)*(-1)^(i-1)*M(i)*V(k,:);
                 end 
                end
            end
            G=double([V(:,1:end-1) G]);
            A=[G(:,1:n) -1.*G(:,n+1:end-1)];
            b=G(:,end);
            S=inv(A)*b;
        end
    end
    end
end