clear
clc
close
disp(' Identification des Processus : La Méthode des Moments généralisés')
disp(' Auteurs : NECHAT Ghiles & MEGHNOUDJ Houssem ');
disp('           ENP, Automatique 2ème année. '),

%initialisation
%     t_fin=10;
%     t=linspace(0,t_fin,1001);
%     a=5;                                  
%     s=tf('s');
t=0:0.01:50;
a=1;
s=tf('s');
K=2; r=1; n=4; T=3;
G=K*exp(-r*s)/(1+s*T)^n;
%système :
    %tu peux mettre un choix pour des exemple déjà présent comme exemple
    %G=2*exp(-2*s)/(1+2*s)^4;
    %G=(1+2*s)/(5+2*s+s^2);
    %G=4*exp(-0.2*s)/(1+0.4*s)^3;
    %G=3/(-4+1*s);
    %G=zpk([],[-1/4 -1/4 -1/4],4/64,'InputDelay',2);
    %G=input('Introduisez votre fonction de transfert :\n G = ');

%Entrée :
    for i=1:length(t)
    u(i,1)=min(exp(t(i)-5),exp(-t(i)+5))-exp(-5);
    end
    %u=[ones(100,1);zeros(1001-100,1)];
    %u=[zeros(1,1);ones(1000,1)];
    fprintf('\nChoisissez l''un des modèles suivants :\n');
    fprintf('\t FP Forme Polynomiale\n')
    fprintf('\t SR STREJC Avec Retard\n')
    choix=input('\nVotre choix : ','s');


%affichage premier
subplot(2,1,1)
grid on
box off
hold on
title('Entr{\''e}e du syst{\`e}me {\`a} identifier','interpret','latex','Fontsize', 13)
plot(t,u,'linewidth',1.1)
plot(t,exp(-a*t'),'r--','linewidth',.01)
ylabel('Input','interprete','latex','Fontsize', 12)
U=exp(-a*t').*u;
hLeg=legend({'$$u_e=f(t)$$','$$e^{-a.t}$$'},'Location','best','Fontsize', 12);
set(hLeg,'Interpreter','latex');

y=lsim(G,u,t);
subplot(2,1,2)
grid on
box off
hold on
title('Sortie du syst{\`e}me {\`a} identifier','interpret','latex','Fontsize', 13)
plot(t,y,'linewidth',1.1)
plot(t,exp(-a*t'),'r--','linewidth',.01)
xlabel('Temps (s)','interprete','latex','Fontsize', 12)
ylabel('Output','interprete','latex','Fontsize', 12)
Y=exp(-a*t').*y;
hLeg=legend({'$$y_s=g(t)$$','$$e^{-a.t}$$'},'Location','best','Fontsize', 12);
set(hLeg,'Interpreter','latex');
pause

clf
%affichage second
subplot(2,1,1)
grid on
box off
hold on
title({'Entr{\''e}e du syst{\`e}me {\''a} identifier','apr{\''e}s ajout des poids'},'interpret','latex','Fontsize', 13)
plot(t,U,'linewidth',1)
ylabel('Input','interprete','latex','Fontsize', 12)
hLeg=legend({'$$u_e=f(t).e^{-a.t}$$'},'Location','best','Fontsize', 12);
set(hLeg,'Interpreter','latex');

y=lsim(G,u,t);
subplot(2,1,2)
grid on
box off
hold on
title({'Sortie du syst{\`e}me {\''a} identifier','apr{\''e}s ajout des poids'},'interpret','latex','Fontsize', 13)
xlabel('Temps (s)','interprete','latex','Fontsize', 12)
ylabel('Output','interprete','latex','Fontsize', 12)
plot(t,Y,'linewidth',1)
hLeg=legend({'$$y_s=g(t).e^{-a.t}$$'},'Location','best','Fontsize', 12);
set(hLeg,'Interpreter','latex');
pause

%coté système
switch choix
    
    case 'FP'
        n=order(G);
        Me=fonctions.CalcMoment(2*n,t,U);
        Ms=fonctions.CalcMoment(2*n,t,Y);
        M =fonctions.CalcMomentGen(Me,Ms);
        D=fonctions.CalcParamPolynome(M,n,a);
        D=fliplr(D');
        Gm=tf(D(n+1:end),[1 D(1:n)]);
        round_a=round(D(1:n),2);
        round_b=round(D(n+1:end),2);
        den=[];
        num=['S^',num2str(n),' + '];
        for h=0:n-1
            if h==0
        den=[den,num2str(round_b(1+h)),'S^',num2str(n-h-1)];
        num=[num,num2str(round_a(1+h)),'S^',num2str(n-h-1)];
            elseif h==n-1
        den=[den,' + ',num2str(round_b(1+h))];
        num=[num,' + ',num2str(round_a(1+h))];
        else
        den=[den,' + ',num2str(round_b(1+h)),'S^',num2str(n-h-1)];
        num=[num,' + ',num2str(round_a(1+h)),'S^',num2str(n-h)];
        end
        end
        Val=['$G_m(s)=\frac{',den,'}{',num,'}$'];
        
        
    case 'SR'
        n=4; 
        Me=fonctions.CalcMoment(n,t,U);      %Moment de l'entrée
        Ms=fonctions.CalcMoment(n,t,Y);      %Moment de la sortie
        M =fonctions.CalcMomentGen(Me,Ms);   %Moment impulsionnel
        m =[M(1) M(2:end)./M(1)];
        
      %calcul des paramètres
        T_0=(2*m(2)^3 - 3*m(3)*m(2) + m(4))/(- 2*m(2)^2 + 2*m(3));
        T=T_0/(1-a*T_0);
        Tau=(m(2)^2*m(3) + m(4)*m(2) - 2*m(3)^2)/(2*m(2)^3 - 3*m(3)*m(2) + m(4));
        n=round((4*(-m(2)^2 + m(3))^3)/(2*m(2)^3 - 3*m(3)*m(2) + m(4))^2);
        K=m(1)*(1+a*T)^n/(exp(-a*Tau));
        Gm=K*exp(-s*Tau)/(1+s*T)^n
        Val=['$G_m(s)=\frac{',num2str(round(K,2)),'  e^{-',num2str(round(Tau,2)),'  t}}{(1+s',num2str(round(T,2)),')^',num2str(n),'}$'];
    
    otherwise
        choix=input('\nVeuillez introduire corréctement votre choix : ');
end


% Partie validation et comparaison
clf
choix=1;     
switch choix
    case 1  %réponse impulsionnel
ym=impulse(Gm,t);
y=impulse(G,t);
J=fonctions.IntSimpson((y-ym).^2,t(2));
y_label=['J = ',num2str(J)];
title({'R{\''e}ponse Impulsionnelle ( Syst{\`e}me + Mod{\`e}le )';Val},'interpret','latex','fontsize',14)
xlabel('Temps (s)','interprete','latex','fontsize',12)
ylabel({'Erreur Quadratique';y_label},'interprete','latex','fontsize',12)
hold on
plot(t,y)
plot(t,ym,'r--','linewidth',1.1)
hold off
grid on
box off
legend({'$$y(t)$$','$$y_m(t)$$'},'Interpreter','latex','Location','best','fontsize',12);

    case 2  %réponse indicielle
ym=step(Gm,t);
y=step(G,t);
J=fonctions.IntSimpson((y-ym).^2,t(2));
y_label=['J = ',num2str(J)];
title({'R{\''e}ponse indicielle ( Syst{\`e}me + Mod{\`e}le )',Val},'interpret','latex','fontsize',14)
xlabel('Temps (s)','interprete','latex','fontsize',12)
ylabel({'Erreur Quadratique';y_label},'interprete','latex','fontsize',12)
hold on
plot(t,y)
plot(t,ym,'r--','linewidth',1.1)
hold off
grid on
box off
legend({'$$y(t)$$','$$y_m(t)$$'},'Interpreter','latex','Location','best','fontsize',12);
end
