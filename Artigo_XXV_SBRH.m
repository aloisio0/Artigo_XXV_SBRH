clear;
clc;

%---------Variáveis de entrada---------
L=1500; % Comprimento da tubulação
d=1; % Diâmetro da tubulação
f=0.02; % Fator de atrito para o regime permanente. Coeficiente de Darcy-Weisbach
zr=100; % Cota do nível d'água no reservatório em relação ao eixo do tubo
g=9.8; % Aceleração da gravidade
cel=1000; % Velocidade da onda - Celeridade
velocidade=1; % Velocidade no Regime Permanente
tt=60; % Período total de simulação
N=15; % Número de segmentos na tubulação - Unidades computacionais
tc=3; % Período de fechamento - Tempo que leva para fechar a válvula
ts=0; % Hora do inínio do fechamento
se=0; % Porcentagem final de abertura
Em=3; % Coeficiente da válvula

%---------Cálculos iniciais---------
area=pi*d^2/4; % Área transversal da tubulação
q0=velocidade*area; % Vazão no Regime Permanente
dx=L/N;% Comprimento dos segmentos da tubulação - Distância entre os nós
dt=dx/cel;%Passo de tempo - Espaçamento entre nós do eixo temporal
tn=ceil(tt/dt)+1; %Total de Passos de tempo

%---------Matrizes de entrada---------
x=zeros(N+1,1);
T=zeros(tn,1);
ope=zeros(tn,1);
Q=zeros(tn,N+1);
H=zeros(tn,N+1);
dVdt=zeros(tn,N+1);
dVdx=zeros(tn,N);
HMAX=zeros(N+1,1);
HMIN=zeros(N+1,1);

%---------Regra de operação da válvula (tau)---------
for i=1:tn
    T(i)=dt*(i-1);
    aux=1-((i*dt)/tc)^Em;
    if (aux>=0)
        ope(i)=aux;
    end
    if (aux<0)
        ope(i)=0;
    end
end

%---------Vetores Q e H no Regime Permanente---------
for j=1:N+1
    x(j)=(j-1)*dx;
    Q(1,j)=q0; % Vazão do regime permanente
    H(1,j)=zr-(8*f*q0^2*x(j))/(g*pi^2*d^5);
end
H(:,1)=zr;
%dVdt(1)=zeros(N+1);
dVdx(1,:)=diff(Q(1,:)/area)/(L/N);

%---------Determinação do coeficiente de rugosidade D-W com base nas condições iniciais---------
tol=1;
v = 1.004e-6; %Viscosidade cinemática
inicial_Re=abs((q0/area)*d/v);
if f~=0 
    RH=max(10^(-1/1.8/sqrt(f))-6.9/inicial_Re,0);
else
    RH=0;
end

%---------Vetores Q e H na Regime Transitório---------
% Condições de contorno
for i=2:tn
    % Montante (Reservatório):

    % Atualização do Fator de atrito e cálculo da perda carga quase estacionária
    Re=abs(Q(i-1,2)/area*d)/v;
    if (Re<tol) 
        hfs=0;
    else
        a=-1.8*log10(6.9/Re+RH);
        f=(1./a)^2;
        hfs=(f/2./g/d*(Q(i-1,2)/area)*abs(Q(i-1,2)/area));
    end

    %Cálculo do coeficiente de decaimento de cisalhamento de Vardy (C)
    if (Re< 2000) %Escoamento laminar
        C = 4.76e-3;
    else %Escoamento de transição ou turbulento
        C = 7.41/Re^(log10(14.3/Re^0.05));
    end

    % Cálculo do coeficiente de atrito de Brunone (ku)
    ku=2*sqrt(C);
    hfu=ku/g*(dVdt(i-1,2)+cel*sign(Q(i-1,2)/area)*abs(dVdx(i-1,1)));
    hf=hfu+hfs;

    % A carga hidráulica no Reservatório é constante e igual ao nível da água. 
    Q(i,1)=(Q(i-1,2)/area+g/cel*(H(i,1)-H(i-1,2))-g*dt*hf)*area; % Cálculo da vazão

    % Jusante (Válvula):

    % Atualização do Fator de atrito e cálculo da perda carga quase estacionária
    Re=abs(Q(i-1,N)/area*d)/v;
    if (Re<tol)
        hfs=0;
    else
        a=-1.8*log10(6.9/Re+RH);
        f=(1./a)^2;
        hfs=(f/2./g/d*(Q(i-1,N)/area)*abs(Q(i-1,N)/area));
    end

    %Cálculo do coeficiente de decaimento de cisalhamento de Vardy (C)
    if (Re< 2000) %Escoamento laminar
        C = 4.76e-3;
    else %Escoamento de transição ou turbulento
        C = 7.41/Re^(log10(14.3/Re^0.05));
    end

    % Cálculo do coeficiente de atrito de Brunone (ku)
    ku=2*sqrt(C);
    hfu=ku/g*(dVdt(i-1,N)+cel*sign(Q(i-1,N)/area)*abs(dVdx(i-1,N)));
    hf=hfu+hfs;
    
    Q(i,N+1)=q0*ope(i-1); % Cálculo da vazão
    c1=Q(i-1,N)/area+g/cel*H(i-1,N)-g*dt*hf;
    H(i,N+1)=(c1-Q(i,N+1)/area)*cel/g; % Cálculo da carga de pressão
    
    % Nós internos:
    for j=2:N

        %Linha característica C+
        % Atualização do Fator de atrito e cálculo da perda carga quase estacionária
        Re=abs(Q(i-1,j-1)/area*d)/v;
        if (Re<tol) 
            hfs=0;
        else
            a=-1.8*log10(6.9/Re+RH);
            f=(1./a)^2;
            hfs=(f/2./g/d*(Q(i-1,j-1)/area)*abs(Q(i-1,j-1)/area));
        end

        %Cálculo do coeficiente de decaimento de cisalhamento de Vardy (C)
        if (Re< 2000) %Escoamento laminar
            C = 4.76e-3;
        else %Escoamento de transição ou turbulento
            C = 7.41/Re^(log10(14.3/Re^0.05));
        end

        % Cálculo do coeficiente de atrito de Brunone (ku)
        ku=2*sqrt(C);
        hfu=ku/g*(dVdt(i-1,j-1)+cel*sign(Q(i-1,j-1)/area)*abs(dVdx(i-1,j-1)));
        hf1=hfu+hfs;

        c1=(Q(i-1,j-1)/area)+g/cel*H(i-1,j-1)-g*dt*hf1;
        
        %Linha característica C-
        % Atualização do Fator de atrito
        Re=abs(Q(i-1,j+1)/area*d)/v;
        if (Re<tol)
            hfs=0;
        else
            a=-1.8*log10(6.9/Re+RH);
            f=(1./a)^2;
            hfs=(f/2./g/d*(Q(i-1,j+1)/area)*abs(Q(i-1,j+1)/area));
        end

        %Cálculo do coeficiente de decaimento de cisalhamento de Vardy (C)
        if (Re< 2000) %Escoamento laminar
            C = 4.76e-3;
        else %Escoamento de transição ou turbulento
            C = 7.41/Re^(log10(14.3/Re^0.05));
        end

        % Cálculo do coeficiente de atrito de Brunone (ku)
        ku=2*sqrt(C);
        hfu=ku/g*(dVdt(i-1,j+1)+cel*sign(Q(i-1,j+1)/area)*abs(dVdx(i-1,j)));
        hf2=hfu+hfs;
        c2=-(Q(i-1,j+1)/area)+g/cel*H(i-1,j+1)+g*dt*hf2;

        H(i,j)=(c1+c2)/(2*g/cel); % Cálculo da carga de pressão
        Q(i,j)=(+g/cel*H(i,j)-c2)*area; % Cálculo da vazão
    end

    dVdt(i,:)=((Q(i,:)-Q(i-1,:))/area)/dt; % Aceleração instantânea convectiva
    dVdx(i,:)=diff(Q(i-1,:)/area)/(L/N); % Aceleração local
end

% Envoltórias Máximas e Mínimas de Pressão
for j=1:N+1
    HMAX(j)=max(H(:,j));
    HMIN(j)=min(H(:,j));
end

% Variação da Vazão e Carga de Pressão no Nó da Válvula durante o regime transiente
HN1=zeros(tn);
QN1=zeros(tn);
for j=1:tn
    HN1(j)=H(j,N+1);
    QN1(j)=Q(j,N+1);
end

% Variação da Vazão e Carga de Pressão no Nó do Reservatório durante o regime transiente
HN0=zeros(tn);
QN0=zeros(tn);
for j=1:tn
    HN0(j)=H(j,1); 
    QN0(j)=Q(j,1);
end

% Animações
result=0;
while result~=4
    %RESULTADOS:
    result=menu('Vizualização dos resultados:','Animação da Carga de Pressão no Nó da Válvula', 'Animação da Vazão no Nó do Reservatório','Animação das Envoltórias de Pressão', 'Sair');
    if result==1
        figure(1)
        curve=animatedline('Colo','b','LineWidth',2);
        set(gca,'Xlim',[0 tt], 'Ylim',[0 220])
        xlabel('Tempo (s)')
        ylabel('Carga de Pressão (m.c.a.)')
        %legend('Carga na Válvula')
        for i=1:length(T)
            addpoints(curve,T(i),HN1(i));
            title(['Carga na Válvula no tempo t = ',num2str(T(i)),'segundos'])
            drawnow
            pause(0.05);
        end
    else
        if result==2
        figure(1)
        curve=animatedline('Colo','b','LineWidth',2);
        set(gca,'Xlim',[0 tt], 'Ylim',[-1 1])
        xlabel('Tempo (s)')
        ylabel('Vazão (m³/s)')
        %legend('Carga na Válvula')
        for i=1:length(T)
            addpoints(curve,T(i),QN0(i));
            title(['Vazão no Reservatório no tempo t = ',num2str(T(i)),'segundos'])
            drawnow
            pause(0.05);
        end
    else 
        if result==3
            for i=1:length(T)
                figure(2)
                plot(x,H(1,:),'b',x,HMAX,'r',x,HMIN,'g',x,H(i,:),'k','LineWidth',2)
                xlabel('Distância (m)')
                ylabel('Carga de Pressão (m.c.a.)')
                title(['Carga no tempo t = ',num2str(T(i)),'segundos'])
                legend('Regime Permanente','Envoltória Máxima', 'Envoltória Mínima','Envoltórias de Pressão','Location','northwest')
                axis([0 L 0 220])
                grid on
                pause(0.01);
            end
        end
        end
    end
end