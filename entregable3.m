%% Especificación del problema

% Índice de refracción del medio a la izquierda del multicapa.
n0=1;
% Índices de refracción de cada una de las capas. [n1 n2 n3 n4 n5 n6].
n=[2.32 1.38 2.32 1.38 2.32 1.38];
% Índice de refracción del medio a la derecha del multicapa.
nF=1.52;
% Grosor (en metros) de cada una de las capas. [d1 d2 d3 d4 d5 d6].
d=10^-6*[0.086 0.145 0.086 0.145 0.086 0.145];

%% Código del programa

% Número de capas que conforman el multicapa.
N=length(d);
% Vector que contiene las longitudes de onda (en metros) que se estudian.
lambda=10^-6*(0.6:0.00001:1.2);
% Valor exacto de la impedancia de onda en el vacío.
etaVacuum=119.9169832*pi;
% Valor de la impedancia de onda del medio a la izquierda del multicapa.
eta0=etaVacuum/n0;
% Vector que contiene las impedancias de onda de cada una de las capas.
eta=etaVacuum./n;
% Valor de la impedancia de onda del medio a la derecha del multicapa.
etaF=etaVacuum/nF;

% Matriz auxiliar para aprovechar la capacidad de vectorización.
aux=zeros(length(lambda),N);
for ii=1:length(lambda)
    aux(ii,:)=2*pi./lambda(ii)*n.*d;
end

%% Método "fórmula del coeficiente de reflexión" (CDR)

% Método recursivo para encontrar los coeficientes de reflexión.

% Se crea un vector que contiene los coeficientes rho.
rho=zeros(1,N);
rho(N)=(etaF-eta(N))/(etaF+eta(N));
for ii=1:N-1
    rho(N-ii)=(eta(N-ii+1)-eta(N-ii))/(eta(N-ii+1)+eta(N-ii));
end
rho0=(eta(1)-eta0)/(eta(1)+eta0);

% Se calculan los coeficientes de reflexión recurrentemente.
CDRgamma=zeros(length(lambda),N);
% gamma_(N+1) es cero porque a la derecha de n_F no hay más capas.
CDRgamma(:,N)=rho(N);
for ii=1:N-1
    CDRgamma(:,N-ii)=(rho(N-ii)+CDRgamma(:,N-ii+1).*...
        exp(-2*1i*aux(:,N-ii+1)))./(1+rho(N-ii)*CDRgamma(:,N-ii+1).*...
        exp(-2*1i*aux(:,N-ii+1)));
end

% Se calcula el coeficiente de reflexión final.
CDRgammaTOTAL=(rho0+CDRgamma(:,1).*exp(-2*1i*aux(:,1)))./...
        (1+rho0*CDRgamma(:,1).*exp(-2*1i*aux(:,1)));

% A partir de este, se calculan las fracciones de potencia reflejada y
% potencia transmitida.
CDRfracPotReflejada=(abs(CDRgammaTOTAL)).^2;
CDRfracPotTransmitida=(1-(abs(CDRgammaTOTAL)).^2);

% Se representa el resultado en un gráfico. Mucho código poco importante.
CDRa0=find(CDRfracPotReflejada==max(CDRfracPotReflejada));
CDRa1=find(CDRfracPotReflejada>CDRfracPotTransmitida);
figure(1)
subplot(2,1,1)
plot(10^6*lambda,CDRfracPotReflejada,'r',10^6*lambda,...
    CDRfracPotTransmitida,'g',10^6*lambda(CDRa0),...
    max(CDRfracPotReflejada),'xk',10^6*lambda(min(CDRa1)),...
    CDRfracPotReflejada(min(CDRa1)),'xk',10^6*lambda(max(CDRa1)),...
    CDRfracPotReflejada(max(CDRa1)),'xk')
grid;xlabel('\lambda (\mum)');axis([10^6*lambda(1) 10^6*lambda(end) 0 1])
legend('Fracción de potencia reflejada','Fracción de potencia transmitida')
text(10^6*lambda(CDRa0),max(CDRfracPotReflejada),['Máximo del ' ...
    num2str(100*max(CDRfracPotReflejada)) '% para \lambda = ' ...
    num2str(10^6*lambda(CDRa0)) ' \mum'])
text(10^6*lambda(min(CDRa1)),CDRfracPotReflejada(min(CDRa1)),...
    ['50% para \lambda = ' num2str(10^6*lambda(min(CDRa1))) ' \mum'])
text(10^6*lambda(max(CDRa1)),CDRfracPotReflejada(max(CDRa1)),...
    ['50% para \lambda = ' num2str(10^6*lambda(max(CDRa1))) ' \mum'])
title('Método "fórmula del coeficiente de reflexión"')

%% Método "fórmula de la impedancia de entrada" (IDE)

% Método recursivo para encontrar las impedancias de entrada.

% Se van calculando impedancias de entrada hasta simplificarlo a una sola.
Zin=zeros(length(lambda),N);
Zin(:,N)=eta(N)*(etaF+1i*eta(N)*tan(aux(:,N)))./(eta(N)+1i*etaF*...
    tan(aux(:,N)));
for ii=1:N-1
    Zin(:,N-ii)=eta(N-ii)*(Zin(:,N-ii+1)+1i*eta(N-ii)*...
        tan(aux(:,N-ii)))./(eta(N-ii)+1i*Zin(:,N-ii+1).*tan(aux(:,N-ii)));
end

% Se calcula el coeficiente de reflexión final mediante la impedancia de
% entrada final.
IDEgammaTOTAL=(Zin(:,1)-eta0)./(Zin(:,1)+eta0);

% A partir de este, se calculan las fracciones de potencia reflejada y
% potencia transmitida.
IDEfracPotReflejada=(abs(IDEgammaTOTAL)).^2;
IDEfracPotTransmitida=(1-(abs(IDEgammaTOTAL)).^2);

% Se representa el resultado en un gráfico. Mucho código poco importante.
IDEa0=find(IDEfracPotReflejada==max(IDEfracPotReflejada));
IDEa1=find(IDEfracPotReflejada>IDEfracPotTransmitida);
figure(1)
subplot(2,1,2)
plot(10^6*lambda,IDEfracPotReflejada,'r',10^6*lambda,...
    IDEfracPotTransmitida,'g',10^6*lambda(IDEa0),...
    max(IDEfracPotReflejada),'xk',10^6*lambda(min(IDEa1)),...
    IDEfracPotReflejada(min(IDEa1)),'xk',10^6*lambda(max(IDEa1)),...
    IDEfracPotReflejada(max(IDEa1)),'xk')
grid;xlabel('\lambda (\mum)');axis([10^6*lambda(1) 10^6*lambda(end) 0 1])
legend('Fracción de potencia reflejada','Fracción de potencia transmitida')
text(10^6*lambda(IDEa0),max(IDEfracPotReflejada),['Máximo del ' ...
    num2str(100*max(IDEfracPotReflejada)) '% para \lambda = ' ...
    num2str(10^6*lambda(IDEa0)) ' \mum'])
text(10^6*lambda(min(IDEa1)),IDEfracPotReflejada(min(IDEa1)),...
    ['50% para \lambda = ' num2str(10^6*lambda(min(IDEa1))) ' \mum'])
text(10^6*lambda(max(IDEa1)),IDEfracPotReflejada(max(IDEa1)),...
    ['50% para \lambda = ' num2str(10^6*lambda(max(IDEa1))) ' \mum'])
title('Método "fórmula de la impedancia de entrada"')

%% Conclusiones

% ¿Para qué longitudes de onda se comporta como un espejo?
% Si la definición de espejo es que refleja más que absorbe...
% Entonces este multicapa se comporta como un espejo de 643.78 nm a 1053.4 nm.
% Nota: entre 400 nm y 750 nm se encuentra el espectro de luz visible.